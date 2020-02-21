from collections import namedtuple, defaultdict, Counter
import numpy as np
import pandas as pd
import yaml

try:
    from .log import log_
except (ImportError, ModuleNotFoundError):
    from log import log_


Probs = namedtuple("Probs", ("O", "N", "P_cont", "alpha", "beta", "lib"))
Probs2 = namedtuple(
    "Probs", ("O", "N", "P_cont", "alpha", "beta", "lib", "alpha_hap", "beta_hap")
)
Pars = namedtuple(
    "Pars", ("alpha0", "trans_mat", "cont", "error", "F", "tau", "gamma_names", "sex")
)  # for returning from varios functions
ParsHD = namedtuple(
    "ParsHD",
    (
        "alpha0",
        "alpha0_hap",
        "trans",
        "trans_hap",
        "cont",
        "error",
        "F",
        "tau",
        "gamma_names",
        "sex",
    ),
)  # for returning from varios functions


class _IX:
    """class to be filled with various indices
    """
    def __init__(self):
        pass


def data2probs(
    df,
    IX,
    state_ids,
    cont_id=None,
    prior=None,
    cont_prior=(1e-8, 1e-8),
    ancestral=None,
    ancestral_prior = 0
):
    """create data structure that holds the reference genetic data

    creates an object of type `Probs` with the following entries:
    O : array[n_obs]: the number of alternative reads
    N : array[n_obs]: the total number of reads
    P_cont : array[n_obs]: the contaminant allele frequency
    lib[n_obs] : the library /read group of the observation
    alpha[n_snps, n_states] : the reference allele beta-prior
    beta[n_snps, n_states] : the alt allele beta-prior


    input:

    df: merged reference and SNP data. has columns tref, talt with the read
    counts at each SNP, and "X_alt, X_ref" for each source pop
    IX: index object, with number of snps, number of reads, etc.
    state_ids: the references to keep
    prior: None for empirical bayes prior, otherwise prior to be added
    ancestral: ancestray allele
    """

    alt_ix = ["%s_alt" % s for s in state_ids]
    ref_ix = ["%s_ref" % s for s in state_ids]
    snp_ix_states = set(alt_ix + ref_ix)

    if cont_id is not None:
        cont = "%s_alt" % cont_id, "%s_ref" % cont_id
        snp_ix_states.update(cont)
    if ancestral is not None:
        anc = "%s_alt" % ancestral, "%s_ref" % ancestral
        snp_ix_states.update(anc)

    snp_df = df[list(snp_ix_states)]
    snp_df = snp_df[~snp_df.index.get_level_values('snp_id').duplicated()]
    #snp_df = df[list(snp_ix_states)].groupby(df.index.names).first()
    n_snps = len(snp_df.index.get_level_values('snp_id'))
    n_states = len(state_ids)


    if prior is None:  # empirical bayes, estimate from data
        alt_prior = np.empty((n_snps, n_states))
        ref_prior = np.empty((n_snps, n_states))
        if cont_id is not None:
            ca, cb = empirical_bayes_prior(snp_df[cont[0]], snp_df[cont[1]])

        if ancestral is None:
            for i, (a, b, s) in enumerate(zip(alt_ix, ref_ix, state_ids)):
                pa, pb = empirical_bayes_prior(snp_df[a], snp_df[b])
                log_.info("[%s]EB prior [a=%.4f, b=%.4f]: " % (s, pa, pb))
                alt_prior[:, i] = snp_df[a] + pa
                ref_prior[:, i] = snp_df[b] + pb
        else:
            anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"

            #set up vectors stating which allele is ancestral
            ref_is_anc = (snp_df[anc_ref] > 0) & (snp_df[anc_alt] == 0)
            alt_is_anc = (snp_df[anc_alt] > 0) & (snp_df[anc_ref] == 0)
            ref_is_der, alt_is_der = alt_is_anc, ref_is_anc
            anc_is_unknown = (1 - alt_is_anc) * (1 - ref_is_anc) == 1

            for i, (alt_col, ref_col, s) in enumerate(zip(alt_ix, ref_ix, state_ids)):

                #1. set up base entries based on observed counts
                alt_prior[:, i] = snp_df[alt_col]
                ref_prior[:, i] = snp_df[ref_col]

                #2. where anc is unknown, add symmetric prior estimated from data
                pa, pb = empirical_bayes_prior(snp_df[alt_col], snp_df[ref_col])
                log_.info("[%s]EB prior0 [anc=%.4f, der=%.4f]: " % (s, pa, pb))
                alt_prior[anc_is_unknown, i] += pa
                ref_prior[anc_is_unknown, i] += pb

                #3. where anc is known, create indices
                m_anc = pd.concat((ref_is_anc, alt_is_anc), 1)
                m_der = pd.concat((ref_is_der, alt_is_der), 1)
                ANC = np.array(snp_df[[ref_col, alt_col]])[m_anc]
                DER = np.array(snp_df[[ref_col, alt_col]])[m_der]

                pder, panc = empirical_bayes_prior(DER, ANC, known_anc=True)
                panc += ancestral_prior
                log_.info("[%s]EB prior1 [anc=%.4f, der=%.4f]: " % (s, panc, pder))
                alt_prior[alt_is_anc, i] += panc 
                alt_prior[alt_is_der, i] += pder
                ref_prior[ref_is_anc, i] += panc
                ref_prior[ref_is_der, i] += pder

        P = Probs2(
            O=np.array(df.talt.values, np.int8),
            N=np.array(df.tref.values + df.talt.values, np.int8),
            P_cont=np.zeros_like(df.talt.values)
            if cont_id is None
            else np.array(
                (df[cont[0]].values + ca) / (df[cont[0]].values + df[cont[1]].values + ca + cb)
            ),
            alpha=alt_prior[IX.diploid_snps],
            beta=ref_prior[IX.diploid_snps],
            alpha_hap=alt_prior[IX.haploid_snps],
            beta_hap=ref_prior[IX.haploid_snps],
            lib=np.array(df.lib),
        )
        return P
    else:

        """ancestral allele contribution to prior
        the ancestral allele adds one pseudocount to the data
        """
        if ancestral is None:
            prior_anc_alt, prior_anc_ref = np.zeros(1), np.zeros(1)
        else:
            anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"
            prior_anc_alt = snp_df[anc_alt] * ancestral_prior
            prior_anc_ref = snp_df[anc_ref] * ancestral_prior

        cont = "%s_alt" % cont_id, "%s_ref" % cont_id
        ca, cb = cont_prior

        alt_prior = snp_df[alt_ix].to_numpy() + prior_anc_alt[:, np.newaxis] + prior
        ref_prior = snp_df[ref_ix].to_numpy() + prior_anc_ref[:, np.newaxis] + prior
        #breakpoint()
        P = Probs2(
            O=np.array(df.talt.values, np.int8),
            N=np.array(df.tref.values + df.talt.values, np.int8),
            P_cont=0.
            if cont_id is None
            else np.array(
                (df[cont[0]] + ca) / (df[cont[0]] + df[cont[1]] + ca + cb)
            ),
            alpha=alt_prior[IX.diploid_snps],
            beta=ref_prior[IX.diploid_snps],
            alpha_hap=alt_prior[IX.haploid_snps],
            beta_hap=ref_prior[IX.haploid_snps],
            lib=np.array(df.lib),
        )
        return P


def bins_from_bed(df, bin_size, sex=None):
    """create a bunch of auxillary data frames for binning

    - bins: columns are chrom_id, chrom, bin_pos, bin_id, map
    - IX: container storing all indices, will need to be cleaned later on
        currently contains:
        - IX.SNP2BIN [n_snps]: array giving bin for each snp
        - IX.OBS2SNP [n_obs]: array giving snp for each obs
        - IX.OBS2BIN [n_obs]: array giving bin for each obs
        - IX.bin_sizes [n_chroms]: number of bins per chromosome
        - IX.RG2OBS [n_libs] : [list] dict giving obs for each readgroup
        - IX.libs : names of all libraries
    """
    IX = _IX()
    IX.libs = np.unique(df.lib)

    obsix = df.index.to_frame(index=False)
    snp = obsix.drop_duplicates()

    chroms = pd.unique(snp.chrom)
    n_snps = len(snp.snp_id.unique())

    haplo_chroms, diplo_chroms = [], []
    if sex is None:
        diplo_chroms = chroms
    else:
        for c in chroms:
            if c[0] in "YyWw":
                haplo_chroms.append(c)
            elif c[0] in "Xx" and sex == "m":
                haplo_chroms.append(c)
            elif c[0] in "Zz" and sex == "f":
                haplo_chroms.append(c)
            elif c.startswith("hap"):
                haplo_chroms.append(c)
            else:
                diplo_chroms.append(c)


    IX.SNP2BIN = np.empty((n_snps), int)
    IX.OBS2SNP = obsix.snp_id.values

    bin_loc = []
    bin0 = 0

    dtype_bin = np.dtype(
        [("chrom", "U2"), ("map", float), ("pos", int), ("id", int), ("haploid", bool)]
    )

    IX.bin_sizes = []
    IX.diploid_snps, IX.haploid_snps = [], []

    for i, chrom in enumerate(chroms):
        map_ = snp.map[snp.chrom == chrom]
        pos = snp.pos[snp.chrom == chrom]

        chrom_start = float(np.floor(map_.head(1) / bin_size) * bin_size)
        chrom_end = float(np.ceil(map_.tail(1) / bin_size) * bin_size)
        chrom_is_hap = chrom in haplo_chroms

        # create bins
        bins = np.arange(chrom_start, chrom_end, bin_size)
        log_.debug("binning chrom %s: %d bins", chrom, len(bins))

        IX.bin_sizes.append(len(bins))
        bin_ids = range(bin0, bin0 + len(bins))
        _bin = np.empty_like(bins, dtype_bin)
        _bin["chrom"] = chrom
        _bin["pos"] = np.interp(bins, map_, pos)
        _bin["id"] = bin_ids
        _bin["map"] = bins
        _bin["haploid"] = chrom_is_hap
        bin_loc.append(_bin)



        # put SNPs in bins
        snp_ids = snp.snp_id[snp.chrom == chrom]
        dig_snp = np.digitize(snp[snp.chrom == chrom].map, bins, right=False) - 1
        IX.SNP2BIN[snp_ids] = dig_snp + bin0


        if chrom_is_hap:
            IX.haploid_snps.extend(snp_ids)
        else:
            IX.diploid_snps.extend(snp_ids)

        bin0 += len(bins)


    bins = np.hstack(bin_loc)

    # for now, assume data is ordered such that diploid chroms come before haploid ones
    assert (
        len(IX.haploid_snps) == 0
        or len(IX.diploid_snps) == 0
        or min(IX.haploid_snps) > max(IX.diploid_snps)
    )

    if len(IX.haploid_snps) > 0:
        IX.haploid_snps = slice(min(IX.haploid_snps), max(IX.haploid_snps) + 1)
    else:
        IX.haploid_snps = slice(0, 0)
    if len(IX.diploid_snps) > 0:
        IX.diploid_snps = slice(min(IX.diploid_snps), max(IX.diploid_snps) + 1)
    else:
        IX.diploid_snps = slice(0, 0)

    IX.RG2OBS = dict((l, np.where(df.lib == l)[0]) for l in IX.libs)
    IX.OBS2BIN = IX.SNP2BIN[IX.OBS2SNP]

    IX.n_chroms = len(chroms)
    IX.n_bins = len(bins)
    IX.n_snps = len(IX.SNP2BIN)
    IX.n_obs = len(IX.OBS2SNP)
    IX.n_reads = np.sum(df.tref + df.talt)

    IX.chroms = chroms
    IX.haplo_chroms = haplo_chroms
    IX.diplo_chroms = diplo_chroms

    log_.debug("done creating bins")
    return bins, IX


def init_pars(
    state_ids,
    sex=None,
    F0=0.001,
    tau0=1,
    e0=1e-2,
    c0=1e-2,
    est_inbreeding=False,
    init_guess=None,
    do_hap=True,
    **kwargs
):
    """initialize parameters

    returns a pars object
    """
    homo = [s for s in state_ids]
    het = []
    hap = ["h%s" % s for s in homo]

    for i, s in enumerate(state_ids):
        for s2 in state_ids[i + 1 :]:
            het.append(s + s2)
    gamma_names = homo + het
    if est_inbreeding:
        gamma_names.extend(hap)

    n_states = len(gamma_names)
    n_homo = len(homo)
    n_het = len(het)
    n_hap = len(hap)

    alpha0 = np.array([1 / n_states] * n_states)
    alpha0_hap = np.array([1 / n_hap] * n_hap)

    trans_mat = np.zeros((n_states, n_states)) + 2e-2
    trans_mat_hap = np.zeros((n_hap, n_hap)) + 2e-2

    np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
    np.fill_diagonal(trans_mat_hap, 1 - (n_hap - 1) * 2e-2)
    cont = defaultdict(lambda: c0)
    error = defaultdict(lambda: e0)

    if init_guess is not None:
        # guess = [i for i, n in enumerate(gamma_names) if init_guess in n]
        guess = [i for i, n in enumerate(gamma_names) if n in init_guess]
        log_.info("starting with guess %s " % guess)
        trans_mat[:, guess] = trans_mat[:, guess] + 1
        trans_mat /= np.sum(trans_mat, 1)[:, np.newaxis]

    try:
        if len(F0) == n_homo:
            F = F0
        elif len(F0) == 1:
            F = F0 * n_homo
        else:
            F = [F0]
    except TypeError:
        F = [F0] * n_homo
    try:
        if len(tau0) == n_homo:
            tau = tau0
        elif len(F0) == 1:
            tau = tau0 * n_homo
        else:
            tau = [tau0]
    except TypeError:
        tau = [tau0] * n_homo
    if do_hap:
        return ParsHD(
            alpha0,
            alpha0_hap,
            trans_mat,
            trans_mat_hap,
            cont,
            error,
            F,
            tau,
            gamma_names,
            sex=sex,
        )
    else:
        return Pars(alpha0, trans_mat, cont, error, F, tau, gamma_names, sex=sex)


def posterior_table(pg, Z, IX, est_inbreeding=False):
    freq = np.array([0, 1, 2, 0, 1]) if est_inbreeding else np.arange(3)
    PG = np.sum(Z[IX.SNP2BIN][:, :, np.newaxis] * pg, 1)  # genotype probs
    mu = np.sum(PG * freq, 1)[:, np.newaxis] / 2
    PG = np.log10(PG + 1e-40)
    PG = np.minimum(0.0, PG)
    return pd.DataFrame(np.hstack((PG, mu)), columns=["G0", "G1", "G2", "p"])


def empirical_bayes_prior(der, anc, known_anc=False):
    """using beta-binomial plug-in estimator
    """

    n = anc + der
    f = np.nanmean(der / n) if known_anc else 0.5

    H = f * (1.0 - f)  # alt formulation
    if H == 0.0:
        return 1e-6, 1e-6

    V = np.nanvar(der / n) if known_anc else np.nanvar(np.hstack((der / n, anc / n)))
    
    ab = (H - V) / (V - H / np.nanmean(n))
    if np.nanmean(n) < ab:
        return 1e-6, 1e-6
    pa = max((f * ab, 1e-5))
    pb = max(((1 - f) * ab, 1e-5))
    return pa, pb


def guess_sex(ref, data, sex_ratio_threshold=0.75):
    """
        guessing the sex of individuals by comparing heterogametic chromosomes.
        By convention, all chromosomes are assumed to be diploid unless they start
        with an `X` or `Z`
    """
    ref["heterogametic"] = [v[0] in "XZxz" for v in ref.index.get_level_values('chrom')]
    data["heterogametic"] = [v[0] in "XZxz" for v in data.index.get_level_values('chrom')]

    n_sites = ref.groupby(ref.heterogametic).apply(lambda df: len(df))
    n_reads = data.groupby(data.heterogametic).apply(lambda df: np.sum(df.tref + df.talt))
    cov = n_reads / n_sites

    del data["heterogametic"]
    del ref['heterogametic']

    #no heteogametic data
    if True not in cov:
        return 'f'


    if cov[True] / cov[False] < sex_ratio_threshold:
        sex = "m"
        log_.info("guessing sex is male, X/A = %.4f/%.4f" % (cov[True], cov[False]))
    else:
        sex = "f"
        log_.info("guessing sex is female, X/A = %.4f/%.4f" % (cov[True], cov[False]))
    return sex


def scale_mat(M):
    """scale a matrix of probabilities such that it's highest value is one

    modifies M and returns log(scaling)
    """
    scaling = np.max(M, 1)[:, np.newaxis]
    M /= scaling
    assert np.allclose(np.max(M, 1), 1)
    log_scaling = np.sum(np.log(scaling))
    return log_scaling


def scale_mat3d(M):
    """scale a matrix of probabilities such that it's highest value is one

    modifies M and returns log(scaling)
    """
    scaling = np.max(M, (1, 2))[:, np.newaxis, np.newaxis]
    M /= scaling
    assert np.allclose(np.max(M, (1, 2)), 1)
    log_scaling = np.sum(np.log(scaling))
    return log_scaling


def parse_state_string(states, ext=[''], operator='+', state_file=None):
    """parse shortcut parse strings

    the basic way states are defined is as 
        --states AFR NEA DEN

    this assmues the columns AFR, NEA and DEN are known and present in the reference
    We can rename and merge stuff here using the following syntax

    -- states AFR=YRI+SAN NEA=VIN DEN=Denisova2

    return rename dict for reference
    
    """
    ext2 = ['_ref', '_alt']
    d1 = [s.split("=") for s in states if s is not None]
    d2 = [(s if len(s)>1 else (s[0],s[0])) for s in d1] 
    state_dict = dict( (f'{k}{ext_}', f'{i}{ext_}') for i, j in d2
              for k in j.split(operator) 
              for ext_ in ext)
    if state_file is not None:
        """logic here is:
            - yaml file contains some groupings (i.e. assign all
            Yorubans to YRI, all San to SAN
            - state_dict contains further groupings / renaming / filtering
                i.e. AFR=YRI+SAN
            - stuff not in state_dict will not be required

            therefore:
            1. load yaml
            2. get all required target pops from state_dict
            3. add all expansions replacements
        """
        Y = yaml.load(open(state_file), Loader=yaml.BaseLoader)
        #D = dict((v_, k) for (k,v) in Y.items() for v_ in v)
        s2 = dict()

        for state, label in state_dict.items():
            if state in Y:
                for ind in Y[state]:
                    s2[ind] = label
            else:
                s2[state] = label
        state_dict = s2
    return state_dict


