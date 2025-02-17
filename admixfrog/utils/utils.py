import logging
from collections import namedtuple, defaultdict, Counter
import numpy as np
import pandas as pd
import yaml
from ..slug.classes import (
    SlugData,
    SlugParsSquare,
    SlugReads,
)
from numba import njit
from scipy.linalg import expm


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
    ancestral_prior=0,
    doalphabeta=True,
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
    doalphabeta: for e.g. SFS mode, alpha and beta don't need to be calcuated
    """

    alt_ix = ["%s_alt" % s for s in state_ids]
    ref_ix = ["%s_ref" % s for s in state_ids]
    snp_ix_states = set(alt_ix + ref_ix)
    ca, cb = cont_prior

    if cont_id is not None:
        cont_ref, cont_alt = "%s_alt" % cont_id, "%s_ref" % cont_id
        snp_ix_states.update([cont_ref, cont_alt])
    if ancestral is not None:
        anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"
        snp_ix_states.update([anc_ref, anc_alt])

    snp_df = df[list(snp_ix_states)]
    snp_df = snp_df[~snp_df.index.get_level_values("snp_id").duplicated()]
    # snp_df = df[list(snp_ix_states)].groupby(df.index.names).first()
    n_snps = len(snp_df.index.get_level_values("snp_id"))
    n_states = len(state_ids)

    if not doalphabeta:
        if prior is None and cont_id is not None:
            ca, cb = empirical_bayes_prior(snp_df[cont_ref], snp_df[cont_alt])

        P = Probs2(
            O=np.array(df.talt.values, np.uint8),
            N=np.array(df.tref.values + df.talt.values, np.uint8),
            P_cont=0.0
            if cont_id is None
            else np.array(
                (df[cont_ref] + ca) / (df[cont_ref] + df[cont_alt] + ca + cb)
            ),
            alpha=[],
            beta=[],
            alpha_hap=[],
            beta_hap=[],
            lib=np.array(df.lib),
        )
        return P

    if prior is None:  # empirical bayes, estimate from data
        alt_prior = np.empty((n_snps, n_states))
        ref_prior = np.empty((n_snps, n_states))
        if cont_id is not None:
            ca, cb = empirical_bayes_prior(snp_df[cont_ref], snp_df[cont_alt])

        if ancestral is None:
            for i, (a, b, s) in enumerate(zip(alt_ix, ref_ix, state_ids)):
                pa, pb = empirical_bayes_prior(snp_df[a], snp_df[b])
                logging.info("[%s]EB prior [a=%.4f, b=%.4f]: " % (s, pa, pb))
                alt_prior[:, i] = snp_df[a] + pa
                ref_prior[:, i] = snp_df[b] + pb
        else:
            anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"

            # set up vectors stating which allele is ancestral
            ref_is_anc = (snp_df[anc_ref] > 0) & (snp_df[anc_alt] == 0)
            alt_is_anc = (snp_df[anc_alt] > 0) & (snp_df[anc_ref] == 0)
            ref_is_der, alt_is_der = alt_is_anc, ref_is_anc
            anc_is_unknown = (1 - alt_is_anc) * (1 - ref_is_anc) == 1

            for i, (alt_col, ref_col, s) in enumerate(zip(alt_ix, ref_ix, state_ids)):

                # 1. set up base entries based on observed counts
                alt_prior[:, i] = snp_df[alt_col]
                ref_prior[:, i] = snp_df[ref_col]

                # 2. where anc is unknown, add symmetric prior estimated from data
                pa, pb = empirical_bayes_prior(snp_df[alt_col], snp_df[ref_col])
                logging.info("[%s]EB prior0 [anc=%.4f, der=%.4f]: " % (s, pa, pb))
                alt_prior[anc_is_unknown, i] += pa
                ref_prior[anc_is_unknown, i] += pb

                # 3. where anc is known, create indices
                m_anc = np.array(pd.concat((ref_is_anc, alt_is_anc), axis=1))
                m_der = np.array(pd.concat((ref_is_der, alt_is_der), axis=1))
                ANC = np.array(snp_df[[ref_col, alt_col]])[m_anc]
                DER = np.array(snp_df[[ref_col, alt_col]])[m_der]

                pder, panc = empirical_bayes_prior(DER, ANC, known_anc=True)
                panc += ancestral_prior
                logging.info("[%s]EB prior1 [anc=%.4f, der=%.4f]: " % (s, panc, pder))
                alt_prior[alt_is_anc, i] += panc
                alt_prior[alt_is_der, i] += pder
                ref_prior[ref_is_anc, i] += panc
                ref_prior[ref_is_der, i] += pder

        assert np.all(df.tref.values + df.talt.values < 256)

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

        alt_prior = snp_df[alt_ix].to_numpy() + prior_anc_alt[:, np.newaxis] + prior
        ref_prior = snp_df[ref_ix].to_numpy() + prior_anc_ref[:, np.newaxis] + prior
        assert np.all(df.tref.values + df.talt.values < 256)

    # create named tuple for return
    P = Probs2(
        O=np.array(df.talt.values, np.uint8),
        N=np.array(df.tref.values + df.talt.values, np.uint8),
        P_cont=0.0
        if cont_id is None
        else np.array((df[cont_ref] + ca) / (df[cont_ref] + df[cont_alt] + ca + cb)),
        alpha=alt_prior[IX.diploid_snps],
        beta=ref_prior[IX.diploid_snps],
        alpha_hap=alt_prior[IX.haploid_snps],
        beta_hap=ref_prior[IX.haploid_snps],
        lib=np.array(df.lib),
    )
    return P


def bins_from_bed(df, bin_size, sex=None, snp_mode=False):
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

        if snp_mode:
            bins = np.arange(bin0, bin0 + np.sum(snp.chrom == chrom))
            logging.debug("binning chrom %s: %d snp bins", chrom, len(bins))

            IX.bin_sizes.append(len(bins))
            bin_ids = bins
            _bin = np.empty_like(bins, dtype_bin)
            _bin["chrom"] = chrom
            _bin["pos"] = pos
            _bin["id"] = bin_ids
            _bin["map"] = bins
            _bin["haploid"] = chrom_is_hap
            bin_loc.append(_bin)

            snp_ids = bins
            IX.SNP2BIN[bins] = bins

        else:
            # create bins
            bins = np.arange(chrom_start, chrom_end, bin_size)
            logging.debug("binning chrom %s: %d bins", chrom, len(bins))

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

    logging.debug("done creating bins")
    return bins, IX


def get_haploid_stuff(snp, chroms, sex):
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

    haploid_snps = snp.snp_id[snp.chrom.isin(haplo_chroms)]
    if len(haploid_snps) > 0:
        haploid_snps = slice(min(haploid_snps), max(haploid_snps) + 1)
    else:
        haploid_snps = slice(0, 0)

    return haplo_chroms, haploid_snps


def make_obs2rg(v_lib):
    libs = np.unique(v_lib)
    libs = dict((l, i) for i, l in enumerate(libs))
    OBS2RG = np.array([libs[i] for i in v_lib], dtype=np.uint16)
    return libs, OBS2RG


def make_reads2rg(v_lib):
    libs = np.unique(v_lib)
    libs = dict((l, i) for i, l in enumerate(libs))
    READS2RG = np.array([libs[i] for i in v_lib], dtype=np.uint16)
    return libs, READS2RG


def make_obs2sfs(snp_states, max_states=None, states=None):
    # create dict [sfs] - > column
    if max_states is None:
        sfs = snp_states.reset_index(drop=True).drop_duplicates().reset_index(drop=True)
        data = snp_states
    else:
        data = pd.DataFrame()
        for s in states:
            freq = snp_states[f"{s}_alt"] / (
                snp_states[f"{s}_alt"] + snp_states[f"{s}_ref"]
            )
            freq = np.nan_to_num(freq, 0)
            m = np.max((snp_states[f"{s}_alt"] + snp_states[f"{s}_ref"]))
            m = min((m, max_states))
            freq = np.round(freq * m).astype(np.uint8)
            data[f"{s}_alt"], data[f"{s}_ref"] = freq, m - freq
            sfs = data.drop_duplicates().reset_index(drop=True)

    sfs_dict = dict((tuple(v.values()), k) for (k, v) in sfs.to_dict("index").items())
    data = np.array(data)
    SNP2SFS = np.array([sfs_dict[tuple(i)] for i in data], dtype=np.uint16)

    return sfs, SNP2SFS


def make_obs2sfs_folded(snp, ix_normal, anc_ref, anc_alt, max_states=None, states=None):
    """create sfs data structure taking ancestral allele into account

    basic strat
    1. create FLIPPED, which is true for SNP that need to be flipped, i.e.
        the ancestral allele is the alt-allele
    2. make dict[state] : index for all possible indices
    4. use dict to create SNP2SFS

    
    """

    """1. create FLIPPED, which is true for SNP that need to be flipped"""
    FLIPPED = (snp[anc_ref] == 0) & (snp[anc_alt] > 0)
    FLIPPED.reset_index(drop=True, inplace=True)
    snp.reset_index(drop=True, inplace=True)

    """2. make data1, data2 which are flipped/non-flipped data"""
    data1 = pd.DataFrame()
    data2 = pd.DataFrame()
    for s in states:
        data1[s] = snp.loc[~FLIPPED, f"{s}_alt"] / (
            snp.loc[~FLIPPED, f"{s}_alt"] + snp.loc[~FLIPPED, f"{s}_ref"]
        )
        data2[s] = snp.loc[FLIPPED, f"{s}_ref"] / (
            snp.loc[FLIPPED, f"{s}_alt"] + snp.loc[FLIPPED, f"{s}_ref"]
        )
        data1[s] = np.nan_to_num(data1[s], 0)
        data2[s] = np.nan_to_num(data2[s], 0)
        m = np.max((snp[f"{s}_alt"] + snp[f"{s}_ref"]))
        m = m if max_states is None else min((m, max_states))
        data1[s] = np.round(data1[s] * m).astype(np.uint8)
        data2[s] = np.round(data2[s] * m).astype(np.uint8)

        data1[f"{s}_alt"], data1[f"{s}_ref"] = data1[s], m - data1[s]
        data2[f"{s}_alt"], data2[f"{s}_ref"] = data2[s], m - data2[s]
        del data1[s], data2[s]

    data = pd.DataFrame(
        np.vstack((data1.to_numpy(), data2.to_numpy())), columns=data1.columns
    )
    data[~FLIPPED] = data1
    data[FLIPPED] = data2
    sfs = data.drop_duplicates().reset_index(drop=True)
    sfs_dict = dict((tuple(v.values()), k) for (k, v) in sfs.to_dict("index").items())
    data = np.array(data)
    """4. use dicts to create SNP2SFS"""
    SNP2SFS = np.array([sfs_dict[tuple(i)] for i in data], dtype=np.uint16)

    return sfs, SNP2SFS, FLIPPED


def make_slug_data(df, states, ancestral=None, cont_id=None, sex=None):
    """ create a SlugData object with the following attributes:
    N: np.ndarray  # number of reads [O x 1]
    O: np.ndarray  # number of obs [O x 1]
    psi: np.ndarray  # freq of contaminant [L x 1]

    OBS2RG: np.ndarray  # which obs is in which rg [O x 1]
    OBS2SNP: np.ndarray  # which obs belongs to which locus [O x 1]
    SNP2SFS: np.ndarray  # which snp belongs to which sfs entry [L x 1]

    haploid_snps : Any = None

    states : Any = None
    rgs : Any = None
    chroms : Any = None
    """
    ref_ix, alt_ix = [f"{s}_ref" for s in states], [f"{s}_alt" for s in states]
    sfs_state_ix = alt_ix + ref_ix  # states used in sfs
    all_state_ix = set(alt_ix + ref_ix)  # states in sfs + contamination + ancestral
    sfs_state_flipped = [
        v.replace("alt", "ref") if "alt" in v else v.replace("ref", "alt")
        for v in sfs_state_ix
    ]

    if cont_id is not None:
        cont_ref, cont_alt = f"{cont_id}_ref", f"{cont_id}_alt"
        all_state_ix.update([cont_ref, cont_alt])
    if ancestral is not None:
        anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"
        all_state_ix.update([anc_ref, anc_alt])

    rgs, OBS2RG = make_obs2rg(df.rg)

    snp = df[list(all_state_ix)].reset_index().drop_duplicates()

    chroms = pd.unique(snp.chrom)
    haplo_chroms, haplo_snps = get_haploid_stuff(snp, chroms, sex)

    if cont_id is None:
        psi = np.zeros_like(SNP2SFS)
    else:
        psi = snp[cont_alt] / (snp[cont_alt] + snp[cont_ref] + 1e-100)

    data = SlugData(
        ALT=df.talt.to_numpy(),
        REF=df.tref.to_numpy(),
        psi=psi,
        OBS2RG=OBS2RG,
        OBS2SNP=df.index.get_level_values("snp_id").to_numpy(),
        SNP2SFS=SNP2SFS,
        haploid_snps=haplo_snps,
        states=sfs_state_ix,
        rgs=rgs,
        sex=sex,
        chroms=chroms,
        haplo_chroms=haplo_chroms,
    )

    logging.debug("done creating data")
    return data, sfs


@njit
def make_full_df(df, n_reads):
    READS = np.empty(n_reads, np.uint8)
    READ2RG = np.empty(n_reads, np.uint32)
    READ2SNP = np.empty(n_reads, np.uint32)

    i = 0
    for snp_id, ref, alt, rg in df:
        for _ in range(ref):
            READS[i] = 0
            READ2RG[i] = rg
            READ2SNP[i] = snp_id
            i += 1
        for _ in range(alt):
            READS[i] = 1
            READ2RG[i] = rg
            READ2SNP[i] = snp_id
            i += 1

    return READS, READ2RG, READ2SNP


def make_slug_reads_data(
    df, states, max_states=8, ancestral=None, cont_id=None, sex=None, flip=True
):
    """ create a SlugReads object with the following attributes:
    """
    ref_ix, alt_ix = [f"{s}_ref" for s in states], [f"{s}_alt" for s in states]
    sfs_state_ix = alt_ix + ref_ix  # states used in sfs
    all_state_ix = set(alt_ix + ref_ix)  # states in sfs + contamination + ancestral

    if cont_id is not None:
        cont_ref, cont_alt = f"{cont_id}_ref", f"{cont_id}_alt"
        all_state_ix.update([cont_ref, cont_alt])
    if ancestral is not None:
        anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"
        all_state_ix.update([anc_ref, anc_alt])


    df2 = df.reset_index()[["snp_id", "tref", "talt", "rg"]]
    rgs = np.unique(df2.rg)
    rg_dict = dict((l, i) for i, l in enumerate(rgs))
    df2["rg"] = [rg_dict[rg] for rg in df2.rg]
    n_reads = np.sum(df2.tref + df2.talt)
    READS, READ2RG, READ2SNP = make_full_df(df2.to_numpy(), n_reads)

    assert np.sum(df.talt) == np.sum(READS == 1)
    assert np.sum(df.tref) == np.sum(READS == 0)

    snp = df[list(all_state_ix)].reset_index().drop_duplicates()

    if flip and ancestral is not None:
        sfs, SNP2SFS, FLIPPED = make_obs2sfs_folded(
            snp, sfs_state_ix, anc_ref, anc_alt, max_states, states
        )
    else:
        sfs, SNP2SFS = make_obs2sfs(snp[sfs_state_ix], max_states, states)
        FLIPPED = np.zeros_like(SNP2SFS, bool)

    chroms = pd.unique(snp.chrom)
    haplo_chroms, haplo_snps = get_haploid_stuff(snp, chroms, sex)

    if cont_id is None:
        psi = np.zeros_like(SNP2SFS)
    else:
        psi = snp[cont_alt] / (snp[cont_alt] + snp[cont_ref] + 1e-100)

    data = SlugReads(
        READS=READS,
        READ2RG=READ2RG,
        READ2SNP=READ2SNP,
        SNP2SFS=SNP2SFS,
        FLIPPED=FLIPPED,
        psi=psi,
        haploid_snps=haplo_snps,
        states=sfs_state_ix,
        rgs=rg_dict,
        sex=sex,
        chroms=chroms,
        haplo_chroms=haplo_chroms,
    )

    logging.debug("done creating data")
    return data, sfs


def init_ftau(n_states, F0=0.5, tau0=0):
    """initializes F and tau, which exist for each homozygous state
    """
    try:
        if len(F0) == n_states:
            F = F0
        elif len(F0) == 1:
            F = F0 * n_states
        else:
            F = [F0]
    except TypeError:
        F = [F0] * n_states
    try:
        if len(tau0) == n_states:
            tau = tau0
        elif len(F0) == 1:
            tau = tau0 * n_states
        else:
            tau = [tau0]
    except TypeError:
        tau = [tau0] * n_states

    return np.array(F), np.array(tau)


def init_ce(c0=0.01, e0=0.001):
    cont = defaultdict(lambda: c0)
    error = defaultdict(lambda: e0)
    return cont, error


def init_pars_sfs(n_sfs, n_rgs, F0, tau0, e0, c0, **kwargs):
    cont, error = init_ce(c0, e0)
    F, tau = init_ftau(n_sfs, F0, tau0)
    pars = SlugParsSquare(
        cont=np.zeros(n_rgs) + c0,
        tau=np.zeros(n_sfs) + tau0,
        F=np.zeros(n_sfs) + F0,
        e=e0,
        b=e0,
    )
    return pars


def trans_mat_hap_to_dip(tmat):
    """given a haploid transition rate matrix tmat, returns a diploid transition
    matrix

    assumptions:
     - only one transition at a time
     - "canonical" order of heterozygous states
     - independence between haplotypes
    """
    n = tmat.shape[0]
    n_homo = n
    n_het = int(n * (n - 1) / 2)
    tmat2 = np.zeros((n_homo + n_het, n_homo + n_het))

    # homo -> homo transition
    # for i in range(n_homo):
    #    for j in range(n_homo):
    #        tmat2[i, j] = tmat[i, j] ** 2 #prob both change at once

    # homo -> het transition
    for i in range(n_homo):
        c = n_homo  # state
        for h1 in range(n_homo):  # first het
            for h2 in range(h1 + 1, n_homo):  # second het
                if i == h1:  # only transition second haplotype
                    tmat2[i, c] = tmat[h1, h2]
                    tmat2[c, i] = tmat[h2, h1]
                elif i == h2:  # only transition first haplotype
                    tmat2[i, c] = tmat[h2, h1]
                    tmat2[c, i] = tmat[h1, h2]
                else:  # transition both
                    pass
                    # tmat2[i, c] = tmat[i, h1] * tmat[i, h2] + tmat[i, h2] * tmat[j, h1]
                    # tmat2[c, i] = tmat[h1,i] * tmat[h2, i] + tmat[h2, i] * tmat[h1, j]
                c += 1

    # het -> het transition
    c1 = n_homo
    for i in range(n_homo):  # first het from
        for j in range(i + 1, n_homo):  # second het from
            c2 = n_homo
            for h1 in range(n_homo):  # first het of target
                for h2 in range(h1 + 1, n_homo):  # second het of target
                    if i == h1 and j == h2:  # no transition
                        continue
                    elif i == h1:  # transition second haplotype from j to h2
                        tmat2[c1, c2] = tmat[j, h2]
                    elif j == h2:  # transition first haplotype from i to h1
                        tmat2[c1, c2] = tmat[i, h1]
                    else:  # transition both
                        pass
                    c2 += 1
            c1 += 1

    s = np.sum(tmat2, 1)
    for i in range(tmat2.shape[0]):
        tmat2[i, i] -= s[i]

    return tmat2


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
    transition_matrix=None,
    bin_size=1.0,
    **kwargs,
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

    if transition_matrix is None:
        trans_mat = np.zeros((n_states, n_states)) + 2e-2
        trans_mat_hap = np.zeros((n_hap, n_hap)) + 2e-2

        np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
        np.fill_diagonal(trans_mat_hap, 1 - (n_hap - 1) * 2e-2)

        if init_guess is not None:
            guess = [i for i, n in enumerate(gamma_names) if n in init_guess]
            logging.info("starting with guess %s " % guess)
            trans_mat[:, guess] = trans_mat[:, guess] + 1
            trans_mat /= np.sum(trans_mat, 1)[:, np.newaxis]
    else:
        trans_mat_hap = pd.read_csv(transition_matrix, header=None).to_numpy()
        trans_mat = trans_mat_hap_to_dip(trans_mat_hap)
        print(f"BS:{bin_size}")
        print(trans_mat)
        print(trans_mat_hap)
        trans_mat_hap = expm(trans_mat_hap * bin_size)
        trans_mat = expm(trans_mat * bin_size)
        print(trans_mat)
        print(trans_mat_hap)

    cont, error = init_ce(c0, e0)
    F, tau = init_ftau(n_homo, F0, tau0)

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
    try:
        random = np.random.binomial(1, np.clip(mu, 0, 1))
    except ValueError:
        breakpoint()
    PG = np.log10(PG + 1e-40)
    PG = np.minimum(0.0, PG)
    return pd.DataFrame(
        np.hstack((PG, mu, random)), columns=["G0", "G1", "G2", "p", "random_read"]
    )


def posterior_table_slug(pg, data, gtll=None):
    mu = np.sum(pg * np.arange(3) / 2.0, 1)
    random = np.random.binomial(1, np.clip(mu, 0, 1))
    log_g = np.log10(pg + 1e-40)
    log_g = np.minimum(0.0, log_g)
    df = np.hstack((log_g, mu[:, np.newaxis], random[:, np.newaxis]))
    df = pd.DataFrame(df, columns=["G0", "G1", "G2", "p", "random_read"])
    df.random_read = df.random_read.astype(np.uint8)
    if gtll is not None:
        log_ll = np.log10(gtll + 1e-40)
        df_ll = pd.DataFrame(log_ll, columns=["L0", "L1", "L2"])
        df = pd.concat((df, df_ll), axis=1)
    return df


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
        with an `X` or `Z` or `h`
    """
    ref["heterogametic"] = [
        v[0] in "XZxzh" for v in ref.index.get_level_values("chrom")
    ]
    data["heterogametic"] = [
        v[0] in "XZxzh" for v in data.index.get_level_values("chrom")
    ]

    n_sites = ref.groupby(ref.heterogametic).apply(lambda df: len(df))
    n_reads = data.groupby(data.heterogametic).apply(
        lambda df: np.sum(df.tref + df.talt)
    )
    cov = n_reads / n_sites
    del data["heterogametic"]
    del ref["heterogametic"]

    # no heteogametic data
    if True not in cov:
        return "f"

    if cov[True] / cov[False] < sex_ratio_threshold:
        sex = "m"
        logging.info("guessing sex is male, X/A = %.4f/%.4f" % (cov[True], cov[False]))
    else:
        sex = "f"
        logging.info(
            "guessing sex is female, X/A = %.4f/%.4f" % (cov[True], cov[False])
        )
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

    # zeros = np.isclose(scaling, 0)[:, 0, 0]
    M /= scaling
    # M[zeros] = 1
    # scaling[zeros] = 1

    assert np.allclose(np.max(M, (1, 2)), 1)
    log_scaling = np.sum(np.log(scaling))
    return log_scaling


def parse_state_string(states, ext=[""], operator="+", state_file=None):
    """parse shortcut parse strings

    the basic way states are defined is as 
        --states AFR NEA DEN

    this assmues the columns AFR, NEA and DEN are known and present in the reference
    We can rename and merge stuff here using the following syntax

    -- states AFR=YRI+SAN NEA=VIN DEN=Denisova2

    return rename dict for reference
    
    """
    ext2 = ["_ref", "_alt"]
    d1 = [s.split("=") for s in states if s is not None]
    d2 = [(s if len(s) > 1 else (s[0], s[0])) for s in d1]
    state_dict = dict(
        (f"{k}{ext_}", f"{i}{ext_}")
        for i, j in d2
        for k in j.split(operator)
        for ext_ in ext
    )
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
        # D = dict((v_, k) for (k,v) in Y.items() for v_ in v)
        s2 = dict()

        for state, label in state_dict.items():
            if state in Y:
                for ind in Y[state]:
                    s2[ind] = label
            else:
                s2[state] = label
        state_dict = s2
    return state_dict
