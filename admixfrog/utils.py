from collections import namedtuple, defaultdict, Counter
import pdb
import numpy as np
from numba import njit
import pandas as pd
try:
    from .log import log_
except (ImportError, ModuleNotFoundError):
    from log import log_


from scipy.stats import binom

Probs = namedtuple("Probs", ("O", "N", "P_cont", "alpha", "beta", "lib"))
Pars = namedtuple(
    "Pars", ("alpha0", "trans_mat", "cont", "e0", "F", "tau", "gamma_names", "sex")
)  # for returning from varios functions
HAPX = (2699520, 155260560)  # start, end of haploid region on X


class _IX:
    def __init__(self):
        pass


def data2probs(
    data,
    ref,
    state_ids,
    cont_id,
    prior=None,
    cont_prior=(1e-8, 1e-8),
    ancestral=None,
):
    n_states = len(state_ids)
    n_snps = ref.shape[0]
    alpha_ix = ["%s_alt" % s for s in state_ids]
    beta_ix = ["%s_ref" % s for s in state_ids]

    cont = "%s_alt" % cont_id, "%s_ref" % cont_id
    if ancestral is None:
        data = data.merge(ref[list(cont) + ['chrom', 'pos', 'map']] )
    else:
        anc = "%s_alt" % ancestral, "%s_ref" % ancestral
        data = data.merge(ref[list(cont) + list(anc) + ['chrom', 'pos', 'map']] )

    if prior is None:    # empirical bayes
        alpha = np.empty((n_snps, n_states))
        beta = np.empty((n_snps, n_states))
        ca, cb = empirical_bayes_prior(ref[cont[0]], ref[cont[1]])

        if ancestral is None:
            for i, (a, b, s) in enumerate(zip(alpha_ix, beta_ix, state_ids)):
                pa, pb = empirical_bayes_prior(ref[a], ref[b])
                log_.info("[%s]EB prior [a=%.4f, b=%.4f]: " % (s, pa, pb))
                alpha[:, i] = ref[a] + pa
                beta[:, i] = ref[b] + pb
        else:
            anc_ref, anc_alt = ancestral + "_ref", ancestral + "_alt"
            ref_is_anc = (ref[anc_ref] == 1) & (ref[anc_alt] == 0)
            alt_is_anc = (ref[anc_alt] == 1) & (ref[anc_ref] == 0)
            ref_is_der, alt_is_der = alt_is_anc, ref_is_anc
            anc_is_unknown = (1 - alt_is_anc) * (1 - ref_is_anc) == 1
            for i, (a, b, s) in enumerate(zip(alpha_ix, beta_ix, state_ids)):
                pa, pb = empirical_bayes_prior(ref[a], ref[b])
                log_.info("[%s]EB prior0 [anc=%.4f, der=%.4f]: " % (s, pa, pb))
                alpha[:, i], beta[:, i] = ref[a], ref[b]
                alpha[anc_is_unknown, i] += pa
                beta[anc_is_unknown, i] += pb

                m_anc = pd.concat((ref_is_anc, alt_is_anc), 1)
                m_der = pd.concat((ref_is_der, alt_is_der), 1)
                ANC = np.array(ref[[b, a]])[m_anc]
                DER = np.array(ref[[b, a]])[m_der]

                pder, panc = empirical_bayes_prior(DER, ANC, True)
                log_.info("[%s]EB prior1 [anc=%.4f, der=%.4f]: " % (s, panc, pder))
                alpha[alt_is_anc, i] += panc
                alpha[alt_is_der, i] += pder
                beta[ref_is_anc, i] += panc
                beta[ref_is_der, i] += pder

        P = Probs(
            O=np.array(data.talt, np.int8),
            N=np.array(data.tref + data.talt, np.int8),
            P_cont=np.array(
                (data[cont[0]] + ca) / (data[cont[0]] + data[cont[1]] + ca + cb)
            ),
            alpha=alpha,
            beta=beta,
            lib=np.array(data.lib),
        )
        return P

    else:
        if ancestral is None:
            pass
        else:
            # anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"
            anc_ref, anc_alt = ancestral + "_ref", ancestral + "_alt"
            pa = data[anc_alt] + prior * (1 - 2 * np.sign(data[anc_alt]))
            pb = data[anc_ref] + prior * (1 - 2 * np.sign(data[anc_ref]))
        cont = "%s_alt" % cont_id, "%s_ref" % cont_id
        ca, cb = cont_prior


        print(alpha_ix)
        P = Probs(
            O=np.array(data.talt, np.int8),
            N=np.array(data.tref + data.talt, np.int8),
            P_cont=np.array(
                (data[cont[0]] + ca) / (data[cont[0]] + data[cont[1]] + ca + cb)
            ),
            alpha=np.array(ref[alpha_ix]) + prior,
            beta=np.array(ref[beta_ix]) + prior,
            lib=np.array(data.lib)
        )
        return P


def bins_from_bed(bed, snp, data, bin_size, sex=None, pos_mode=False,
                  ld_weighting=False):
    """create a bunch of auxillary data frames for binning

    - bins: columns are chrom_id, chrom, bin_pos, bin_id, map
    - IX: container storing all indices, will need to be cleaned later on
    """
    IX = _IX()
    libs = np.unique(data.lib)
    if pos_mode:
        bed.map = bed.pos
    chroms = pd.unique(bed.chrom)

    if sex == "m":
        snp.loc[
            (data.chrom == "X") & (HAPX[0] < data.pos) & (data.pos < HAPX[1]), "hap"
        ] = True
    n_snps = snp.shape[0]


    IX.SNP2BIN = np.empty((n_snps), int)
    IX.OBS2SNP = np.array(data["snp_id"])

    bin_loc = []
    bin0 = 0

    dtype_bin = np.dtype(
        [
            ("chrom", "U2"),
            ("map", float),
            ("pos", int),
            ("id", int),
        ]
    )

    IX.bin_sizes = []

    for i, chrom in enumerate(chroms):
        log_.debug("binning chrom %s", chrom)
        map_ = bed.map[bed.chrom == chrom]
        pos = bed.pos[bed.chrom == chrom]
        map_data = data.map[data.chrom == chrom]

        chrom_start = float(np.floor(map_.head(1) / bin_size) * bin_size)
        chrom_end = float(np.ceil(map_.tail(1) / bin_size) * bin_size)

        bins = np.arange(chrom_start, chrom_end, bin_size)
        IX.bin_sizes.append(len(bins))
        bin_ids = range(bin0, bin0 + len(bins))
        _bin = np.empty_like(bins, dtype_bin)
        _bin["chrom"] = chrom
        _bin["pos"] = np.interp(bins, map_, pos)
        _bin["id"] = bin_ids
        _bin["map"] = bins

        #if chrom == "X" and sex == "m":
        #    _bin["hap"] = True
        #else:
        #    _bin["hap"] = False
        bin_loc.append(_bin)

        snp_ids = snp.snp_id[snp.chrom == chrom]
        dig_snp = np.digitize(snp[snp.chrom == chrom].map, bins, right=False) - 1
        IX.SNP2BIN[snp_ids] = dig_snp + bin0

        bin0 += len(bins)

    bins = np.hstack(bin_loc)

    IX.RG2OBS = dict((l, np.where(data.lib == l)[0]) for l in libs)
    IX.OBS2BIN = IX.SNP2BIN[IX.OBS2SNP]
    IX.HAPSNP = []
    IX.HAPBIN = []

    #IX.HAPOBS = np.where(data.hap)[0]
    #IX.HAPSNP = np.unique(IX.OBS2SNP[IX.HAPOBS])
    #IX.DIPOBS = np.where(np.logical_not(data.hap))[0]
    #IX.DIPSNP = np.unique(IX.OBS2SNP[IX.DIPOBS])
    #IX.HAPBIN = bins["id"][bins["hap"]]
    #assert all(x in IX.HAPBIN for x in IX.SNP2BIN[IX.HAPSNP])
    #assert all(x in IX.HAPBIN for x in IX.OBS2BIN[IX.HAPOBS])

    IX.n_chroms = len(chroms)
    IX.n_bins = len(bins)
    IX.n_snps = len(IX.SNP2BIN)
    IX.n_obs = len(IX.OBS2SNP)
    IX.n_reads = np.sum(data.tref + data.talt)

    # ld weighting:
    # IX.snp_weight give for each SNP how it is supposed to be downweighted
    IX.snp_weight = np.ones(IX.n_snps)
    if ld_weighting:
        ctr = Counter(IX.SNP2BIN)  # counter[bin] : n_snp
        c = 0
        for i in range(IX.n_bins):
            for j in range(ctr[i]):
                IX.snp_weight[c] = 1. / ctr[i]
                c += 1
        log_.debug("mean ld-weight %s" % np.mean(IX.snp_weight))

    log_.debug("done creating bins")
    return bins, IX  # , data_bin


def init_pars(
    state_ids, sex=None, F0=0.001, tau0=1, e0=1e-2, c0=1e-2, est_inbreeding=False
):
    homo = [s for s in state_ids]
    het = []
    for i, s in enumerate(state_ids):
        for s2 in state_ids[i + 1 :]:
            het.append(s + s2)
    gamma_names = homo + het
    if est_inbreeding:
        gamma_names.extend(["h%s" % s for s in homo])

    n_states = len(gamma_names)
    n_homo = len(state_ids)

    alpha0 = np.array([1 / n_states] * n_states)
    trans_mat = np.zeros((n_states, n_states)) + 2e-2
    np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
    cont = defaultdict(lambda: c0)
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
    return Pars(alpha0, trans_mat, cont, e0, F, tau, gamma_names, sex=sex)


def filter_ref(ref, states, filter_delta = None, filter_pos = None):
    n_states = len(states)


    if filter_delta is not None:
        kp = np.zeros(ref.shape[0], np.bool)
        for i, s1 in enumerate(states):
            for j in range(i+1, n_states):
                s2 = states[j]
                f1 = np.nan_to_num(ref[s1 + "_alt"] / (ref[s1 + "_alt"] + ref[s1 + "_ref"]))
                f2 = np.nan_to_num(ref[s2 + "_alt"] / (ref[s2 + "_alt"] + ref[s2 + "_ref"]))
                delta = np.abs(f1 -f2)
                kp = np.logical_or(kp, delta >= filter_delta)

        log_.info("filtering %s SNP due to delta", np.sum(1-kp))
        ref = ref[kp]

    if filter_pos is not None:
        chrom = pd.factorize(ref.chrom)[0]
        pos = np.array(ref.pos)
        kp = nfp(chrom ,pos, ref.shape[0], filter_pos)
        log_.info("filtering %s SNP due to pos filter", np.sum(1-kp))
        ref = ref[kp]

    return ref

@njit
def nfp(chrom, pos, n_snps, filter_pos):
    kp = np.ones(n_snps, np.bool_)
    prev_chrom, prev_pos = -1, -10000000
    for i in range(n_snps):
        if prev_chrom != chrom[i]:
            prev_chrom, prev_pos = chrom[i], pos[i]
            continue
        if pos[i] - prev_pos < filter_pos:
            kp[i] = False
        else:
            prev_pos = pos[i]

    return kp


def load_ref(
    ref_file, state_ids, cont_id, prior=0, ancestral=None, autosomes_only=False
):
    states = list(set(list(state_ids)))
    if ancestral is not None:
        states = list(set(list(states) + [ancestral]))
    if cont_id is not None:
        states = list(set(list(states) + [cont_id]))

    dtype_ = dict(chrom="category")
    ref = pd.read_csv(ref_file, dtype=dtype_)
    ref.chrom.cat.reorder_categories(pd.unique(ref.chrom), inplace=True)

    if "UNIF" in states:
        ref["UNIF_ref"] = 1 - prior
        ref["UNIF_alt"] = 1 - prior
    if "REF" in states:
        ref["REF_ref"] = 1
        ref["REF_alt"] = 0
    if "NRE" in states:
        ref["NRE_ref"] = 0
        ref["NRE_alt"] = 1
    if "ZERO" in states:
        ref["ZERO_ref"] = 1e-7 - prior
        ref["ZERO_alt"] = 1e-7 - prior
    if "SFS" in states:
        ref["SFS_ref"] = prior
        ref["SFS_alt"] = prior
    if "HALF" in states:
        ref["SFS_ref"] = 0.5 - prior
        ref["SFS_alt"] = 0.5 - prior
    if "PAN" in states:
        ref["PAN_ref"] /= 2
        ref["PAN_alt"] /= 2

    ix = list(ref.columns[:5])
    suffixes = ["_alt", "_ref"]
    cols = ix + [s + x for s in states for x in suffixes]
    ref = ref[cols].dropna()
    if autosomes_only:
        ref = ref[ref.chrom != "X"]
        ref = ref[ref.chrom != "Y"]
    return ref


def load_read_data(infile, split_lib=True, downsample=1):
    dtype_ = dict(chrom="category")
    data = pd.read_csv(infile, dtype=dtype_).dropna()
    data.chrom.cat.reorder_categories(pd.unique(data.chrom), inplace=True)

    if "lib" not in data or (not split_lib):
        data = data.groupby(["chrom", "pos"], as_index=False).agg(
            {"tref": sum, "talt": sum}
        )
        data["lib"] = "lib0"

    if downsample < 1:
        data.tref = binom.rvs(data.tref, downsample, size=len(data.tref))
        data.talt = binom.rvs(data.talt, downsample, size=len(data.talt))

    # rm sites with extremely high coverage
    data = data[data.tref + data.talt > 0]
    q = np.quantile(data.tref + data.talt, 0.999)
    data = data[data.tref + data.talt <= q]
    return data


def load_gt_data(infile):
    dtype_ = dict(chrom="category")
    data = pd.read_csv(infile, dtype=dtype_).dropna()
    data.chrom.cat.reorder_categories(pd.unique(data.chrom), inplace=True)

    # rm sites with extremely high coverage
    data = data[data.tref + data.talt > 0]
    return data


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
    # H = ref * alt / n / n
    H = f * (1.0 - f)  # alt formulation
    if H == 0.0:
        return 1e-6, 1e-6

    V = np.nanvar(der / n) if known_anc else np.nanvar(np.hstack((der / n, anc / n)))
    ab = (H - V) / (V - H / np.nanmean(n))
    pa = max((f * ab, 1e-5))
    pb = max(((1-f) * ab, 1e-5))
    return pa, pb


def guess_sex(data):
    cov = data.groupby(data.chrom == "X").apply(
        lambda df: np.sum(df.tref + df.talt)
    )
    cov = cov.astype(float)
    cov[True] /= np.sum(data.chrom == "X")
    cov[False] /= np.sum(data.chrom != "X")

    if cov[True] / cov[False] < 0.8:
        sex = "m"
        log_.info("guessing sex is male, %.4f/%.4f" % (cov[True], cov[False]))
    else:
        sex = "f"
        log_.info("guessing sex is female, %.4f/%.4f" % (cov[True], cov[False]))
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

