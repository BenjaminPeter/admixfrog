from collections import namedtuple, defaultdict
import numpy as np
import pandas as pd
from scipy.stats import binom

Probs = namedtuple("Probs", ("O", "N", "P_cont", "alpha", "beta", "lib"))
Pars = namedtuple(
    "Pars", ("alpha0", "trans_mat", "cont", "e0", "F", "gamma_names", "sex")
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
    state_priors=(1e-5, 1e-5),
    cont_prior=(1e-8, 1e-8),
    ancestral=None,
):
    alpha_ix = ["%s_alt" % s for s in state_ids]
    beta_ix = ["%s_ref" % s for s in state_ids]
    if ancestral is None:
        pa, pb = state_priors
    else:
        # anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"
        anc_ref, anc_alt = ancestral + "_ref", "ancestral" + "_alt"
        pa = data[anc_alt] + state_priors[0] * (1 - 2 * np.sign(data[anc_alt]))
        pb = data[anc_ref] + state_priors[1] * (1 - 2 * np.sign(data[anc_ref]))
    cont = "%s_alt" % cont_id, "%s_ref" % cont_id
    ca, cb = cont_prior

    print(alpha_ix, beta_ix)

    P = Probs(
        O=np.array(data.talt),
        N=np.array(data.tref + data.talt),
        P_cont=np.array(
            (data[cont[0]] + ca) / (data[cont[0]] + data[cont[1]] + ca + cb)
        ),
        alpha=np.array(ref[alpha_ix]) + pa,
        beta=np.array(ref[beta_ix]) + pb,
        lib=np.array(data.lib),
    )
    return P


def bins_from_bed(bed, data, bin_size, sex=None, pos_mode=False):
    """create a bunch of auxillary data frames for binning

    - bins: columns are chrom_id, chrom, bin_pos, bin_id, map
    - IX: container storing all indices, will need to be cleaned later on
    """
    IX = _IX()
    libs = np.unique(data.lib)
    if pos_mode:
        bed.map = bed.pos
    chroms = pd.unique(bed.chrom)
    snp = data[["chrom", "pos", "map"]].drop_duplicates()
    n_snps = snp.shape[0]
    snp["snp_id"] = range(n_snps)
    snp["hap"] = False

    if sex == "m":
        snp.loc[
            (data.chrom == "X") & (HAPX[0] < data.pos) & (data.pos < HAPX[1]), "hap"
        ] = True

    IX.SNP2CHROMBIN = np.empty((n_snps, 2), int)
    IX.SNP2BIN = np.empty((n_snps), int)
    IX.hapsnp = np.zeros(n_snps, bool)

    data = data.merge(snp)
    IX.OBS2SNP = np.array(data["snp_id"])

    bin_loc = []
    bin0 = 0

    dtype_bin = np.dtype(
        [
            ("chrom", "U2"),
            ("map", float),
            ("pos", int),
            ("id", int),
            ("chrom_id", int),
            ("hap", bool),
        ]
    )

    IX.bin_sizes = []

    for i, chrom in enumerate(chroms):
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
        _bin["chrom_id"] = i

        if chrom == "X" and sex == "m":
            _bin["hap"] = (HAPX[0] < _bin["pos"]) & (HAPX[1] > _bin["pos"])
        else:
            _bin["hap"] = False
        bin_loc.append(_bin)

        snp_ids = snp.snp_id[snp.chrom == chrom]
        dig_snp = np.digitize(snp[snp.chrom == chrom].map, bins, right=False) - 1
        IX.SNP2CHROMBIN[snp_ids, 0] = i
        IX.SNP2CHROMBIN[snp_ids, 1] = dig_snp
        IX.SNP2BIN[snp_ids] = dig_snp + bin0

        bin0 += len(bins)

    bins = np.hstack(bin_loc)

    IX.RG2OBS = dict((l, np.where(data.lib == l)[0]) for l in libs)
    IX.RG2SNP = dict((k, IX.OBS2SNP[v]) for k, v in IX.RG2OBS.items())
    IX.RG2BIN = dict((k, IX.SNP2BIN[v]) for k, v in IX.RG2SNP.items())
    IX.OBS2RG = np.array(data.lib)
    IX.OBS2BIN = IX.SNP2BIN[IX.OBS2SNP]
    IX.OBS2CHROMBIN = IX.SNP2CHROMBIN[IX.OBS2SNP]

    IX.HAPOBS = np.where(data.hap)[0]
    IX.HAPSNP = np.unique(IX.OBS2SNP[IX.HAPOBS])
    IX.DIPOBS = np.where(np.logical_not(data.hap))[0]
    IX.DIPSNP = np.unique(IX.OBS2SNP[IX.DIPOBS])
    IX.HAPBIN = bins["id"][bins["hap"]]
    assert all(x in IX.HAPBIN for x in IX.SNP2BIN[IX.HAPSNP])
    assert all(x in IX.HAPBIN for x in IX.OBS2BIN[IX.HAPOBS])

    IX.n_chroms = len(chroms)
    IX.n_bins = len(bins)
    IX.n_snps = len(IX.SNP2BIN)
    IX.n_obs = len(IX.OBS2SNP)
    IX.n_reads = np.sum(data.tref + data.talt)

    return bins, IX  # , data_bin


def init_pars(state_ids, sex=None, F0=0.001, e0=1e-2, c0=1e-2):
    homo = [s for s in state_ids]
    het = []
    for i, s in enumerate(state_ids):
        for s2 in state_ids[i + 1 :]:
            het.append(s + s2)
    gamma_names = homo + het
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
    return Pars(alpha0, trans_mat, cont, e0, F, gamma_names, sex=sex)


def load_ref(
    ref_file, state_ids, cont_id, prior=0, ancestral=None, autosomes_only=False
):
    if ancestral is None:
        states = list(set(list(state_ids) + [cont_id]))
    else:
        states = list(set(list(state_ids) + [cont_id, ancestral]))
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


def load_data(infile, split_lib=True, downsample=1):
    dtype_ = dict(chrom="category")
    data = pd.read_csv(infile, dtype=dtype_).dropna()
    data.chrom.cat.reorder_categories(pd.unique(data.chrom), inplace=True)
    if "lib" not in data or (not split_lib):
        data = data.groupby(["chrom", "pos"], as_index=False).agg(
            {"tref": sum, "talt": sum}
        )
        data["lib"] = "lib0"

    # rm sites with extremely high coverage
    if downsample < 1:
        data.tref = binom.rvs(data.tref, downsample, size=len(data.tref))
        data.talt = binom.rvs(data.talt, downsample, size=len(data.talt))

    data = data[data.tref + data.talt > 0]
    q = np.quantile(data.tref + data.talt, 0.999)
    data = data[data.tref + data.talt <= q]
    return data


def posterior_table(pg, Z, IX):
    PG = np.sum(Z[IX.SNP2BIN][:, :, np.newaxis] * pg, 1)  # genotype probs
    mu = np.sum(PG * np.arange(3), 1)[:, np.newaxis] / 2
    # musq = np.sum(PG * (np.arange(3))**2,1)[:, np.newaxis]
    # ahat = (2 *mu - musq) / (2 * (mu / musq - mu - 1) + mu)
    # bhat = (2 - mu)*(2- mu/musq) / (2 * (mu / musq - mu - 1) + mu)
    PG = np.log10(PG + 1e-40)
    PG = np.minimum(0.0, PG)
    return pd.DataFrame(np.hstack((PG, mu)), columns=["G0", "G1", "G2", "p"])
