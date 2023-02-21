import logging
from collections import namedtuple, defaultdict, Counter
import numpy as np
import pandas as pd
import yaml
from .classes import SlugData, FrogData
from .pars import FrogPars, SlugPars
from numba import njit
from scipy.linalg import expm




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
            freq = np.nan_to_num(freq)
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
        data1[s] = np.nan_to_num(data1[s])
        data2[s] = np.nan_to_num(data2[s])
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


def make_slug_data(
    df, states, max_states=8, ancestral=None, cont_id=None, sex=None, flip=True
):
    """create a SlugReads object with the following attributes:"""
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

    snp = (
        df[list(all_state_ix)]
        .reset_index("rg", drop=True)
        .reset_index()
        .drop_duplicates()
    )

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

    data = SlugData(
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


def init_pars_sfs(n_sfs, n_rgs, F0, tau0, e0, c0, **kwargs):
    pars = SlugPars(
        cont=np.zeros(n_rgs) + c0,
        tau=np.zeros(n_sfs) + tau0,
        F=np.zeros(n_sfs) + F0,
        e=np.zeros(1) + e0,
        b=np.zeros(1) + e0,
    )
    return pars



def posterior_table(pg, Z, P, est_inbreeding=False):
    freq = np.array([0, 1, 2, 0, 1]) if est_inbreeding else np.arange(3)
    PG = np.sum(Z[P.SNP2BIN][:, :, np.newaxis] * pg, 1)  # genotype probs
    mu = np.sum(PG * freq, 1)[:, np.newaxis] / 2
    random = np.random.binomial(1, np.clip(mu, 0, 1))
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
    M /= scaling

    assert np.allclose(np.max(M, (1, 2)), 1)
    log_scaling = np.sum(np.log(scaling))
    return log_scaling


def filter_ref(
    ref,
    states,
    filter_delta=None,
    filter_pos=None,
    filter_map=None,
    filter_ancestral=False,
    filter_cont=True,
    **kwargs,
):

    if filter_ancestral and states.ancestral is not None:
        no_ancestral_call = ref[f"{states.ancestral}_ref"] + ref[f"{states.ancestral}_alt"] == 0
        ref = ref.loc[~no_ancestral_call]
        logging.info(
            "filtering %s SNP due to missing ancestral call", np.sum(no_ancestral_call)
        )

    if filter_cont and states.contamination is not None:
        no_cont_data = ref[f"{states.contamination}_ref"] + ref[f"{states.contamination}_alt"] == 0
        ref = ref.loc[~no_cont_data]
        logging.info(
            "filtering %s SNP due to missing contaminant call", np.sum(no_cont_data)
        )

    if filter_delta is not None:
        kp = np.zeros(ref.shape[0], bool)
        for i, s1 in enumerate(states):
            for j in range(i + 1, states.n_raw_states):
                s2 = states[j]
                f1 = np.nan_to_num(
                    ref[s1 + "_alt"] / (ref[s1 + "_alt"] + ref[s1 + "_ref"])
                )
                f2 = np.nan_to_num(
                    ref[s2 + "_alt"] / (ref[s2 + "_alt"] + ref[s2 + "_ref"])
                )
                delta = np.abs(f1 - f2)
                kp = np.logical_or(kp, delta >= filter_delta)

        logging.info("filtering %s SNP due to delta", np.sum(1 - kp))
        ref = ref[kp]

    if filter_pos is not None and filter_pos >= 0:
        chrom = ref.index.get_level_values("chrom").factorize()[0]
        pos = ref.index.get_level_values("pos").values
        kp = nfp(chrom, pos, ref.shape[0], filter_pos)
        logging.info("filtering %s SNP due to pos filter", np.sum(1 - kp))
        ref = ref[kp]

    if filter_map is not None and filter_map >= 0:
        chrom = ref.index.get_level_values("chrom").factorize()[0]
        pos = ref.index.get_level_values("map").values
        kp = nfp(chrom, pos, ref.shape[0], filter_map)
        logging.info("filtering %s SNP due to map filter", np.sum(1 - kp))
        ref = ref[kp]

    return ref


@njit
def nfp(chrom, pos, n_snps, filter_pos):
    kp = np.ones(n_snps, np.bool_)
    prev_chrom, prev_pos = -1, -10_000_000
    for i in range(n_snps):
        if prev_chrom != chrom[i]:
            prev_chrom, prev_pos = chrom[i], pos[i]
            continue
        if pos[i] - prev_pos <= filter_pos:
            kp[i] = False
        else:
            prev_pos = pos[i]

    return kp
