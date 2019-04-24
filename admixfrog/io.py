from numba import njit
from scipy.stats import binom
import pandas as pd
import numpy as np
try:
    from .log import log_
except (ImportError, ModuleNotFoundError):
    from log import log_

"""reading and writing files"""

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

def filter_ref(ref, states, 
               filter_delta = None, 
               filter_pos = None,
               filter_map=None):
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

    if filter_map is not None:
        chrom = pd.factorize(ref.chrom)[0]
        pos = np.array(ref.map)
        kp = nfp(chrom ,pos, ref.shape[0], filter_map)
        log_.info("filtering %s SNP due to map filter", np.sum(1-kp))
        ref = ref[kp]

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
    """
    load genotype data. Currently program just uses load_read_data
    """
    dtype_ = dict(chrom="category")
    data = pd.read_csv(infile, dtype=dtype_).dropna()
    data.chrom.cat.reorder_categories(pd.unique(data.chrom), inplace=True)

    # rm sites with no coverage
    data = data[data.tref + data.talt > 0]
    return data


@njit
def nfp(chrom, pos, n_snps, filter_pos):
    kp = np.ones(n_snps, np.bool_)
    prev_chrom, prev_pos = -1, -10000000
    for i in range(n_snps):
        if prev_chrom != chrom[i]:
            prev_chrom, prev_pos = chrom[i], pos[i]
            continue
        if pos[i] - prev_pos <= filter_pos:
            kp[i] = False
        else:
            prev_pos = pos[i]

    return kp
