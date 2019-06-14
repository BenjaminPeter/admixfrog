import yaml
from numba import njit
from collections import Counter
from scipy.stats import binom
import pandas as pd
import numpy as np

try:
    from .utils import posterior_table
    from .log import log_
except (ImportError, ModuleNotFoundError):
    from utils import posterior_table
    from log import log_

"""reading and writing files"""


def load_ref(
    ref_files, state_ids, cont_id, prior=0, ancestral=None, autosomes_only=False
):
    states = list(set(list(state_ids)))
    if ancestral is not None:
        states = list(set(list(states) + [ancestral]))
    if cont_id is not None:
        states = list(set(list(states) + [cont_id]))

    basic_cols = ["chrom", "pos", "ref", "alt"]  # required in every ref
    map_col = ["map"]
    ref_cols = [f"{s}_ref" for s in states]
    alt_cols = [f"{s}_alt" for s in states]
    cols = map_col + ref_cols + alt_cols
    file_ix = [None for i in cols]

    # read headers of each ref
    headers = list(list(pd.read_csv(r, nrows=0).columns) for r in ref_files)
    for i, col in enumerate(cols):
        for j, h in enumerate(headers):
            if col in h and file_ix[i] is None:
                file_ix[i] = j
                log_.debug(f"found col {col} in header {j}")
                break

    if None in file_ix:
        s = [c for i, c in zip(file_ix, cols) if i is None]
        raise ValueError("columns not found in reference: " + ", ".join(s))

    dtype_ = dict(chrom="category")
    for i, ref_file in enumerate(ref_files):
        cols0 = basic_cols + [col for ix, col in zip(file_ix, cols) if ix == i]
        ref0 = pd.read_csv(ref_file, dtype=dtype_, usecols=cols0)
        ref0.chrom.cat.reorder_categories(pd.unique(ref0.chrom), inplace=True)
        ref0.drop_duplicates(inplace=True, subset=basic_cols)
        if i == 0:
            ref = ref0
        else:
            ref = ref.merge(ref0, on=basic_cols, how="outer")

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
        ref["HALF_ref"] = 0.5 - prior
        ref["HALF_alt"] = 0.5 - prior

    ref.dropna(inplace=True, subset=basic_cols + map_col)
    ref.iloc[:, 5:].fillna(value=0, inplace=True)

    if autosomes_only:
        ref = ref[ref.chrom != "X"]
        ref = ref[ref.chrom != "Y"]
        ref = ref[ref.chrom != "mt"]
    return ref


def filter_ref(ref, states, filter_delta=None, filter_pos=None, filter_map=None):
    n_states = len(states)

    if filter_delta is not None:
        kp = np.zeros(ref.shape[0], np.bool)
        for i, s1 in enumerate(states):
            for j in range(i + 1, n_states):
                s2 = states[j]
                f1 = np.nan_to_num(
                    ref[s1 + "_alt"] / (ref[s1 + "_alt"] + ref[s1 + "_ref"])
                )
                f2 = np.nan_to_num(
                    ref[s2 + "_alt"] / (ref[s2 + "_alt"] + ref[s2 + "_ref"])
                )
                delta = np.abs(f1 - f2)
                kp = np.logical_or(kp, delta >= filter_delta)

        log_.info("filtering %s SNP due to delta", np.sum(1 - kp))
        ref = ref[kp]

    if filter_pos is not None:
        chrom = pd.factorize(ref.chrom)[0]
        pos = np.array(ref.pos)
        kp = nfp(chrom, pos, ref.shape[0], filter_pos)
        log_.info("filtering %s SNP due to pos filter", np.sum(1 - kp))
        ref = ref[kp]

    if filter_map is not None:
        chrom = pd.factorize(ref.chrom)[0]
        pos = np.array(ref.map)
        kp = nfp(chrom, pos, ref.shape[0], filter_map)
        log_.info("filtering %s SNP due to map filter", np.sum(1 - kp))
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


class IndentDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(IndentDumper, self).increase_indent(flow, False)


def write_pars_table(pars, outname=None):
    P = dict(pars._asdict())
    for k, v in P.items():
        try:
            P[k] = v.tolist()
        except:
            pass

    s = yaml.dump(P, Dumper=IndentDumper, default_flow_style=False, indent=4)

    if outname is not None:
        with open(outname, "wt") as f:
            f.write(s)

    return s


def write_cont_table(data, cont, error, outname=None):

    df_libs = pd.DataFrame(cont.items(), columns=["lib", "cont"])
    df_error = pd.DataFrame(error.items(), columns=["lib", "error"])
    df_libs = df_libs.merge(df_error)

    rgs, deams, len_bins = [], [], []
    for l in df_libs.lib:
        try:
            rg, len_bin, deam = l.split("_")
        except ValueError:
            rg, len_bin, deam = l, 0, "NA"
        rgs.append(rg)
        deams.append(deam)
        len_bins.append(len_bin)
    df_libs["rg"] = rgs
    df_libs["len_bin"] = len_bins
    df_libs["deam"] = deams

    CC = data.groupby(["lib"]).agg(({"tref": sum, "talt": sum})).reset_index()
    CC["n_snps"] = CC.tref + CC.talt
    del CC["tref"]
    del CC["talt"]

    df_libs = df_libs.merge(CC)
    df_libs.sort_values("n_snps", ascending=False)

    if outname is not None:
        df_libs.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return df_libs


def write_bin_table(Z, bins, viterbi_df, gamma_names, IX, outname=None):
    df_bin = pd.DataFrame(Z, columns=gamma_names)
    CC = Counter(IX.SNP2BIN)
    snp = pd.DataFrame([CC[i] for i in range(len(df_bin))], columns=["n_snps"])
    df_bin = pd.concat((pd.DataFrame(bins), viterbi_df, snp, df_bin), axis=1)

    if outname is not None:
        df_bin.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return df_bin


def write_snp_table(data, G, Z, IX, gt_mode=False, outname=None):
    D = (
        data.groupby(["chrom", "pos", "map"])
        .agg({"tref": sum, "talt": sum})
        .reset_index()
    )
    if gt_mode:
        snp_df = pd.concat((D, pd.DataFrame(IX.SNP2BIN, columns=["bin"])), axis=1)
    else:
        T = posterior_table(G, Z, IX)
        snp_df = pd.concat((D, T, pd.DataFrame(IX.SNP2BIN, columns=["bin"])), axis=1)
    if outname is not None:
        snp_df.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return snp_df


def write_out_ref(data, G, Z, IX, gt_mode=False, outname=None):
    D = (
        data.groupby(["chrom", "pos", "map"])
        .agg({"tref": sum, "talt": sum})
        .reset_index()
    )
    if gt_mode:
        snp_df = pd.concat((D, pd.DataFrame(IX.SNP2BIN, columns=["bin"])), axis=1)
    else:
        T = posterior_table(G, Z, IX)
        snp_df = pd.concat((D, T, pd.DataFrame(IX.SNP2BIN, columns=["bin"])), axis=1)
    if outname is not None:
        snp_df.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return snp_df


def write_est_runs(df, outname=None):
    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False, compression="xz")


def write_sim_runs(df, outname=None):
    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False, compression="xz")
