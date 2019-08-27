import yaml
from numba import njit
from collections import Counter
from scipy.stats import binom
import pandas as pd
import numpy as np
import itertools


try:
    from .utils import posterior_table, parse_state_string
    from .log import log_
except (ImportError, ModuleNotFoundError):
    from utils import posterior_table, parse_state_string
    from log import log_

"""reading and writing files"""


def load_ref(
    ref_files, state_dict, cont_id, ancestral=None, autosomes_only=False,
    map_col = 'map'
):
    """loads reference in custom (csv) format
    ref_files: paths to files
    state_dict: a dict D[label] : pop. All populatios in sources with label are added to pop. This is 
        expected to be generated with utils.parse_state_string
    """

    #1. get list of states we care about
    label_states = list(state_dict.keys())
    EXT = ['_ref', '_alt']
    D = dict(((k+e), (v+e)) for ((k, v), e) in itertools.product(state_dict.items(), EXT))


    if ancestral is not None:
        label_states = list(set(list(label_states) + [ancestral]))
    if cont_id is not None:
        label_states = list(set(list(label_states) + [cont_id]))

    #2. required in every ref
    basic_cols = ["chrom", "pos", "ref", "alt"]  # required in every ref

    #target states
    ref_cols = [f"{s}_ref" for s in label_states]
    alt_cols = [f"{s}_alt" for s in label_states]
    data_cols = ref_cols + alt_cols
    
    #which file a column is in
    file_ix = [None for i in data_cols]
    map_ix = None

    # read headers of each ref
    headers = list(list(pd.read_csv(r, nrows=0).columns) for r in ref_files)
    map_file = [i for i, h in enumerate(headers) if map_col in h]
    map_file = map_file[0] if len(map_file) > 0 else None


    for i, col in enumerate(data_cols):
        for j, h in enumerate(headers):
            if col in h and file_ix[i] is None:
                file_ix[i] = j
                log_.debug(f"found col {col} in header {j}")
                break

    if None in file_ix:
        s = [c for i, c in zip(file_ix, data_cols) if i is None]
        raise ValueError("columns not found in reference: " + ", ".join(s))

    #read correct cols from each file
    dtype_ = dict(chrom="category")
    for i, ref_file in enumerate(ref_files):

        cols0 = basic_cols + [col for ix, col in zip(file_ix, data_cols) if ix == i]
        if map_file == i:
            cols0 = cols0 + [map_col] 
            ix_cols = basic_cols  + [map_col]
        else:
            ix_cols = basic_cols

        ref0 = pd.read_csv(ref_file, dtype=dtype_, 
                           usecols=cols0,
                           index_col=ix_cols)
        ref0.index.rename('map', level=map_col, inplace=True)
        #ref0.chrom.cat.reorder_categories(pd.unique(ref0.chrom), inplace=True)
        ref0 = ref0.loc[~ref0.index.duplicated()]
        ref0 = ref0[~np.isnan(ref0.reset_index('map').map.values)]

        if i == 0:
            ref = ref0
        else:
            ref = ref.join(ref0)

    ref.fillna(value=0, inplace=True)


    #aggregate different labels
    ref = ref.rename(D, axis=1).groupby(level=0, axis=1).agg(sum)

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
        chrom = ref.index.get_level_values('chrom').factorize()[0]
        pos = ref.index.get_level_values('pos').values
        kp = nfp(chrom, pos, ref.shape[0], filter_pos)
        log_.info("filtering %s SNP due to pos filter", np.sum(1 - kp))
        ref = ref[kp]

    if filter_map is not None:
        chrom = ref.index.get_level_values('chrom').factorize()[0]
        pos = ref.index.get_level_values('map').values
        kp = nfp(chrom, pos, ref.shape[0], filter_map)
        log_.info("filtering %s SNP due to map filter", np.sum(1 - kp))
        ref = ref[kp]

    return ref


def load_read_data(infile, split_lib=True, downsample=1, high_cov_filter=0.001):
    dtype_ = dict(chrom="category")
    data = pd.read_csv(infile, dtype=dtype_, index_col=['chrom', 'pos']).dropna()
    #data.chrom.cat.reorder_categories(pd.unique(data.chrom), inplace=True)

    if "lib" not in data:
        data["lib"] = "lib0"
    elif not split_lib:
            data = data.groupby(data.index.names).agg(sum)

    if downsample < 1:
        data.tref = binom.rvs(data.tref, downsample, size=len(data.tref))
        data.talt = binom.rvs(data.talt, downsample, size=len(data.talt))

    # rm sites with extremely high coverage
    data = data[data.tref + data.talt > 0]
    q = np.quantile(data.tref + data.talt, 1 - high_cov_filter)
    data = data[data.tref + data.talt <= q]
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
