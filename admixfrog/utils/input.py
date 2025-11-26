import logging
import yaml
from numba import njit
from collections import Counter
from scipy.stats import binom
import pandas as pd
import numpy as np
import itertools

from .utils import posterior_table, posterior_table_slug, parse_chroms


def load_ref(
    ref_files,
    states,
    cont_id,
    ancestral=None,
    autosomes_only=False,
    map_col="map",
    large_ref=True,
    sex_chroms=["X", "Y", "Z", "W", "mt"],
):
    """loads reference in custom (csv) format
    ref_files: paths to files
    """

    # 1. get list of states we care about
    label_states = list(states.state_dict.keys())
    EXT = ["_ref", "_alt"]
    D = dict(
        ((k + e), (v + e))
        for ((k, v), e) in itertools.product(states.state_dict.items(), EXT)
    )

    if ancestral is not None:
        label_states = list(set(list(label_states) + [ancestral]))
    if cont_id is not None:
        label_states = list(set(list(label_states) + [cont_id]))

    # 2. required in every ref
    basic_cols = ["chrom", "pos", "ref", "alt"]  # required in every ref

    # target states
    ref_cols = [f"{s}_ref" for s in label_states]
    alt_cols = [f"{s}_alt" for s in label_states]
    data_cols = ref_cols + alt_cols

    dtype_ = dict(chrom="category")
    for col in data_cols:
        dtype_[col] = np.uint16 if large_ref else np.uint8

    # which file a column is in
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
                logging.debug(f"found col {col} in header {j}")
                break

    if None in file_ix:
        s = [c for i, c in zip(file_ix, data_cols) if i is None]
        raise ValueError("columns not found in reference: " + ", ".join(s))

    # read correct cols from each file
    for i, ref_file in enumerate(ref_files):

        cols0 = basic_cols + [col for ix, col in zip(file_ix, data_cols) if ix == i]
        if map_file == i:
            cols0 = cols0 + [map_col]
            ix_cols = basic_cols + [map_col]
        else:
            ix_cols = basic_cols

        ref0 = pd.read_csv(ref_file, dtype=dtype_, usecols=cols0)
        ref0.set_index(ix_cols, inplace=True)
        ref0.index.rename("map", level=map_col, inplace=True)

        ref0 = ref0.loc[~ref0.index.duplicated()]
        ref0 = ref0[~np.isnan(ref0.reset_index("map")["map"].values)]

        if i == 0:
            ref = ref0
        else:
            ref = ref.join(ref0)

    ref.fillna(value=0, inplace=True)

    # aggregate different labels
    ref = ref.rename(D, axis=1).T.groupby(level=0).sum().T

    if autosomes_only:
        if type(sex_chroms) is not list:
            sex_chroms = parse_chroms(sex_chroms)
        ref = ref[~ref.index.get_level_values("chrom").isin(sex_chroms)]
    return ref


def filter_ref(
    ref,
    states,
    ancestral=None,
    cont=None,
    filter_delta=None,
    filter_pos=None,
    filter_map=None,
    filter_ancestral=False,
    filter_cont=True,
    **kwargs,
):

    if filter_ancestral and ancestral is not None:
        no_ancestral_call = ref[f"{ancestral}_ref"] + ref[f"{ancestral}_alt"] == 0
        ref = ref.loc[~no_ancestral_call]
        logging.info(
            "filtering %s SNP due to missing ancestral call", np.sum(no_ancestral_call)
        )

    if filter_cont and cont is not None:
        no_cont_data = ref[f"{cont}_ref"] + ref[f"{cont}_alt"] == 0
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


def qbin(x, bin_size=1000, neg_is_nan=True, short_threshold=100_000):
    if bin_size == -1:
        """case when we only want deam in first 3 pos vs everything else"""
        res = pd.Series(x, dtype=np.uint8)
        res.loc[(x >= 0) & (x <= short_threshold)] = 0
        res.loc[x > short_threshold] = 255
        res.loc[x < 0] = 255
        return res
    if neg_is_nan:
        y_short = x.loc[(x >= 0) & (x <= short_threshold)]
        y_long = x.loc[(x > short_threshold)]

        n_short_bins = min((254, max((1, int(y_short.shape[0] // bin_size)))))
        short_cuts = pd.qcut(
            y_short, q=n_short_bins, precision=0, labels=False, duplicates="drop"
        )
        if len(short_cuts) == 0:
            ymax = 1
        else:
            ymax = max(short_cuts)

        n_long_bins = min((254 - ymax, max((1, int(y_long.shape[0] // bin_size)))))
        long_cuts = pd.qcut(
            y_long, q=n_long_bins, precision=0, labels=False, duplicates="drop"
        )

        res = pd.Series(x, dtype=np.uint8)
        res.loc[(x >= 0) & (x <= short_threshold)] = short_cuts
        res.loc[x > short_threshold] = long_cuts + ymax + 1
        res.loc[x < 0] = 255
        return res
    else:
        n_bins = min((254, max((1, int(x.shape[0] // bin_size)))))
        cuts = pd.qcut(
            x, q=n_bins, precision=0, labels=False, duplicates="drop"
        )  # .astype(np.uint8)
        print(n_bins, "x", x.shape, min(Counter(cuts).values()))
        return cuts


def bin_reads(data, deam_bin_size=10000, len_bin_size=1000, short_threshold=2):
    """Bin reads into bins with approximately the same number of reads"""
    data.reset_index(inplace=True)
    data["n"] = data["tref"] + data["talt"]
    data["deam_bin"] = data.groupby(["lib"], group_keys=False).deam.apply(
        qbin, bin_size=deam_bin_size, short_threshold=short_threshold
    )
    data["len_bin"] = data.groupby(["lib", "deam_bin"], group_keys=False).len.apply(
        qbin, bin_size=len_bin_size, neg_is_nan=False
    )
    data["rg"] = [
        f"{lib}_{d}_{l}" for lib, d, l in zip(data.lib, data.deam_bin, data.len_bin)
    ]
    data.set_index(["chrom", "pos"], inplace=True)

    ix = (
        data[["rg", "lib", "deam_bin", "deam", "len_bin", "len"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # df that count reads per bin and for each data atom
    # n_bin = data[['rg','n']].groupby('rg').sum()
    n_exact = data[["rg", "deam", "len", "n"]].groupby(["rg", "deam", "len"]).sum()
    ix = ix.set_index(["rg", "deam", "len"]).join(n_exact, how="left").reset_index()
    ix.rename(columns={"n": "n_exact"}, inplace=True)

    return ix


def load_gt_data(infile):
    dtype_mandatory = dict(
        chrom="category", pos=np.uint32, talt=np.uint8, tref=np.uint8
    )
    data = pd.read_csv(infile, dtype=dtype_mandatory)
    data.set_index(["chrom", "pos"], inplace=True)
    return data


def load_read_data(
    infile,
    split_lib=True,
    downsample=1,
    deam_bin_size=10000,
    len_bin_size=1000,
    high_cov_filter=0.001,
    make_bins=False,
):
    """loading the read data, second generation of input file.
    The idea here (and chief difference to the previous vesion) is that the
    binning of reads/data is done AFTER reading it in.

    Returns:
    data: a pandas df indexed by chromosome and position, with data columns
        - 'rg' (a string identifier of reads grouped together'
        - 'tref' number of ref alleles at site
        - 'talt' number of alt alleles at site

    ix: a pandas df that gives relationship between rg and read properties
        - deam: at which positionr a read has the first deamination (-1 = no
          deamination)
        - len: length (in base pairs) of the read
        - lib: the library ID (string)
        - deam_bin: bin id for deaminated sites (no-deamination is assigned to
          bin 255, max 256 bins)
        - len_bin: bin id for length (max 256 bins)
        - n_bin/n_exact (number of reads in each bin) or each category
    """
    dtype_mandatory = dict(
        chrom="category", pos=np.uint32, talt=np.uint8, tref=np.uint8
    )

    dtype_optional = dict(
        lib=str, rg=str, score=int, deam=np.int16, len=np.uint8, dmgpos=bool
    )

    data0 = pd.read_csv(infile, dtype=dtype_mandatory, nrows=1)
    for c in data0.columns:
        if c in dtype_optional:
            dtype_mandatory[c] = dtype_optional[c]

    data = pd.read_csv(
        infile, dtype=dtype_mandatory, usecols=dtype_mandatory.keys()
    ).dropna()
    data.set_index(["chrom", "pos"], inplace=True)

    if "rg" not in data:
        if "lib" in data:
            data["rg"] = data["lib"]
        else:
            data["rg"] = "lib0"
            data["lib"] = data["rg"]
    elif not split_lib:
        data = data.groupby(data.index.names).agg(sum)

    if downsample < 1:
        data.tref = binom.rvs(data.tref, downsample, size=len(data.tref))
        data.talt = binom.rvs(data.talt, downsample, size=len(data.talt))

    # rm sites with no or extremely high coverage
    data = data[data.tref + data.talt > 0]
    q = np.quantile(data.tref + data.talt, 1 - high_cov_filter)
    to_remove = data.tref + data.talt > q
    logging.info(f"high-cov filter at {q}, removing {sum(to_remove)} sites")
    data = data[~to_remove]

    if make_bins:
        ix = bin_reads(data, deam_bin_size, len_bin_size)
        data = (
            data[["rg", "tref", "talt"]].reset_index().set_index(["chrom", "pos", "rg"])
        )
        data = data.groupby(data.index.names, observed=True).sum()
        return data, ix
    else:
        return data, None


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
