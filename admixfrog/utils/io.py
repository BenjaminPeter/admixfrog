import logging
import yaml
from numba import njit
from collections import Counter
from scipy.stats import binom
import pandas as pd
import numpy as np
import itertools

from .utils import posterior_table, parse_state_string, posterior_table_slug


"""reading and writing files"""


def load_ref(
    ref_files,
    state_dict,
    cont_id,
    ancestral=None,
    autosomes_only=False,
    map_col="map",
    large_ref=True,
):
    """loads reference in custom (csv) format
    ref_files: paths to files
    state_dict: a dict D[label] : pop. All populatios in sources with label are added to pop. This is 
        expected to be generated with utils.parse_state_string
    """

    # 1. get list of states we care about
    label_states = list(state_dict.keys())
    EXT = ["_ref", "_alt"]
    D = dict(
        ((k + e), (v + e)) for ((k, v), e) in itertools.product(state_dict.items(), EXT)
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
        ref0.chrom = pd.Categorical(ref0.chrom, ordered=True)
        ref0.set_index(ix_cols, inplace=True)
        ref0.index.rename("map", level=map_col, inplace=True)

        # ref0.chrom.cat.reorder_categories(pd.unique(ref0.chrom), inplace=True)
        ref0 = ref0.loc[~ref0.index.duplicated()]
        ref0 = ref0[~np.isnan(ref0.reset_index("map")['map'].values)]

        if i == 0:
            ref = ref0
        else:
            ref = ref.join(ref0)

    ref.fillna(value=0, inplace=True)

    # aggregate different labels
    ref = ref.rename(D, axis=1).groupby(level=0, axis=1).agg(sum)

    if autosomes_only:
        ref = ref[ref.chrom != "X"]
        ref = ref[ref.chrom != "Y"]
        ref = ref[ref.chrom != "mt"]

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
    n_states = len(states)

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


def qbin_old(x, bin_size=1000, neg_is_nan=True, short_threshold=100_000):
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
    data.reset_index(inplace=True)
    data["deam_bin"] = data.groupby(["lib"]).deam.apply(
        qbin, bin_size=deam_bin_size, short_threshold=short_threshold
    )
    data["len_bin"] = data.groupby(["lib", "deam_bin"]).len.apply(
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
    ix = data.groupby("rg").agg(({"tref": sum, "talt": sum})).reset_index().merge(ix)
    ix["n_bin"] = ix.tref + ix.talt
    del ix["tref"], ix["talt"]

    ix = (
        data.groupby(["rg", "deam", "len"])
        .agg(({"tref": sum, "talt": sum}))
        .reset_index()
        .merge(ix)
    )
    ix["n_exact"] = ix.tref + ix.talt
    ix.n_exact = ix.n_exact.astype(int)
    ix.n_bin = ix.n_bin.astype(int)
    return ix


def load_read_data(
    infile,
    split_lib=True,
    downsample=1,
    make_bins=True,
    deam_bin_size=10000,
    len_bin_size=1000,
    high_cov_filter=0.001,
):
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
    data.chrom = pd.Categorical(data.chrom, ordered=True)
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
        return data, ix

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


# yaml writer with indent
class IndentDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(IndentDumper, self).increase_indent(flow, False)


def write_pars_table(pars, outname=None):
    try:
        P = pars.__dict__.copy()
        for k, v in P.items():
            try:
                P[k] = v.tolist()
            except:
                pass

    except AttributeError:
        P = pars
    s = yaml.dump(P, Dumper=IndentDumper, default_flow_style=False, indent=4)

    if outname is not None:
        with open(outname, "wt") as f:
            f.write(s)

    return s


def write_cont_table(df, cont, error, tot_n_snps, outname=None):
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

    CC = df.groupby(["lib"]).agg(({"tref": sum, "talt": sum})).reset_index()
    CC["n_reads"] = CC.tref + CC.talt
    del CC["tref"]
    del CC["talt"]

    df_libs = df_libs.merge(CC)
    df_libs.sort_values("n_reads", ascending=False)
    df_libs["tot_n_snps"] = tot_n_snps

    if outname is not None:
        df_libs.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return df_libs


def write_cont_table_slug(ix, rgs, cont, n_reads, tot_n_snps, se=None, outname=None):
    df_rg = pd.DataFrame([k for k in rgs], columns=["rg"])
    df_cont = pd.DataFrame(cont, columns=["cont"])
    df_n_reads = pd.DataFrame(n_reads, columns=["n_reads"])
    df_libs = pd.concat((df_rg, df_cont, df_n_reads), axis=1)
    df_libs["n_sites"] = tot_n_snps
    if se is not None:
        df_libs["se_cont"] = se
        df_libs["l_cont"] = np.clip(cont - 1.96 * se, 0, 1)
        df_libs["h_cont"] = np.clip(cont + 1.96 * se, 0, 1)

    if ix is None:
        df = df_libs
    else:
        df = ix.merge(df_libs, how="left")
        df.sort_values(["len", "deam", "lib"], inplace=True)

    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return df


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
        data.reset_index(drop=False)
        .groupby("snp_id")
        .agg(
            {
                "tref": sum,
                "talt": sum,
                "chrom" : lambda x: x.iloc[0],
                "pos": min,
                "map": min,
            }
        )
        .reset_index()
    )

    if gt_mode:
        snp_df = pd.concat((D, pd.DataFrame(IX.SNP2BIN, columns=["bin"])), axis=1)
        snp_df = snp_df[['chrom', 'pos', 'map', 'snp_id', 'bin', 'tref', 'talt']]
    else:
        T = posterior_table(G, Z, IX)
        snp_df = pd.concat((D, T, pd.DataFrame(IX.SNP2BIN, columns=["bin"])), axis=1)

    snp_df.sort_values(['chrom', 'pos'], inplace=True)
    if outname is not None:
        snp_df.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return snp_df


def write_snp_table_slug(df, posterior_gt, data, outname=None):
    D = (
        df.groupby(["chrom", "pos", "map", "ref", "alt"])
        .agg({"tref": sum, "talt": sum})
        .reset_index()
    )

    T = posterior_table_slug(pg=posterior_gt, data=data)
    snp_df = pd.concat((D, T, pd.DataFrame(data.SNP2SFS, columns=["sfs"])), axis=1)
    if outname is not None:
        snp_df.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return snp_df


def write_snp_table_slug2(df, posterior_gt, data, outname=None):
    D = (
        df.groupby(["chrom", "pos", "map", "ref", "alt"])
        .agg({"tref": sum, "talt": sum})
        .reset_index()
    )

    T = posterior_table_slug(pg=posterior_gt, data=data)
    snp_df = pd.concat((D, T, pd.DataFrame(data.SNP2SFS, columns=["sfs"])), axis=1)
    if outname is not None:
        snp_df.to_csv(outname, float_format="%.6f", index=False, compression="xz")

    return snp_df


def write_vcf_header():
    s = """##fileformat=VCFv4.2\n"""
    s += '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\n'
    s += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    s += '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
    s += '##FORMAT=<ID=GL,Number=1,Type=Float,Description="Genotype Likelihood">\n'
    s += '##FORMAT=<ID=GP,Number=1,Type=Float,Description="Genotype Probability">\n'
    return s


def write_vcf_line(l, flipped):
    aa = l.alt if flipped else l.ref
    meta = l.chrom, l.pos, l.ref, l.alt, aa
    gts = l.random_read, (l.L0, l.L1, l.L2), (l.G0, l.G1, l.G2), l.tref + l.talt
    CHROM, POS, REF, ALT, ANC = meta
    random_read, (l0, l1, l2), (g0, g1, g2), (depth) = gts
    s = f"{CHROM}\t{POS}\t{CHROM}_{POS}\t{REF}\t{ALT}\t.\t.\tAA={ANC}\tGT:GL:GP:DP\t"
    s += f"{random_read}:{l0:.4f},{l1:.4f},{l2:.4f}:{g0:.4f},{g1:.4f},{g2:.4f}:{int(depth)}"
    s += "\n"
    return s


def write_vcf_chroms(chroms):
    s = ""
    for chrom in chroms:
        s += f"##contig=<ID={chrom}>\n"
    return s


def write_vcf(df, data, posterior_gt, genotype_ll, sample_name="test", outname=None):
    D = (
        df.groupby(["chrom", "pos", "map", "ref", "alt"])
        .agg({"tref": sum, "talt": sum})
        .reset_index()
    )
    T = posterior_table_slug(pg=posterior_gt, data=data, gtll=genotype_ll)
    snp_df = pd.concat((D, T, pd.DataFrame(data.SNP2SFS, columns=["sfs"])), axis=1)

    with open(outname, "wt") as f:
        f.write(write_vcf_header())
        f.write(write_vcf_chroms(data.chroms))
        f.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER\tINFO\tFORMAT\t")
        f.write(f"{sample_name}\n")
        for flipped, (_, l) in zip(data.FLIPPED, snp_df.iterrows()):
            f.write(write_vcf_line(l, flipped))


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


def write_sfs(sfs, pars, data, outname=None):
    df2 = pd.DataFrame(sfs.keys(), columns=data.states)
    index = pd.Series(list(sfs.values()), name="sfs")
    df2 = df2.reindex(sorted(df2.columns), axis=1)
    df2 = pd.concat((index, df2), 1)
    n_snps = pd.DataFrame(pd.Series(dict(Counter(data.SNP2SFS))), columns=["n_snps"])

    n_reads, n_endo = Counter(), Counter()

    for (sfs_, rg, n) in zip(data.OBS2SFS, data.OBS2RG, data.N):
        n_reads[sfs_] += n
        n_endo[sfs_] += n * (1 - pars.cont[rg])

    n_reads = pd.Series((n_reads[i] for i in index), dtype=int, name="n_reads")
    n_endo = pd.Series((n_endo[i] for i in index), dtype=float, name="n_endo")
    F = pd.DataFrame(pars.F, columns=["F"])
    tau = pd.DataFrame(pars.tau, columns=["tau"])

    sfs_df = pd.concat((df2, F, tau, n_snps, n_reads, n_endo), axis=1)

    if outname is not None:
        sfs_df.to_csv(outname, float_format="%5f", index=False, compression="xz")


def write_sfs2(sfs, pars, data, se_tau=None, se_F=None, outname=None):
    n_snps = pd.DataFrame(pd.Series(dict(Counter(data.SNP2SFS))), columns=["n_snps"])

    n_reads, n_endo = Counter(), Counter()

    for (sfs_, rg) in zip(data.READ2SFS, data.READ2RG):
        n_reads[sfs_] += 1
        n_endo[sfs_] += 1 * (1 - pars.cont[rg])

    n_anc, n_der = Counter(), Counter()
    for (sfs_, read_, flipped) in zip(
        data.READ2SFS, data.READS, data.FLIPPED[data.READ2SNP]
    ):
        if flipped:
            n_anc[sfs_] += read_
            n_der[sfs_] += 1 - read_
        else:
            n_anc[sfs_] += 1 - read_  # 0 = anc -> 1-0 = 1 anc allele
            n_der[sfs_] += read_  # normal means 1==derived

    n_reads = pd.Series((n_reads[i] for i in sfs.index), dtype=int, name="n_reads")
    n_endo = pd.Series((n_endo[i] for i in sfs.index), dtype=float, name="n_endo")
    n_anc = pd.Series((n_anc[i] for i in sfs.index), dtype=float, name="n_anc")
    n_der = pd.Series((n_der[i] for i in sfs.index), dtype=float, name="n_der")
    F = pd.DataFrame(pars.F, columns=["F"])
    tau = pd.DataFrame(pars.tau, columns=["tau"])

    sfs_df = pd.concat((sfs, F, tau, n_snps, n_reads, n_endo), axis=1)
    np.nan_to_num(sfs_df.n_snps,copy=False, nan=0)
    sfs_df["read_ratio"] = n_der / (n_anc + n_der + 1e-400)
    sfs_df["cont_est"] = 1 - sfs_df["n_endo"] / sfs_df["n_reads"]
    sfs_df["psi"] = sfs_df["tau"] + (sfs_df["read_ratio"] - sfs_df["tau"]) / (
        sfs_df["cont_est"] + 1e-400
    )

    if se_tau is not None:
        tau = pars.tau
        sfs_df["se_tau"] = se_tau
        sfs_df["l_tau"] = np.clip(tau - 1.96 * se_tau, 0, 1)
        sfs_df["h_tau"] = np.clip(tau + 1.96 * se_tau, 0, 1)

    if se_F is not None:
        F = pars.F
        sfs_df["se_F"] = se_F
        sfs_df["l_F"] = np.clip(F - 1.96 * se_F, 0, 1)
        sfs_df["h_F"] = np.clip(F + 1.96 * se_F, 0, 1)

    if outname is not None:
        sfs_df.to_csv(outname, float_format="%5f", index=False, compression="xz")

    return sfs_df


def write_f3_table(df, outname=None):
    df = df[['X', 'A', 'B', 'f3', 'rep']] 
    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False, compression="xz")
    return df

def write_f4_table(df, outname=None):
    df = df[['A', 'B', 'C', 'D', 'f4', 'rep']] 
    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False, compression="xz")
    return df



