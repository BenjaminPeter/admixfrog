import logging
import yaml
from numba import njit
from collections import Counter
from scipy.stats import binom
import pandas as pd
import numpy as np
import itertools

from .utils import posterior_table, posterior_table_slug


def write_cont_table_frog(df, rgs, cont, error, tot_n_snps, outname=None):
    df_libs = pd.DataFrame(zip(rgs, cont, error), columns=["rg", "cont", "error"])

    libs, deams, len_bins = [], [], []
    for l in df_libs.rg:
        try:
            lib, len_bin, deam = l.split("_")
        except ValueError:
            lib, len_bin, deam = l, 0, "NA"
        libs.append(lib)
        deams.append(deam)
        len_bins.append(len_bin)
    df_libs["lib"] = libs
    df_libs["len_bin"] = len_bins
    df_libs["deam"] = deams

    CC = df.groupby(["rg"]).agg(({"tref": sum, "talt": sum})).reset_index()
    CC["n_reads"] = CC.tref + CC.talt
    del CC["tref"]
    del CC["talt"]

    df_libs = df_libs.merge(CC)
    df_libs.sort_values("n_reads", ascending=False)
    df_libs["tot_n_snps"] = tot_n_snps

    if outname is not None:
        df_libs.to_csv(outname, float_format="%.6f", index=False)

    return df_libs


def write_bin_table(Z, bins, viterbi_df, gamma_names, P, outname=None):
    df_bin = pd.DataFrame(Z, columns=gamma_names)
    CC = Counter(P.SNP2BIN)
    snp = pd.DataFrame([CC[i] for i in range(len(df_bin))], columns=["n_snps"])
    df_bin = pd.concat((pd.DataFrame(bins), viterbi_df, snp, df_bin), axis=1)

    if outname is not None:
        df_bin.to_csv(outname, float_format="%.6f", index=False)

    return df_bin


def write_snp_table(data, G, Z, P, gt_mode=False, outname=None):
    D = (
        data.reset_index(drop=False)
        .groupby("snp_id")
        .agg(
            {
                "tref": sum,
                "talt": sum,
                "chrom": lambda x: x.iloc[0],
                "pos": min,
                "map": min,
            }
        )
        .reset_index()
    )

    if gt_mode:
        snp_df = pd.concat((D, pd.DataFrame(P.SNP2BIN, columns=["bin"])), axis=1)
        snp_df = snp_df[["chrom", "pos", "map", "snp_id", "bin", "tref", "talt"]]
    else:
        T = posterior_table(G, Z, P)
        snp_df = pd.concat((D, T, pd.DataFrame(P.SNP2BIN, columns=["bin"])), axis=1)
        snp_df["random_read"] = snp_df["random_read"].astype(int)

    snp_df.sort_values(["chrom", "pos"], inplace=True)
    if outname is not None:
        snp_df.to_csv(outname, float_format="%.6f", index=False)

    return snp_df


def write_est_runs(df, outname=None):
    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False)


def write_sim_runs(df, outname=None):
    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False)
