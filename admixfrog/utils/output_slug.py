import logging
import yaml
from numba import njit
from collections import Counter
from scipy.stats import binom
import pandas as pd
import numpy as np
import itertools

from .utils import posterior_table, posterior_table_slug


"""reading and writing files"""


def write_cont_table_slug(ix, rgs, cont, tot_n_snps, se=None, outname=None):
    df = np.empty_like(cont, dtype=object)
    for rg, i in rgs.items():
        df[i] = rg
    df = pd.DataFrame(df, columns=["rg"])

    df["cont"] = cont
    df = df.set_index("rg").join(ix[["rg", "n_exact"]].groupby("rg").sum())
    df.rename({"n_exact": "n"}, inplace=True)
    df["n_sites"] = tot_n_snps

    if se is not None:
        df["se_cont"] = se
        df["l_cont"] = np.clip(cont - 1.96 * se, 0, 1)
        df["h_cont"] = np.clip(cont + 1.96 * se, 0, 1)

    if outname is not None:
        df.reset_index().to_csv(
            outname, float_format="%.6f", index=False
        )

    return df


def write_snp_table_slug(df, posterior_gt, data, outname=None):
    D = (
        df.groupby(["chrom", "pos", "map", "ref", "alt"])
        .agg({"tref": sum, "talt": sum})
        .reset_index()
    )

    T = posterior_table_slug(pg=posterior_gt, data=data)
    snp_df = pd.concat((D, T, pd.DataFrame(data.SNP2SFS, columns=["sfs"])), axis=1)
    if outname is not None:
        snp_df.to_csv(outname, float_format="%.6f", index=False)

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
    np.nan_to_num(sfs_df.n_snps, copy=False, nan=0)
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
        sfs_df.to_csv(outname, float_format="%5f", index=False)

    return sfs_df


def write_f3_table(df, outname=None):
    df = df[["X", "A", "B", "f3", "rep"]]
    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False)
    return df


def write_f4_table(df, outname=None):
    df = df[["A", "B", "C", "D", "f4", "rep"]]
    if outname is not None:
        df.to_csv(outname, float_format="%.6f", index=False)
    return df
