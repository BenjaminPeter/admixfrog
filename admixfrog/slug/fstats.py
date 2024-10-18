from copy import deepcopy
import itertools
import numpy as np
import pandas as pd


def calc_fstats(sfs, pop_list=None, name="XXX"):
    """should calculate a table with all f3/f4 between source states and target"""
    SAMPLE_NAME = name
    freqs, pops = sfs_to_freq(sfs)

    f3s, f4s = [], []
    if pop_list is not None:  # only calc subset of stats
        pops = list(set(pops).intersection(pop_list))
    if "rep" not in freqs:
        freqs["rep"] = 0.0
    if "n_snps" not in freqs:
        freqs["n_snps"] = 1
    if "tau" not in freqs:
        freqs["tau"] = np.nan
    pis = freq_to_pi(freqs, pops, SAMPLE_NAME)

    # f3s
    for X in pops + [f"{SAMPLE_NAME}"]:
        for A, B in itertools.combinations(pops + [f"{SAMPLE_NAME}"], 2):
            if A != X and B != X:
                f3s.append(single_f3(pis, X, A, B))
    f3s = pd.concat(f3s).reset_index()

    # f4s without target
    if len(pops >= 4):
        for A, B in itertools.combinations(pops, 2):
            for C, D in itertools.combinations(pops, 2):
                if len(set((A, B, C, D))) == 4:
                    f4s.append(single_f4(pis, A, B, C, D))

    # f4s with target
    A = f"{SAMPLE_NAME}"
    if len(pops >= 4):
        for B in pops:
            for C, D in itertools.combinations(pops, 2):
                if len(set((A, B, C, D))) == 4:
                    f4s.append(single_f4(pis, A, B, C, D))

    f4s = pd.concat(f4s).reset_index()
    return f3s, f4s, pis


def sfs_to_freq(sfs, ref_freq=False):
    """calculate frequencies for each row in an SFS table

    by default, alt allele frequency is calculated, ref_freq gives inverse
    """
    sfs = deepcopy(sfs)
    pops = [p[:-4] for p in sfs.columns if p.endswith("_alt")]
    for pop in pops:
        f = sfs[f"{pop}_alt"] / (sfs[f"{pop}_alt"] + sfs[f"{pop}_ref"])
        sfs[f"within|{pop}"] = 2 * f * sfs[f"{pop}_ref"]
        sfs[f"within|{pop}"] /= sfs[f"{pop}_alt"] + sfs[f"{pop}_ref"] - 1
        del sfs[f"{pop}_alt"], sfs[f"{pop}_ref"]
        sfs[pop] = 1 - f if ref_freq else f
    return sfs, pops


def f_calc_pi(pop1, pop2):
    def f(df):
        pi = df[pop1] * (1 - df[pop2]) + df[pop2] * (1 - df[pop1])
        nans = np.isnan(pi)
        if np.all(nans):
            return np.nan
        return np.average(pi[~nans], weights=df.n_snps[~nans])

    return f


def f_calc_pi_within(pop):
    def f(df):
        return np.average(df[pop], weights=df["n_snps"])

    return f


def f_jk_sd(stat="f3"):
    def f(df):
        n = df.shape[0]
        m = np.mean(df[stat])
        return np.sqrt((n - 1) / n * np.sum((m - df[stat]) ** 2))

    return f


def freq_to_pi(freqs, pops, SAMPLE_NAME="XXX"):
    df = pd.DataFrame()
    fg = freqs.groupby("rep")

    for p1, p2 in itertools.combinations(pops, 2):
        df[f"{p1}|{p2}"] = fg.apply(f_calc_pi(p1, p2))
        df[f"{p2}|{p1}"] = df[f"{p1}|{p2}"]

    # pw to target and within
    for pop in pops:
        df[f"{SAMPLE_NAME}|{pop}"] = fg.apply(f_calc_pi("tau", pop))
        df[f"{pop}|{SAMPLE_NAME}"] = df[f"{SAMPLE_NAME}|{pop}"]
        df[f"within|{pop}"] = fg.apply(f_calc_pi_within(pop))

    df[f"within|{SAMPLE_NAME}"] = 0

    return df


def single_f3(pis, X, A, B):
    f3 = pis[f"{X}|{A}"] + pis[f"{X}|{B}"] - pis[f"{A}|{B}"] - pis[f"within|{X}"]
    f3.name = "f3"
    f3 = pd.DataFrame(f3)
    f3["X"], f3["A"], f3["B"] = X, A, B
    return f3


def single_f4(pis, A, B, C, D):
    f4 = pis[f"{A}|{D}"] + pis[f"{B}|{C}"] - pis[f"{A}|{C}"] - pis[f"{B}|{D}"]
    f4.name = "f4"
    f4 = pd.DataFrame(f4)
    f4["A"], f4["B"], f4["C"], f4["D"] = A, B, C, D
    return f4


def summarize_f3(df):
    f3s = summarize_f(df, stat="f3", pops=["X", "A", "B"])
    f = f3s[["X", "B", "A", "f3", "sd"]]
    f.columns = f3s.columns
    return pd.concat((f3s, f))


def summarize_f4(df):
    f4s = summarize_f(df, stat="f4", pops=["A", "B", "C", "D"])
    f = f4s[["A", "B", "D", "C", "f4", "sd"]]
    f.columns = f4s.columns
    f["f4"] = -f["f4"]
    return pd.concat((f4s, f))


def summarize_f(df, stat="f3", pops=["X", "A", "B"]):
    fg = df.groupby(pops)
    m = fg[stat].mean()
    m.name = stat
    sd = fg.apply(f_jk_sd(stat))
    sd.name = "sd"
    return pd.concat((m, sd), axis=1).reset_index(drop=False)
