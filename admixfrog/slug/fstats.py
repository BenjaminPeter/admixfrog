from copy import deepcopy
import itertools
import numpy as np
import pandas as pd


def calc_fstats(sfs, pop_list=None, name="XXX"):
    """should calculate a table with all f2/f3/f4 between source states and target"""


    freqs, pops = sfs_to_freq(sfs)

    breakpoint()

    f2s, f3s, f4s = [], [], []
    if pop_list is not None:  # only calc subset of stats
        pops = list(set(pops).intersection(pop_list))
    if "rep" not in freqs:
        freqs["rep"] = 0.0
    if "n_snps" not in freqs:
        freqs["n_snps"] = 1
    if "tau" not in freqs:
        freqs["tau"] = np.nan
    pis = freq_to_pi(freqs, pops, name)

    # f2s
    if len(pops) >= 1:
        for A, B in itertools.combinations(pops + [f"{name}"], 2):
            if A != B:
                f2s.append(single_f2(pis, A, B))
    if len(f2s) > 0:
        f2s = pd.concat(f2s).reset_index()
    else:
        f2s = pd.DataFrame(columns=["rep", "f2", "A", "B"])

    # f3s
    if len(pops) >= 2:
        for X in pops + [f"{name}"]:
            for A, B in itertools.combinations(pops + [f"{name}"], 2):
                if A != X and B != X:
                    f3s.append(single_f3(pis, X, A, B))
    if len(f3s) > 0:
        f3s = pd.concat(f3s).reset_index()
    else:
        f3s = pd.DataFrame(columns=["rep", "f3", "X", "A", "B"])

    pops2 = [*pops, f"{name}"]
    if len(pops2) >= 4:
        for A, B in itertools.combinations(pops2, 2):
            for C, D in itertools.combinations(pops2, 2):
                if len(set((A, B, C, D))) == 4:
                    f4s.append(single_f4_permuted(pis, A, B, C, D))

    if len(f4s) > 0:
        f4s = pd.concat(f4s).reset_index()
    else:
        f4s = pd.DataFrame(columns=["rep", "f4", "A", "B", "C", "D"])

    """reformating pis"""
    pis_data = pis.T.reset_index(drop=True)
    pairs = pd.DataFrame(
        pis.T.index.str.split("|").tolist(),
        columns=["pop1", "pop2"],
    )
    pairs.loc[pairs.pop1 == "within", "pop1"] = pairs.pop2[pairs.pop1 == "within"]
    pairs["is_between"] = pairs.pop1 != pairs.pop2
    pis_data.set_index(pairs.pop1, inplace=True)
    pis_data.set_index(pairs.pop2, inplace=True, append=True)
    pis_data.set_index(pairs.is_between, inplace=True, append=True)
    pis = pis_data.melt(ignore_index=False, value_name="pi").reset_index()

    # add all the permutations
    return f2s, f3s, f4s, pis


def sfs_to_freq(sfs):
    """calculate frequencies for each row in an SFS table

    by default, alt allele frequency is calculated, ref_freq gives inverse
    """

    breakpoint()

    sfs = deepcopy(sfs)
    pops = [p[:-4] for p in sfs.columns if p.endswith("_der")]
    if len(pops) == 0:
        raise ValueError("no reference populations")

    for pop in pops:
        f = sfs[f"{pop}_der"] / (sfs[f"{pop}_der"] + sfs[f"{pop}_anc"])
        sfs[f"within|{pop}"] = 2 * f * sfs[f"{pop}_anc"]
        sfs[f"within|{pop}"] /= sfs[f"{pop}_der"] + sfs[f"{pop}_anc"] - 1
        del sfs[f"{pop}_der"], sfs[f"{pop}_anc"]
    return sfs, pops


def f_calc_pi(pop1, pop2):
    def f(df):
        if df.n_snps.sum() == 0:  # no data
            return np.nan
        pi = df[pop1] * (1 - df[pop2]) + df[pop2] * (1 - df[pop1])
        nans = np.isnan(pi)
        if np.all(nans):
            return np.nan
        return np.average(pi[~nans], weights=df.n_snps[~nans])

    return f


def f_calc_pi_within(pop):
    def f(df):
        if df.n_snps.sum() == 0:  # no data
            return np.nan
        pi = np.average(df[f"within|{pop}"], weights=df["n_snps"])
        return pi if not np.isnan(pi) else 0

    return f


def f_jk_sd(stat="f3"):
    def f(df):
        n = df.shape[0]
        m = np.mean(df[stat])
        return np.sqrt((n - 1) / n * np.sum((m - df[stat]) ** 2))

    return f


def freq_to_pi(freqs, pops, name="XXX"):
    df = pd.DataFrame()
    fg = freqs.groupby(["rep", "is_sex_chrom"])

    for p1, p2 in itertools.combinations(pops, 2):
        df[f"{p1}|{p2}"] = fg.apply(f_calc_pi(p1, p2), include_groups=False)
        df[f"{p2}|{p1}"] = df[f"{p1}|{p2}"]

    # pw to target and within
    for pop in pops:
        df[f"{name}|{pop}"] = fg.apply(f_calc_pi("tau", pop), include_groups=False)
        df[f"{pop}|{name}"] = df[f"{name}|{pop}"]
        df[f"within|{pop}"] = fg.apply(f_calc_pi_within(pop), include_groups=False)

    df[f"within|{name}"] = 0

    return df


def single_f2(pis, A, B):
    f2 = 2 * pis[f"{A}|{B}"] - pis[f"within|{A}"] - pis[f"within|{B}"]
    f2.name = "f2"
    f2 = pd.DataFrame(f2)
    f2["A"], f2["B"] = A, B
    return f2


def single_f3(pis, X, A, B):
    f3 = pis[f"{X}|{A}"] + pis[f"{X}|{B}"] - pis[f"{A}|{B}"] - pis[f"within|{X}"]
    f3.name = "f3"
    f3 = pd.DataFrame(f3)
    f3["X"], f3["A"], f3["B"] = X, A, B
    return f3


def f4_rot(df, pops):
    x = df.loc[:, ["f4", *pops]].reset_index().to_numpy()
    return x


def single_f4_permuted(*args):
    f4 = single_f4(*args)
    pops0 = [
        ["A", "B", "C", "D"],
        ["B", "A", "D", "C"],
        ["C", "D", "A", "B"],
        ["D", "C", "B", "A"],
    ]
    pops1 = [
        ["A", "B", "D", "C"],
        ["B", "A", "C", "D"],
        ["C", "D", "B", "A"],
        ["D", "C", "A", "B"],
    ]

    rots0 = [f4_rot(f4, p) for p in pops0]

    # flipped sign combinations
    f4.f4 = -f4.f4
    rots1 = [f4_rot(f4, p) for p in pops1]

    f4_full = pd.DataFrame(
        np.vstack([*rots0, *rots1]), columns=f4.reset_index().columns
    )
    f4_full.set_index(f4.index.names, inplace=True)
    return f4_full


def single_f4(pis, A, B, C, D):
    f4 = pis[f"{A}|{D}"] + pis[f"{B}|{C}"] - pis[f"{A}|{C}"] - pis[f"{B}|{D}"]
    f4.name = "f4"
    f4 = pd.DataFrame(f4)
    f4["A"], f4["B"], f4["C"], f4["D"] = A, B, C, D
    return f4


def summarize_f2(df):
    cols = ["is_sex_chrom", "A", "B", "f2", "sd"]
    if len(df) == 0:  # empty case
        return pd.DataFrame(columns=cols)
    f2s = summarize_f(
        df,
        stat="f2",
        pops=[
            "A",
            "B",
        ],
    )
    return f2s[cols]


def summarize_f3(df):
    cols = ["is_sex_chrom", "X", "A", "B", "f3", "sd"]
    if len(df) == 0:  # empty case
        return pd.DataFrame(columns=cols)
    f3s = summarize_f(df, stat="f3", pops=["X", "A", "B"])
    return f3s[cols]


def summarize_pi(pis):
    ix = ["pop1", "pop2", "is_between"]
    pis = summarize_f(pis, stat="pi", pops=ix)
    return pis


def summarize_f4(df):
    cols = ["is_sex_chrom", "A", "B", "C", "D", "f4", "sd"]
    if len(df) == 0:  # empty case
        return pd.DataFrame(columns=cols)
    f4s = summarize_f(df, stat="f4", pops=["A", "B", "C", "D"])
    return f4s[cols]


def summarize_f(df, stat, pops):
    fg = df.groupby(["is_sex_chrom", *pops])
    m = fg[stat].mean()
    m.name = stat
    sd = fg.apply(f_jk_sd(stat), include_groups=False)
    sd.name = "sd"

    return pd.concat((m, sd), axis=1).reset_index(drop=False)
