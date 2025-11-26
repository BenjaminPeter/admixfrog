"""
read Reich geno format
"""

import numpy as np
from math import ceil
import itertools
import pandas as pd
from collections.abc import Mapping
from pprint import pprint
from .input import filter_ref
import logging
import warnings


def row_length(n_ind):
    return max(ceil(n_ind / 4), 48)


def read_geno(fname, pops=None, target_ind=None, guess_ploidy=True):
    """read binary Reich-lab format data into a pandas data frame
    See https://github.com/DReichLab/AdmixTools/ for specification

    in general, because pandas is limited with recognized data types,
    this might be somewhat memory inefficient.

    MISSING DATA WARNING:
        - the used data type does not allow for missing data.
        - missing data is thus stored as '3', and will need to
            be handled

    fname: reich-lab format without extensions: requires
            {fname}.ind/.geno/.snp to be present
    pops: filter for populations to be used
        supported are:
            - str  :  single pop / state string
            - list(str) : subset of state strings
            - dict[str1] : [str2, str3,...]  : decoded state string

        populations requested but not present will raise error
    target_ind : for admixfrog stuff, single out an individual as a target
    guess_ploidy: if set to true, tries to estimate ploidy of sample from data
        i.e. if there is a value equal to 2, assume diploid
        if false, set ploidy to diploid always
    """

    with open(f"{fname}.geno", "rb") as f:
        head = f.read(48).strip(b"\x00").split()
        assert head[0] == b"GENO"
        n_ind, n_snp = int(head[1]), int(head[2])
        hash_ind, hash_snp = head[3:]
        rlen = row_length(n_ind)

    z = np.fromfile(f"{fname}.geno", dtype=np.uint8).reshape(n_snp + 1, rlen)[1:]
    X = np.unpackbits(z, axis=1).reshape(n_snp, rlen * 4, 2)
    del z

    ind = pd.read_csv(
        f"{fname}.ind", sep=" ", names=["id", "sex", "pop"], skipinitialspace=True
    )

    snp = pd.read_csv(
        f"{fname}.snp",
        sep=" ",
        skipinitialspace=True,
        names=["snp", "chrom", "map", "pos", "ref", "alt"],
    )

    # chromosome formatting
    snp.chrom = snp.chrom.astype(str)
    snp.loc[snp.chrom == "23", "chrom"] = "X"
    snp.loc[snp.chrom == "24", "chrom"] = "Y"
    snp.loc[snp.chrom == "90", "chrom"] = "mt"
    snp["map"] *= 100  # from Morgan to cM

    ix = pd.MultiIndex.from_frame(ind[["pop", "sex", "id"]])

    Y = pd.DataFrame(X[:, :n_ind, 0] + 2 * X[:, :n_ind, 1], columns=ix)
    Y.index = pd.MultiIndex.from_frame(snp[["chrom", "pos", "map", "ref", "alt"]])
    del X

    # set filtering for populations
    if pops is None:
        pops = {d: d for d in Y.columns.unique(level="pop")}
    logging.debug(pprint(pops))
    assert isinstance(pops, Mapping)

    if target_ind is not None:
        sex = Y.xs(target_ind, level="id", axis=1).columns.get_level_values("sex")[0]
        Y["TARGET", sex, target_ind] = Y.xs(target_ind, level="id", axis=1)
        target_has_data = ((Y["TARGET"] <= 2) & (Y["TARGET"] >= 0)).values
        Y = Y[target_has_data]
        pops["TARGET"] = "TARGET"

    Y.rename(pops, axis=1, inplace=True)

    if pops is not None:
        Y = Y[list(set(pops.values()))]

    if guess_ploidy:
        ploidy = Y.agg(lambda x: np.max(x * (x < 3)))
        Y = Y.T.set_index(pd.Index(ploidy, name="ploidy"), append=True).T
    else:
        Y = Y.T.set_index(pd.Index([2] * Y.shape[1], name="ploidy"), append=True).T

    return Y


def ref_alt(Y, copy=False):
    """create a ref/alt table from a geno data-file Y

    Notes:
        - still a bit experimental
        - will flip alleles inplace (i.e. input array Y WILL BE CHANGED
            unless copy=True is set)
        - recognizes `Y` and `mt` as haploid
        - recognizes `X` as haploid if sex is male
        - the population named TARGET is protected
    """
    if copy:
        Y = Y.copy()

    if "Y" in Y.index:
        Y.loc["Y"] = Y.loc["Y"].transform(lambda x: np.where(x > 1, 3, x)).values
    if "mt" in Y.index:
        Y.loc["mt"] = Y.loc["mt"].transform(lambda x: np.where(x > 1, 3, x)).values

    Y[Y > 2] = np.nan
    with warnings.catch_warnings():
        # TODO: pandas axis=1 is deprecated, but alternative is very, very slow
        warnings.filterwarnings("ignore", category=FutureWarning)
        ref = Y.groupby("pop", axis=1).sum()
    ref.rename(columns=lambda x: f"{x}_ref", inplace=True)

    CHROMS = pd.unique(Y.index.get_level_values("chrom"))
    AUTOSOMES = [c for c in CHROMS if c not in ["mt", "Y", "X"]]

    # now transform Y s.t. we get non-ref counts
    Y.loc[AUTOSOMES] = Y.loc[AUTOSOMES].transform(lambda x: x.name[3] - x).values

    if "X" in Y.index:
        Y.loc["X"] = (
            Y.loc["X"]
            .transform(lambda x: (x.name[3] if x.name[1] != "M" else 1) - x)
            .values
        )

    if "Y" in Y.index:  #  always haploid
        Y.loc["Y"] = Y.loc["Y"].transform(lambda x: 1 - x).values

    if "mt" in Y.index:  # always haploid
        Y.loc["mt"] = Y.loc["mt"].transform(lambda x: 1 - x).values

    with warnings.catch_warnings():
        # TODO: pandas axis=1 is deprecated, but alternative is very, very slow
        warnings.filterwarnings("ignore", category=FutureWarning)
        alt = Y.groupby("pop", axis=1).sum()
    alt.rename(columns=lambda x: f"{x}_alt", inplace=True)

    df = ref.merge(alt, on=["chrom", "pos", "map", "ref", "alt"])
    df.rename(columns={"TARGET_ref": "tref", "TARGET_alt": "talt"}, inplace=True)

    return df


def ref_count(x):
    """count number of reference alleles"""
    v = np.sum(x * (x < 3), axis=0)
    return v


def read_geno_ref(*args, **kwargs):
    Y = read_geno(*args, **kwargs)
    return ref_alt(Y)
