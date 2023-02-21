import logging
import yaml
from collections import Counter, defaultdict
from scipy.stats import binom
import pandas as pd
import numpy as np
import itertools
from ..utils.geno_io import read_geno_ref, read_geno
from .utils import posterior_table, posterior_table_slug
from .utils import guess_sex, filter_ref
from .classes import FrogData


def load_ref(
    ref_files,
    states,
    autosomes_only=False,
    map_col="map",
    large_ref=True,
):
    """loads reference in custom (csv) format
    ref_files: paths to files
    """

    cont_id, ancestral = states.contamination, states.ancestral

    # 1. get list of states we care about
    label_states = list(states.state_dict.keys())
    EXT = ["_ref", "_alt"]
    D = dict(
        ((k + e), (v + e))
        for ((k, v), e) in itertools.product(states.state_dict.items(), EXT)
    )

    if states.ancestral is not None:
        label_states = list(set(list(label_states) + [states.ancestral]))
    if states.contamination is not None:
        label_states = list(set(list(label_states) + [states.contamination]))

    # 2. required in every ref
    basic_cols = ["chrom", "pos", "ref", "alt"]  # required in every ref

    # target states
    ref_cols = [f"{s}_ref" for s in label_states]
    alt_cols = [f"{s}_alt" for s in label_states]
    data_cols = ref_cols + alt_cols

    dtype_ = dict(chrom="object")
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
        ref0 = ref0[~np.isnan(ref0.reset_index("map").map.values)]

        if i == 0:
            ref = ref0
        else:
            ref = ref.join(ref0)

    ref.fillna(value=0, inplace=True)

    # aggregate different labels
    ref = ref.rename(D, axis=1).groupby(level=0, axis=1).agg(sum)

    if autosomes_only:
        ref = ref[ref.index.get_level_values("chrom") != "X"]
        ref = ref[ref.index.get_level_values("chrom") != "Y"]
        ref = ref[ref.index.get_level_values("chrom") != "mt"]

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
    dtype_mandatory = dict(chrom="object", pos=np.uint32, talt=np.uint8, tref=np.uint8)

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


def load_admixfrog_data_geno(
    geno_file,
    states,
    ancestral,
    filter=filter,
    target_ind=None,
    guess_ploidy=False,
    pos_mode=False,
):
    df = read_geno_ref(
        fname=geno_file,
        pops=states.state_dict,
        target_ind=target_ind,
        guess_ploidy=guess_ploidy,
    )
    df = filter_ref(df, states, **filter)
    df["rg"] = "rg0"
    if pos_mode:
        df.reset_index("map", inplace=True)
        df.map = df.index.get_level_values("pos")
        df.set_index("map", append=True, inplace=True)

    cats = pd.unique(df.index.get_level_values("chrom"))
    chrom_dtype = pd.CategoricalDtype(cats, ordered=True)

    if not df.index.is_unique:
        dups = df.index.duplicated()
        logging.warning(f"\033[91mWARNING: {np.sum(dups)} duplicate sites found\033[0m")
        logging.warning(" ==> strongly consider re-filtering input file")
        df = df[~dups]
    return df


def load_admixfrog_data(
    states,
    target=None,
    target_file=None,
    gt_mode=False,
    filter=defaultdict(lambda: None),
    ref_files=None,
    geno_file=None,
    split_lib=True,
    pos_mode=False,
    downsample=1,
    guess_ploidy=True,
    sex=None,
    map_col="map",
    fake_contamination=0,
    autosomes_only=False,
    bin_reads=False,
    deam_bin_size=100000,
    len_bin_size=10000,
    prior=None,
    cont_prior=(1e-8, 1e-8),
    ancestral_prior=0,
):
    """
    we have the following possible input files
    1. only a geno file (for target and reference)
    2. custom ref/target (from standalone csv)
    3. geno ref (and target from csv)
    """
    tot_n_snps = 0
    cont_ref, cont_alt = f'{states.contamination}_ref', f'{states.contamination}_alt', 

    "1. only geno file"
    if ref_files is None and target_file is None and geno_file and target:
        df = load_admixfrog_data_geno(
            geno_file=geno_file,
            states=states,
            ancestral=ancestral,
            filter=filter,
            target_ind=target,
            guess_ploidy=guess_ploidy,
            pos_mode=pos_mode,
        )
        ix = None

    elif ref_files and target_file and (geno_file is None) and (target is None):
        "2. standard input"
        hcf = filter.pop("filter_high_cov")

        # load reference first
        ref = load_ref(ref_files, states, autosomes_only, map_col=map_col)
        ref = filter_ref(ref, states, **filter)

        if gt_mode:  # gt mode does not do read emissions, assumes genotypes are known
            data, ix = load_read_data(target_file, make_bins=False)
            assert np.max(data.tref + data.talt) <= 2
            assert np.min(data.tref + data.talt) >= 0
            if not data.index.is_unique:
                dups = data.index.duplicated()
                logging.warning(
                    f"\033[91mWARNING: {np.sum(dups)} duplicate sites found\033[0m"
                )
                logging.warning(" ==> strongly consider re-filtering input file")
                data = data[~dups]
        else:
            data, ix = load_read_data(
                target_file,
                split_lib,
                downsample,
                make_bins=bin_reads,
                len_bin_size=len_bin_size,
                deam_bin_size=deam_bin_size,
                high_cov_filter=hcf,
            )
            # data = data[["rg", "tref", "talt"]]

        # sexing stuff
        if sex is None:
            sex = guess_sex(ref, data)

        if fake_contamination and samples.contamination:
            """filter for SNP with fake ref data"""
            ref = ref[ref[cont_ref] + ref[cont_alt] > 0]

        tot_n_snps = ref.shape[0]

        if pos_mode:
            ref.reset_index("map", inplace=True)
            ref.map = ref.index.get_level_values("pos")
            ref.set_index("map", append=True, inplace=True)
        ref = ref.loc[~ref.index.duplicated()]

        df = ref.join(data, how="inner")

    elif geno_file and target_file and ref_file is None:
        "4. geno ref, standard target"
        raise NotImplementedError("")
    else:
        raise NotImplementedError("ambiguous input")

    # sort by chromosome name. PANDAS is VERY DUMB WTF
    cats = pd.unique(df.index.get_level_values("chrom"))
    chrom_dtype = pd.CategoricalDtype(cats, ordered=True)

    df.index = df.index.set_levels(df.index.levels[0].astype(chrom_dtype), level=0)
    df.reset_index(inplace=True)
    df.sort_values(["chrom", "pos"], inplace=True)
    df.set_index(["chrom", "pos", "ref", "alt", "map"], inplace=True)

    # get ids of unique snps
    snp_ids = df[~df.index.duplicated()].groupby(df.index.names, observed=True).ngroup()
    snp_ids.rename("snp_id", inplace=True)

    snp_ids = pd.DataFrame(snp_ids)
    snp_ids.set_index("snp_id", append=True, inplace=True)


    df = snp_ids.join(df)

    if fake_contamination and samples.contamination:
        mean_endo_cov = np.mean(df.tref + df.talt)
        """
            C = x / (e+x);
            Ce + Cx = x
            Ce = x - Cx
            Ce = x(1-C)
            Ce / ( 1- C) = x
        """

        prop_cont = fake_contamination
        target_cont_cov = prop_cont * mean_endo_cov / (1 - prop_cont)
        f_cont = df[cont_alt] / (df[cont_ref] + df[cont_alt])

        logging.debug(f"endogenous cov: {mean_endo_cov}")
        logging.debug(f"fake contamination cov: {target_cont_cov}")

        c_ref = np.random.poisson((1 - f_cont) * target_cont_cov)
        c_alt = np.random.poisson(f_cont * target_cont_cov)
        logging.debug(f"Added cont. reads with ref allele: {np.sum(c_ref)}")
        logging.debug(f"Added cont. reads with alt allele: {np.sum(c_alt)}")
        df.tref += c_ref
        df.talt += c_alt


    alt_prior, ref_prior = get_prior(prior, df[~df.index.duplicated()], states,  
                                     cont_prior, ancestral_prior)

    make_bins(snp = snp_ids.index.to_frame(index=False),
              chroms = snp_ids.index.get_level_values('chrom'),
              snp_mode,


        snp, chroms, bin_size,
                                                          n_snps,
                                                          snp_mode,
                                                          haplo_chroms,
                                                          diplo_chroms)

    breakpoint()
    P = FrogData(
        O=np.array(df.talt.values, np.uint8),
        N=np.array(df.tref.values + df.talt.values, np.uint8),
        psi=np.array(0.0)
        if states.contamination is None
        else np.array((df[cont_ref] + ca) / (df[cont_ref] + df[cont_alt] + ca + cb)),
        alt_prior = alt_prior,
        ref_prior = ref_prior,
        snp
        states=states,
        bin_size = bin_size,
        rgs=rgs,
        sex=sex,
    )


    return df, ix, sex, tot_n_snps

def empirical_bayes_prior(der, anc, known_anc=False):
    """using beta-binomial plug-in estimator"""

    n = anc + der
    f = np.nanmean(der / n) if known_anc else 0.5

    H = f * (1.0 - f)  # alt formulation
    if H == 0.0:
        return 1e-6, 1e-6

    V = np.nanvar(der / n) if known_anc else np.nanvar(np.hstack((der / n, anc / n)))

    ab = (H - V) / (V - H / np.nanmean(n))
    if np.nanmean(n) < ab:
        return 1e-6, 1e-6
    pa = max((f * ab, 1e-5))
    pb = max(((1 - f) * ab, 1e-5))
    return pa, pb

def get_prior(prior, snp_df, states, cont_prior, 
              ancestral_prior):
    alt_ix = ["%s_alt" % s for s in states]
    ref_ix = ["%s_ref" % s for s in states]
    snp_ix_states = set(alt_ix + ref_ix)
    ca, cb = cont_prior
    n_snps = snp_df.shape[0]


    if states.contamination is not None:
        cont_ref, cont_alt = f"{states.contamination}_ref",f"{states.contamination}_alt" 
        snp_ix_states.update([cont_ref, cont_alt])
    if states.ancestral is not None:
        anc_ref, anc_alt = f"{states.ancestral}_ref", f"{states.ancestral}_alt"
        snp_ix_states.update([anc_ref, anc_alt])
    if prior is None:  # empirical bayes, estimate from data
        alt_prior = np.empty((n_snps, states.n_raw_states))
        ref_prior = np.empty((n_snps, states.n_raw_states))

        if states.contamination is not None:
            ca, cb = empirical_bayes_prior(snp_df[cont_ref], snp_df[cont_alt])

        if states.ancestral is None:
            for i, (a, b, s) in enumerate(zip(alt_ix, ref_ix, states)):
                pa, pb = empirical_bayes_prior(snp_df[a], snp_df[b])
                logging.info("[%s]EB prior [a=%.4f, b=%.4f]: " % (s, pa, pb))
                alt_prior[:, i] = snp_df[a] + pa
                ref_prior[:, i] = snp_df[b] + pb
        else:
            anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"

            # set up vectors stating which allele is ancestral
            ref_is_anc = (snp_df[anc_ref] > 0) & (snp_df[anc_alt] == 0)
            alt_is_anc = (snp_df[anc_alt] > 0) & (snp_df[anc_ref] == 0)
            ref_is_der, alt_is_der = alt_is_anc, ref_is_anc
            anc_is_unknown = (1 - alt_is_anc) * (1 - ref_is_anc) == 1

            for i, (alt_col, ref_col, s) in enumerate(zip(alt_ix, ref_ix, states)):

                # 1. set up base entries based on observed counts
                alt_prior[:, i] = snp_df[alt_col]
                ref_prior[:, i] = snp_df[ref_col]

                # 2. where anc is unknown, add symmetric prior estimated from data
                pa, pb = empirical_bayes_prior(snp_df[alt_col], snp_df[ref_col])
                logging.info("[%s]EB prior0 [anc=%.4f, der=%.4f]: " % (s, pa, pb))
                alt_prior[anc_is_unknown, i] += pa
                ref_prior[anc_is_unknown, i] += pb

                # 3. where anc is known, create indices
                m_anc = np.array(pd.concat((ref_is_anc, alt_is_anc), axis=1))
                m_der = np.array(pd.concat((ref_is_der, alt_is_der), axis=1))
                ANC = np.array(snp_df[[ref_col, alt_col]])[m_anc]
                DER = np.array(snp_df[[ref_col, alt_col]])[m_der]

                pder, panc = empirical_bayes_prior(DER, ANC, known_anc=True)
                panc += ancestral_prior
                logging.info("[%s]EB prior1 [anc=%.4f, der=%.4f]: " % (s, panc, pder))
                alt_prior[alt_is_anc, i] += panc
                alt_prior[alt_is_der, i] += pder
                ref_prior[ref_is_anc, i] += panc
                ref_prior[ref_is_der, i] += pder

    else:
        """ancestral allele contribution to prior
        the ancestral allele adds one pseudocount to the data
        """
        if ancestral is None:
            prior_anc_alt, prior_anc_ref = np.zeros(1), np.zeros(1)
        else:
            anc_ref, anc_alt = f"{ancestral}_ref", f"{ancestral}_alt"
            prior_anc_alt = snp_df[anc_alt] * ancestral_prior
            prior_anc_ref = snp_df[anc_ref] * ancestral_prior

        alt_prior = (
            snp_df[alt_ix].to_numpy() + np.array(prior_anc_alt)[:, np.newaxis] + prior
        )
        ref_prior = (
            snp_df[ref_ix].to_numpy() + np.array(prior_anc_ref)[:, np.newaxis] + prior
        )
        assert np.all(df.tref.values + df.talt.values < 256)

    return ref_prior, alt_prior


