import pandas as pd
import numpy as np
import logging
from ..utils.classes import FrogData

dtype_bin = np.dtype(
    [("chrom", "U2"), ("map", float), ("pos", int), ("id", int), ("haploid", bool)]
)


def init_pars(
    states,
    n_rgs,
    homo_ids=None,
    het_ids=None,
    F0=0.001,
    tau0=1,
    e0=1e-2,
    c0=1e-2,
    init_guess=None,
    transition_matrix=None,
    bin_size=1.0,
    **kwargs,
):
    """initialize parameters

    returns a pars object
    """

    n_states, n_hap = states.n_states, states.n_hap

    alpha0 = np.ones(n_states) / n_states
    alpha0_hap = np.ones(n_hap) / n_hap

    if transition_matrix is None:
        trans_mat = np.zeros((n_states, n_states)) + 2e-2
        trans_mat_hap = np.zeros((n_hap, n_hap)) + 2e-2

        np.fill_diagonal(trans_mat, 1 - (n_states - 1) * 2e-2)
        np.fill_diagonal(trans_mat_hap, 1 - (n_hap - 1) * 2e-2)

        if init_guess is not None:
            guess = [i for i, n in enumerate(states.state_names) if n in init_guess]
            logging.info("starting with guess %s " % guess)
            trans_mat[:, guess] = trans_mat[:, guess] + 1
            trans_mat /= np.sum(trans_mat, 1)[:, np.newaxis]
    else:
        trans_mat_hap = pd.read_csv(transition_matrix, header=None).to_numpy()
        trans_mat = trans_mat_hap_to_dip(trans_mat_hap)
        trans_mat_hap = expm(trans_mat_hap * bin_size)
        trans_mat = expm(trans_mat * bin_size)

    return FrogPars(
        alpha0=alpha0,
        alpha0_hap=alpha0_hap,
        trans=trans_mat,
        trans_hap=trans_mat_hap,
        cont=np.zeros(n_rgs) + c0,
        error=np.zeros(n_rgs) + e0,
        F=np.zeros(states.n_homo) + F0,
        tau=np.zeros(states.n_homo) + tau0,
    )


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


def get_prior(prior, snp_df, states, n_snps, cont_prior, ancestral_prior):
    alt_ix = ["%s_alt" % s for s in states]
    ref_ix = ["%s_ref" % s for s in states]
    snp_ix_states = set(alt_ix + ref_ix)
    ca, cb = cont_prior

    if states.contamination is not None:
        cont_ref, cont_alt = (
            f"{states.contamination}_ref",
            f"{states.contamination}_alt",
        )
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

        assert np.all(df.tref.values + df.talt.values < 256)
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


def trans_mat_hap_to_dip(tmat):
    """given a haploid transition rate matrix tmat, returns a diploid transition
    matrix

    assumptions:
     - only one transition at a time
     - "canonical" order of heterozygous states
     - independence between haplotypes
    """
    n = tmat.shape[0]
    n_homo = n
    n_het = int(n * (n - 1) / 2)
    tmat2 = np.zeros((n_homo + n_het, n_homo + n_het))

    # homo -> homo transition
    # for i in range(n_homo):
    #    for j in range(n_homo):
    #        tmat2[i, j] = tmat[i, j] ** 2 #prob both change at once

    # homo -> het transition
    for i in range(n_homo):
        c = n_homo  # state
        for h1 in range(n_homo):  # first het
            for h2 in range(h1 + 1, n_homo):  # second het
                if i == h1:  # only transition second haplotype
                    tmat2[i, c] = tmat[h1, h2]
                    tmat2[c, i] = tmat[h2, h1]
                elif i == h2:  # only transition first haplotype
                    tmat2[i, c] = tmat[h2, h1]
                    tmat2[c, i] = tmat[h1, h2]
                else:  # transition both
                    pass
                    # tmat2[i, c] = tmat[i, h1] * tmat[i, h2] + tmat[i, h2] * tmat[j, h1]
                    # tmat2[c, i] = tmat[h1,i] * tmat[h2, i] + tmat[h2, i] * tmat[h1, j]
                c += 1

    # het -> het transition
    c1 = n_homo
    for i in range(n_homo):  # first het from
        for j in range(i + 1, n_homo):  # second het from
            c2 = n_homo
            for h1 in range(n_homo):  # first het of target
                for h2 in range(h1 + 1, n_homo):  # second het of target
                    if i == h1 and j == h2:  # no transition
                        continue
                    elif i == h1:  # transition second haplotype from j to h2
                        tmat2[c1, c2] = tmat[j, h2]
                    elif j == h2:  # transition first haplotype from i to h1
                        tmat2[c1, c2] = tmat[i, h1]
                    else:  # transition both
                        pass
                    c2 += 1
            c1 += 1

    s = np.sum(tmat2, 1)
    for i in range(tmat2.shape[0]):
        tmat2[i, i] -= s[i]

    return tmat2


def make_bins(snp, chroms, bin_size, n_snps, snp_mode, haplo_chroms, diplo_chroms):
    SNP2BIN = np.empty((n_snps), int)

    bin_loc = []
    bin0 = 0
    n_bins_by_chr = []
    diploid_snps, haploid_snps = [], []

    for i, chrom in enumerate(chroms):
        on_chrom = snp.chrom == chrom
        map_ = snp.map[on_chrom]
        pos = snp.pos[on_chrom]

        chrom_start = float(np.floor(map_.head(1) / bin_size) * bin_size)
        chrom_end = float(np.ceil(map_.tail(1) / bin_size) * bin_size)
        chrom_is_hap = chrom in haplo_chroms

        if snp_mode:
            bins = np.arange(bin0, bin0 + np.sum(on_chrom))
            logging.debug("binning chrom %s: %d snp bins", chrom, len(bins))

            n_bins_by_chr.append(len(bins))
            bin_ids = bins
            _bin = np.empty_like(bins, dtype_bin)
            _bin["chrom"] = chrom
            _bin["pos"] = pos
            _bin["id"] = bin_ids
            _bin["map"] = bins
            _bin["haploid"] = chrom_is_hap
            bin_loc.append(_bin)

            # assign each SNP to its bin
            snp_ids = bins
            SNP2BIN[bins] = bins

        else:
            # create bins
            bins = np.arange(chrom_start, chrom_end, bin_size)
            logging.debug("binning chrom %s: %d bins", chrom, len(bins))

            n_bins_by_chr.append(len(bins))
            bin_ids = range(bin0, bin0 + len(bins))
            _bin = np.empty_like(bins, dtype_bin)
            _bin["chrom"] = chrom
            _bin["pos"] = np.interp(bins, map_, pos)
            _bin["id"] = bin_ids
            _bin["map"] = bins
            _bin["haploid"] = chrom_is_hap
            bin_loc.append(_bin)

            # put SNPs in bins
            snp_ids = snp.snp_id[on_chrom]
            dig_snp = np.digitize(snp.map[on_chrom], bins, right=False) - 1
            SNP2BIN[snp_ids] = dig_snp + bin0

        if chrom_is_hap:
            haploid_snps.extend(snp_ids)
        else:
            diploid_snps.extend(snp_ids)

        bin0 += len(bins)

    bins = pd.DataFrame(np.hstack(bin_loc))
    logging.debug("done creating bins")

    # for now, assume data is ordered such that diploid chroms come before haploid ones
    assert (
        len(haploid_snps) == 0
        or len(diploid_snps) == 0
        or min(haploid_snps) > max(diploid_snps)
    )

    if len(haploid_snps) > 0:
        haploid_snps = slice(min(haploid_snps), max(haploid_snps) + 1)
    else:
        haploid_snps = slice(0, 0)
    if len(diploid_snps) > 0:
        diploid_snps = slice(min(diploid_snps), max(diploid_snps) + 1)
    else:
        diploid_snps = slice(0, 0)
    return bins, n_bins_by_chr, SNP2BIN, diploid_snps, haploid_snps


def data2probs(
    df,
    states,
    bin_size,
    sex="F",
    prior=None,
    cont_id=None,
    cont_prior=(1e-8, 1e-8),
    ancestral_prior=0,
    ancestral=None,
    snp_mode=False,
):
    """create data structure that holds the reference genetic data

    creates an object of type `FrogData` defined in utils/classes.py
    It will contain both the target and reference data, as well as
    data structures that assign each SNP to bins and read groups

    input:
    df: merged reference and SNP data. has columns tref, talt with the read
    counts at each SNP, and "X_alt, X_ref" for each source pop
    states: state object with the required states
    sex : 'm' or 'f', for male or female
    bin_size: size of each bin
    prior: None for empirical bayes prior, otherwise prior to be added
    ancestral: ancestray allele
    ancestral_prior=0,
    cont_id, cont_prior: contamination population and prior
    """

    OBS2SNP = df.index.get_level_values("snp_id")
    n_snps = len(pd.unique(OBS2SNP))

    chroms = pd.unique(df.index.get_level_values("chrom"))
    haplo_chroms, diplo_chroms = [], []

    if sex is None:
        diplo_chroms = chroms
    else:
        for c in chroms:
            if c[0] in "YyWw":
                haplo_chroms.append(c)
            elif c[0] in "Xx" and sex == "m":
                haplo_chroms.append(c)
            elif c[0] in "Zz" and sex == "f":
                haplo_chroms.append(c)
            elif c.startswith("hap"):
                haplo_chroms.append(c)
            else:
                diplo_chroms.append(c)

    snp = df.index.to_frame(index=False).drop_duplicates()
    bins, n_bins_by_chr, SNP2BIN, diploid_snps, haploid_snps = make_bins(
        snp, chroms, bin_size, n_snps, snp_mode, haplo_chroms, diplo_chroms
    )

    rgs = pd.unique(df.rg)
    RG2OBS = dict((l, np.where(df.rg == l)[0]) for l in rgs)

    alt_prior, ref_prior = get_prior(
        prior, df, states, n_snps, cont_id, cont_prior, ancestral, ancestral_prior
    )

    P = FrogData(
        O=np.array(df.talt.values, np.uint8),
        N=np.array(df.tref.values + df.talt.values, np.uint8),
        psi=np.array(0.0)
        if cont_id is None
        else np.array((df[cont_ref] + ca) / (df[cont_ref] + df[cont_alt] + ca + cb)),
        alpha=alt_prior[diploid_snps],
        beta=ref_prior[diploid_snps],
        alpha_hap=alt_prior[haploid_snps],
        beta_hap=ref_prior[haploid_snps],
        states=states,
        n_bins=bins.shape[0],
        n_bins_by_chr=n_bins_by_chr,
        rgs=rgs,
        RG2OBS=RG2OBS,
        OBS2SNP=OBS2SNP,
        SNP2BIN=SNP2BIN,
        chroms=chroms,
        haplo_chroms=haplo_chroms,
        diplo_chroms=diplo_chroms,
        haploid_snps=haploid_snps,
        diploid_snps=diploid_snps,
        sex=sex,
    )
    return P, bins
