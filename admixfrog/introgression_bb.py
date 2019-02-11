import numpy as np
import pandas as pd
from scipy.special import betaln
import sys
from collections import namedtuple, defaultdict, Counter
from scipy.stats import binom
from scipy.optimize import minimize

try:
    from utils import bins_from_bed, data2probs, init_pars, Pars, load_data, load_ref
    from utils import posterior_table
    from distributions import dbetabinom
    from hmmbb import get_emissions_bb_cy
    from hmm_updates import update_contamination
    from fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions
    from posterior_geno import post_geno_py, update_tau
except (ModuleNotFoundError, ImportError):
    from .utils import bins_from_bed, data2probs, init_pars, Pars, load_data, load_ref
    from .utils import posterior_table
    from .distributions import dbetabinom
    from .hmmbb import get_emissions_bb_cy
    from .hmm_updates import update_contamination
    from .fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions
    from .posterior_geno import post_geno_py, update_tau

np.set_printoptions(suppress=True, precision=4)


def p_reads_given_gt(P, n_obs, c, error):
    read_emissions = np.ones((n_obs, 3))
    for g in range(3):
        p = c * P.P_cont + (1 - c) * g / 2
        p = p * (1 - error) + (1 - p) * error
        read_emissions[:, g] = binom.pmf(P.O, P.N, p)

    return read_emissions


def get_po_given_c_py(c, e, O, N, P_cont, G, pg, IX):

    ll = 0.0
    n_states = pg.shape[1]
    for i in range(IX.n_obs):
        chrom, bin_ = IX.OBS2BIN[i]
        snp = IX.OBS2SNP[i]
        for s in range(n_states):
            for g in range(3):
                p = c * P_cont[i] + (1.0 - c) * g / 2.0
                p = p * (1 - e) + (1 - p) * e
                ll += (
                    G[chrom][bin_, s]
                    * pg[snp, s, g]
                    * (O[i] * log(p) + (N[i] - O[i]) * log(1 - p))
                )
    return ll


def update_emissions(E, P, IX, cont, tau, error, bad_bin_cutoff=1e-200):
    n_homo_states = P.alpha.shape[1]
    n_het_states = int(n_homo_states * (n_homo_states - 1) / 2)
    n_states = n_homo_states + n_het_states
    n_obs = P.O.shape[0]
    n_snps = P.alpha.shape[0]
    c = np.array([cont[l] for l in P.lib])

    read_emissions = p_reads_given_gt(P, n_obs, c, error)
    read_emissions[IX.HAPOBS, 1] = 0.0  # convention is that middle gt is set to zero
    # P(SNP | GT)
    gt_emissions = np.ones((n_snps, 3))
    for i, row in enumerate(IX.OBS2SNP):
        gt_emissions[row] *= read_emissions[i]

    GT = np.zeros((n_snps, n_states, 3)) + 300
    # P(GT | Z)
    for s in range(n_homo_states):
        for g in range(3):
            GT[:, s, g] = dbetabinom(
                np.zeros(n_snps, int) + g,
                np.zeros(n_snps, int) + 2,
                tau[s] * (P.alpha[:, s]),
                tau[s] * (P.beta[:, s]),
            )

    s = n_homo_states
    fa, fb = (P.alpha).T, (P.beta).T
    pp = fa / (fa + fb)
    qq = 1 - pp
    for s1 in range(n_homo_states):
        for s2 in range(s1 + 1, n_homo_states):
            for g in range(3):
                GT[:, s, 0] = qq[s1] * qq[s2]
                GT[:, s, 2] = pp[s1] * pp[s2]
                GT[:, s, 1] = 1 - GT[:, s, 0] - GT[:, s, 2]
            s += 1
    assert np.allclose(np.sum(GT, 2), 1)

    GT[IX.HAPSNP, :, 1] = 0.0  # no het emissions
    GT[IX.HAPSNP, n_homo_states:] = 0.0  # no het hidden state
    for s in range(n_homo_states):
        a, b = P.alpha[IX.HAPSNP, s], P.beta[IX.HAPSNP, s]
        GT[IX.HAPSNP, s, 0] = b / (a + b)
        GT[IX.HAPSNP, s, 2] = a / (a + b)

    snp_emissions = np.sum(GT * gt_emissions[:, np.newaxis, :], 2)

    E[:] = 1  # reset
    for bin_, snp in zip(IX.SNP2BIN, snp_emissions):
        E[bin_] *= snp

    E[IX.HAPBIN, n_homo_states:] = 0

    bad_bins = np.sum(E, 1) < bad_bin_cutoff
    print("bad bins", sum(bad_bins))
    E[bad_bins] = bad_bin_cutoff / E.shape[1]


def bw_bb(
    P, IX, pars, max_iter=1000, ll_tol=1e-1, est_contamination=True, est_tau=True,
    freq_contamination =1, freq_tau = 1,
):

    alpha0, trans_mat, cont, error, tau, gamma_names, sex = pars

    libs = np.unique(P.lib)
    ll = -np.inf
    n_states = len(alpha0)

    # create posterior states, and view for each chromosome
    Z, E = (
        np.zeros((sum(IX.bin_sizes), n_states)),
        np.ones((sum(IX.bin_sizes), n_states)),
    )
    gamma, emissions = [], []
    row0 = 0
    for r in IX.bin_sizes:
        gamma.append(Z[row0 : (row0 + r)])
        emissions.append(E[row0 : (row0 + r)])
        row0 += r

    update_emissions(E, P, IX, cont, tau, error)
    n_seqs = len(emissions)

    for it in range(max_iter):

        alpha, beta, n = fwd_bwd_algorithm(alpha0, emissions, trans_mat, gamma)
        ll, old_ll = np.sum([np.sum(np.log(n_i)) for n_i in n]), ll
        assert np.allclose(np.sum(Z, 1), 1)
        if ll - old_ll < ll_tol:
            pg, _, _ = post_geno_py(P, cont, tau, IX, error)
            break
        tpl = (
            IX.n_reads / 1000,
            IX.n_obs / 1000,
            IX.n_snps / 1000,
            IX.n_bins / 1000,
            it,
            np.mean(np.max(Z, 1) >= 0.95),
            ll,
            ll - old_ll,
        )
        print("[%dk|%dk|%dk|%dk]: iter:%d |p95:%.3f\tLL:%.4f\tΔLL:%.4f" % tpl)

        # update stuff
        trans_mat = update_transitions(trans_mat, alpha, beta, gamma, emissions, n, sex)
        alpha0 = np.linalg.matrix_power(trans_mat, 10000)[0]

        if gamma_names is not None:
            print(*gamma_names, sep="\t")
        print(*["%.3f" % a for a in alpha0], sep="\t")


        cond_cont = (est_contamination and (it % freq_contamination == 0 or it < 3))
        cond_tau = (est_tau and (it % freq_tau == 0 or it < 3))
        if  cond_tau or cond_cont:
            pg, x, y = post_geno_py(P, cont, tau, IX, error)
            assert np.all(pg >= 0)
            assert np.all(pg <= 1)
            assert np.allclose(np.sum(pg, 2), 1)
        if cond_cont:
            delta = update_contamination(cont, error, P, Z, pg, IX, libs)
            if delta < 1e-5: #when we converged, do not update contamination
                est_contamination, cond_contamination = False, False
                print("stopping contamination updates")
        if cond_tau:
            delta = update_tau(tau, Z, pg, P, IX)
            if delta < 1e-5: #when we converged, do not update tau
                est_tau, cond_tau = False, False
                print("stopping tau updates")
        if  cond_tau or cond_cont:
            update_emissions(E, P, IX, cont, tau, error)

    pars = Pars(alpha0, trans_mat, dict(cont), error, tau, gamma_names, sex)
    return Z, pg, pars, ll, emissions


def run_hmm_bb(
    infile,
    ref_file,
    state_ids=("AFR", "VIN", "DEN"),
    cont_id="AFR",
    split_lib=True,
    bin_size=1e4,
    prior=0.5,
    sex=None,
    pos_mode=False,
    autosomes_only=False,
    downsample=1,
    tau0=1,
    e0 = 1e-2,
    c0 = 1e-2,
    **kwargs
):

    bin_size = bin_size if pos_mode else bin_size * 1e-6

    data = load_data(infile, split_lib, downsample)
    ref = load_ref(ref_file, state_ids, cont_id, prior, autosomes_only)
    if pos_mode:
        data.map = data.pos
        ref.map = ref.pos

    # sexing stuff
    if "Y" in data.chrom.values:
        sex = "m"
    if sex is None and "X" in data.chrom.values:
        """guess sex"""
        cov = data.groupby(data.chrom == "X").apply(
            lambda df: np.sum(df.tref + df.talt)
        )
        cov = cov.astype(float)
        cov[True] /= np.sum(ref.chrom == "X")
        cov[False] /= np.sum(ref.chrom != "X")

        if cov[True] / cov[False] < 0.8:
            sex = "m"
            print("guessing sex is male, %.4f/%.4f" % (cov[True], cov[False]))
        else:
            sex = "f"
            print("guessing sex is female, %.4f/%.4f" % (cov[True], cov[False]))

    # merge. This is a bit overkill
    ref = ref.merge(data.iloc[:, :3].drop_duplicates()).drop_duplicates()
    data = data.merge(ref)
    ref = ref.sort_values(["chrom", "map", "pos"])
    data = data.sort_values(["chrom", "map", "pos"])
    bins, IX = bins_from_bed(
        bed=ref.iloc[:, :5], data=data, bin_size=bin_size, pos_mode=pos_mode, sex=sex
    )
    P = data2probs(data, ref, state_ids, cont_id, (prior, prior))
    assert ref.shape[0] == P.alpha.shape[0]

    pars = init_pars(state_ids, sex, tau0, e0, c0)

    print("done loading data")

    Z, G, pars, ll, emissions = bw_bb(P, IX, pars, **kwargs)

    viterbi_path = viterbi(pars, emissions)

    # output formating from here
    V = np.array(pars.gamma_names)[np.hstack(viterbi_path)]
    viterbi_df = pd.Series(V, name="viterbi")
    df = pd.DataFrame(Z, columns=pars.gamma_names)
    CC = Counter(IX.SNP2BIN)
    snp = pd.DataFrame([CC[i] for i in range(len(df))], columns=["n_snps"])
    df = pd.concat((pd.DataFrame(bins), viterbi_df, snp, df), axis=1)

    D = (
        data.groupby(["chrom", "pos", "map"])
        .agg({"tref": sum, "talt": sum})
        .reset_index()
    )
    T = posterior_table(G, Z, IX)
    snp_df = pd.concat((D, T, pd.DataFrame(IX.SNP2BIN, columns=["bin"])), axis=1)

    df_libs = pd.DataFrame(pars.cont.items(), columns=["lib", "cont"])
    rgs, deams = [], []
    for l in df_libs.lib:
        try:
            rg, deam = l.split("_")
        except ValueError:
            rg, deam = l, "NA"
        rgs.append(rg)
        deams.append(deam)
    df_libs["rg"] = rgs
    df_libs["deam"] = deams
    CC = Counter(data.lib)
    snp = pd.DataFrame([CC[l] for l in df_libs.lib], columns=["n_snps"])
    df_libs = pd.concat((df_libs, snp), axis=1)
    df_libs.sort_values("n_snps", ascending=False)

    df_pars = pd.DataFrame(pars.trans_mat, columns=pars.gamma_names)
    df_pars["alpha0"] = pars.alpha0
    df_pars["state"] = pars.gamma_names
    df_pars["tau"] = 0
    df_pars["ll"] = ll
    for i in range(len(pars.tau)):
        df_pars.loc[i, "tau"] = pars.tau[i]

    return df, snp_df, df_libs, df_pars, ll
