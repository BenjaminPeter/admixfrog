import numpy as np
import pickle
import pandas as pd
from collections import defaultdict, Counter
import pdb
from .utils import bins_from_bed, data2probs, init_pars, Pars
from .utils import posterior_table, load_data, load_ref
from .read_emissions import update_contamination
from .fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions
from .genotype_emissions import update_post_geno, update_F, update_snp_prob
from .genotype_emissions import update_emissions
from .rle import get_rle
from .decode import pred_sims

np.set_printoptions(suppress=True, precision=4)


def bw_bb(
    P,
    IX,
    pars,
    max_iter=1000,
    ll_tol=1e-1,
    est_contamination=True,
    est_F=True,
    freq_contamination=1,
    freq_F=1,
    haploid=False
):
    n_gt = 3

    alpha0, trans_mat, cont, error, F, gamma_names, sex = pars

    libs = np.unique(P.lib)
    ll = -np.inf
    n_states = len(alpha0)

    # create arrays for posterior, emissions
    Z = np.zeros((sum(IX.bin_sizes), n_states))  # P(Z | O)
    E = np.ones((sum(IX.bin_sizes), n_states))  # P(O | Z)
    # GT = np.ones((IX.n_snps, n_states, n_gt))      # P(G | Z)
    # OBS = np.ones((IX.n_snps, n_gt))               # P(O | G)
    SNP = np.zeros((IX.n_snps, n_states, n_gt))  # P(O, G | Z)
    PG = np.zeros((IX.n_snps, n_states, n_gt))  # P(G Z | O)

    gamma, emissions = [], []
    row0 = 0
    for r in IX.bin_sizes:
        gamma.append(Z[row0 : (row0 + r)])
        emissions.append(E[row0 : (row0 + r)])
        row0 += r

    update_snp_prob(SNP, P, IX, cont, error, F, haploid)  # P(O, G | Z)
    e_scaling = update_emissions(E, SNP, P, IX, haploid)  # P(O | Z)

    for it in range(max_iter):

        alpha, beta, n = fwd_bwd_algorithm(alpha0, emissions, trans_mat, gamma)
        if sex == "m":
            print("male ll doens't take x into account", end="\t")
            ll, old_ll = (
                np.sum([np.sum(np.log(n_i)) for _, n_i in zip(range(22), n)])
                + e_scaling,
                ll,
            )
        else:
            ll, old_ll = np.sum([np.sum(np.log(n_i)) for n_i in n]) + e_scaling, ll
        assert np.allclose(np.sum(Z, 1), 1)
        if np.isnan(ll):
            pdb.set_trace()
        assert not np.isnan(ll)
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
        if ll - old_ll < ll_tol:
            break

        if np.any(np.isnan(Z)):
            raise ValueError("nan observed in state posterior")
        if np.any(np.isnan(E)):
            raise ValueError("nan observed in emissions")

        # update stuff
        trans_mat = update_transitions(trans_mat, alpha, beta, gamma, emissions,
                                       n, sex, haploid)
        alpha0 = np.linalg.matrix_power(trans_mat, 10000)[0]

        if gamma_names is not None:
            print(*gamma_names, sep="\t")
        print(*["%.3f" % a for a in alpha0], sep="\t")

        # updating parameters
        cond_cont = est_contamination and (it % freq_contamination == 0 or it < 3)
        cond_F = est_F and (it % freq_F == 0 or it < 3)
        if cond_F or cond_cont:
            update_post_geno(PG, SNP, Z, IX)
        if cond_cont:
            # need P(G, Z|O') =  P(Z | O') P(G | Z, O')
            delta = update_contamination(cont, error, P, PG, IX, libs)
            if delta < 1e-5:  # when we converged, do not update contamination
                est_contamination, cond_cont = False, False
                print("stopping contamination updates")
        if cond_F:
            # need P(G, Z | O') =  P(Z| O') P(G | Z, O')
            delta = update_F(F, PG, P, IX)
            if delta < 1e-5:  # when we converged, do not update F
                est_F, cond_F = False, False
                print("stopping F updates")
        if cond_F or cond_cont:
            update_snp_prob(SNP, P, IX, cont, error, F, haploid)  # P(O, G | Z)
            e_scaling = update_emissions(E, SNP, P, IX, haploid)  # P(O | Z)
            print("e-scaling:", e_scaling)

    pars = Pars(alpha0, trans_mat, dict(cont), error, F, gamma_names, sex,
                haploid)
    return Z, PG, pars, ll, emissions, (alpha, beta, n)


def run_admixfrog(
    infile,
    ref_file,
    state_ids=("AFR", "VIN", "DEN"),
    cont_id="AFR",
    split_lib=True,
    bin_size=1e4,
    prior=1e-5,
    ancestral=None,
    sex=None,
    pos_mode=False,
    autosomes_only=False,
    downsample=1,
    F0=0,
    e0=1e-2,
    c0=1e-2,
    run_penalty=0.9,
    n_post_replicates=100,
    haploid=False,
    **kwargs
):

    bin_size = bin_size if pos_mode else bin_size * 1e-6

    data = load_data(infile, split_lib, downsample)
    ref = load_ref(ref_file, state_ids, cont_id, prior, ancestral, autosomes_only)
    if pos_mode:
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
    ref = ref.merge(data.iloc[:, :2].drop_duplicates()).drop_duplicates()
    print(ref.shape)
    data = data.merge(ref)
    print(data.shape)
    ref = ref.sort_values(["chrom", "map", "pos"])
    data = data.sort_values(["chrom", "map", "pos"])
    bins, IX = bins_from_bed(
        bed=ref.iloc[:, :5], data=data, bin_size=bin_size, pos_mode=pos_mode, sex=sex
    )
    P = data2probs(data, ref, state_ids, cont_id, (prior, prior))
    assert ref.shape[0] == P.alpha.shape[0]
    del ref

    pars = init_pars(state_ids, sex, F0, e0, c0, haploid)
    print(pars.alpha0)
    print("done loading data")

    Z, G, pars, ll, emissions, (alpha, beta, n) = bw_bb(P, IX, pars,
                                                        haploid=haploid, **kwargs)

    pickle.dump((alpha, beta, n, emissions, pars), open("dump.pickle", "wb"))

    viterbi_path = viterbi(pars, emissions)

    df_pred = pred_sims(
        trans=pars.trans_mat,
        emissions=emissions,
        beta=beta,
        alpha0=pars.alpha0,
        n=n,
        n_homo=len(state_ids),
        n_sims=n_post_replicates,
    )

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
    CC = Counter(data.lib)
    snp = pd.DataFrame([CC[l] for l in df_libs.lib], columns=["n_snps"])
    df_libs = pd.concat((df_libs, snp), axis=1)
    df_libs.sort_values("n_snps", ascending=False)

    df_pars = pd.DataFrame(pars.trans_mat, columns=pars.gamma_names)
    df_pars["alpha0"] = pars.alpha0
    df_pars["state"] = pars.gamma_names
    df_pars["F"] = 0
    df_pars["ll"] = ll
    for i in range(len(pars.F)):
        df_pars.loc[i, "F"] = pars.F[i]

    df_rle = get_rle(df, state_ids, run_penalty)

    return df, snp_df, df_libs, df_pars, df_rle, df_pred
