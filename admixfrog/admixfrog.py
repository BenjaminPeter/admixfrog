import numpy as np
import pickle
import pandas as pd
from collections import Counter
from .io import load_read_data, load_gt_data, load_ref, filter_ref
from .utils import bins_from_bed, data2probs, init_pars, Pars, ParsHD
from .utils import posterior_table,  guess_sex
from .fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions
from .genotype_emissions import update_post_geno, update_F, update_snp_prob
from .genotype_emissions import update_emissions, update_Ftau, update_tau
from .gtmode_emissions import update_geno_emissions_gt
from .read_emissions import update_contamination
from .rle import get_rle
from .decode import pred_sims
from .log import log_, setup_log
import itertools

COORDS = ['chrom', 'map', 'pos']



def baum_welch(
    P,
    IX,
    pars,
    max_iter=1000,
    ll_tol=1e-1,
    est_contamination=True,
    est_F=True,
    est_tau=True,
    freq_contamination=1,
    freq_F=1,
    est_inbreeding=False,
    gt_mode=False
):
    breakpoint()
    alpha0, alpha0_hap, trans, trans_hap, cont, error, F, tau, gamma_names, sex = pars
    gll_mode = not gt_mode
    ll = -np.inf
    n_states = len(alpha0)
    n_hap_states = len(alpha0_hap)
    n_gt = 3 if gt_mode else 3

    # create arrays for posterior, emissions
    Z = np.zeros((sum(IX.bin_sizes), n_states))  # P(Z | O)
    E = np.ones((sum(IX.bin_sizes), n_states))  # P(O | Z)
    # P(O, G | Z), scaled such that max for each row is 1
    #if gll_mode:
    SNP = np.zeros((IX.n_snps, n_states, n_gt))  # P(O, G | Z)
    PG = np.zeros((IX.n_snps, n_states, n_gt))  # P(G Z | O) 
    #else:
    #    SNP = np.zeros((IX.n_snps, n_states))  # P(G | Z)

    gamma, emissions = [], []
    hap_gamma, hap_emissions = [], []
    row0 = 0
    for r, chrom in zip(IX.bin_sizes, IX.chroms):
        if chrom in IX.haplo_chroms:
            hap_gamma.append(Z[row0 : (row0 + r), :n_hap_states])
            hap_emissions.append(E[row0 : (row0 + r), :n_hap_states])
        else:
            gamma.append(Z[row0 : (row0 + r)])
            emissions.append(E[row0 : (row0 + r)])
        row0 += r

    breakpoint()
    #if gll_mode:
    s_scaling = update_snp_prob(
        SNP, P, IX, cont, error, F, tau, est_inbreeding, gt_mode
    )  # P(O, G | Z)
    #else:
    #    s_scaling = update_geno_emissions_gt(
    #        SNP, P, IX, F, tau, est_inbreeding
    #    )  # P(O, G | Z)

    e_scaling = update_emissions(E, SNP, P, IX)  # P(O | Z)
    scaling = e_scaling + s_scaling

    for it in range(max_iter):
        alpha, beta, n = fwd_bwd_algorithm(alpha0, emissions, trans, gamma)
        alpha_hap, beta_hap, n_hap = fwd_bwd_algorithm(alpha0_hap, 
                                                       hap_emissions, 
                                                       trans_hap, 
                                                       hap_gamma)

        ll, old_ll = np.sum([np.sum(np.log(n_i)) for n_i in itertools.chain(n, n_hap)]) + scaling, ll
        assert np.allclose(np.sum(Z, 1), 1)
        if np.isnan(ll):
            pass
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
        log_.info("[%dk|%dk|%dk|%dk]: iter:%d |p95:%.3f\tLL:%.4f\tΔLL:%.4f" % tpl)
        if ll - old_ll < ll_tol:
            break

        if np.any(np.isnan(Z)):
            raise ValueError("nan observed in state posterior")
        if np.any(np.isnan(E)):
            raise ValueError("nan observed in emissions")

        # update stuff
        trans = update_transitions(
            trans,
            alpha,
            beta,
            gamma,
            emissions,
            n,
            est_inbreeding=est_inbreeding,
        )
        trans_hap = update_transitions(
            trans_hap,
            alpha_hap,
            beta_hap,
            hap_gamma,
            hap_emissions,
            n_hap,
            est_inbreeding=False
        )
        alpha0 = np.linalg.matrix_power(trans, 10000)[0]
        alpha0_hap = np.linalg.matrix_power(trans_hap, 10000)[0]

        if gamma_names is not None:
            log_.info("\t".join(gamma_names))
        log_.info("\t".join(["%.3f" % a for a in alpha0]))

        scaling, est_contamination, est_F, est_tau = update_emission_stuff(
            it,
            E, P, PG, SNP, Z, IX, 
            cont, error, F, tau, scaling,
            est_contamination, est_F, est_tau, est_inbreeding,
            freq_contamination, freq_F,
            gt_mode)


    pars = ParsHD(alpha0, alpha0_hap, trans, trans_hap, dict(cont), error, F, tau, gamma_names, sex)
    return Z, PG, pars, ll, emissions, hap_emissions, (alpha, beta, n), (alpha_hap, beta_hap, n_hap)

def update_emission_stuff(it, 
                          E, P, PG, SNP, Z, IX, 
                          cont, error, F, tau, scaling,
                          est_contamination, est_F, est_tau,  est_inbreeding,
                          freq_contamination, freq_F,
                          gt_mode):
    cond_cont = est_contamination and (it % freq_contamination == 0 or it < 3)

    cond_F = est_F and (it % freq_F == 0 or it < 3)
    cond_tau = est_tau and (it % freq_F == 0 or it < 3)
    cond_Ftau = cond_F and cond_tau

    if (cond_F or cond_cont or cond_tau): # and gll_mode:
        update_post_geno(PG, SNP, Z, IX)
    #elif gt_mode:
    #    PG = Z
    if cond_cont and not gt_mode:
        # need P(G, Z|O') =  P(Z | O') P(G | Z, O')
        delta = update_contamination(cont, error, P, PG, IX, IX.libs)
        if delta < 1e-5:  # when we converged, do not update contamination
            est_contamination, cond_cont = False, False
            log_.info("stopping contamination updates")
    if cond_Ftau:
        delta = update_Ftau(F, tau, PG, P, IX)
        if delta < 1e-5:  # when we converged, do not update F
            est_F, est_tau = False, False
            cond_Ftau, cond_F, cond_tau = False, False, False
            log_.info("stopping Ftau updates")
    elif cond_F:
        # need P(G, Z | O') =  P(Z| O') P(G | Z, O')
        delta = update_F(F, tau, PG, P, IX)
        if delta < 1e-5:  # when we converged, do not update F
            est_F, cond_F = False, False
            log_.info("stopping F updates")
    elif cond_tau:
        delta = update_tau(F, tau, PG, P, IX)
        if delta < 1e-5:  # when we converged, do not update F
            est_tau, cond_tau = False, False
            log_.info("stopping F updates")
    if cond_F or cond_cont or cond_tau:
        #    if gll_mode:
        s_scaling = update_snp_prob(
            SNP, P, IX, cont, error, F, tau, est_inbreeding=est_inbreeding,
            gt_mode=gt_mode
        )  # P(O, G | Z)
        #    else:
        #        s_scaling = update_geno_emissions_gt(
        #            SNP, P, IX, F, tau, est_inbreeding
    #        )  # P(O, G | Z)
        e_scaling = update_emissions(
            E, SNP, P, IX
        )  # P(O | Z)
        log_.info("e-scaling: %s", e_scaling)
        log_.info("s-scaling: %s", s_scaling)
        scaling = e_scaling + s_scaling

    return scaling, est_contamination, est_tau, est_F


def run_admixfrog(
    infile,
    ref_file,
    state_ids=("AFR", "VIN", "DEN"),
    cont_id="AFR",
    split_lib=True,
    bin_size=1e4,
    prior=None,
    ancestral=None,
    sex=None,
    pos_mode=False,
    autosomes_only=False,
    downsample=1,
    F0=0,
    e0=1e-2,
    c0=1e-2,
    tau0=1,
    run_penalty=0.9,
    n_post_replicates=100,
    ld_weighting=False, #to be removed
    est_inbreeding=False,
    est_contamination=True,
    filter_delta=None,
    filter_pos=None,
    filter_map=None,
    init_guess=None,
    gt_mode=False,
    **kwargs
):

    #numpy config
    np.set_printoptions(suppress=True, precision=4)
    np.seterr(divide='ignore', invalid='ignore')


    #by default, bin size is scaled by 10^6 - could be changed
    bin_size = bin_size if pos_mode else bin_size * 1e-6


    #loading data and reference
    if gt_mode:   # gt mode does not do read emissions, assumes genotypes are known
        data = load_read_data(infile)
    else:
        data = load_read_data(infile, split_lib, downsample)

    if not est_contamination and c0 ==0:
        cont_id = None

    ref = load_ref(ref_file, state_ids, cont_id, prior, ancestral, autosomes_only)
    ref = filter_ref(ref, state_ids, filter_delta, filter_pos, filter_map)
    ref = ref.drop_duplicates(COORDS)
    if pos_mode:
        ref.map = ref.pos

    # sexing stuff
    if sex is None:
        sex = guess_sex(data)


    log_.debug(ref.shape)
    data = data.merge(ref[COORDS], how='inner').dropna()
    log_.debug(data.shape)

    ref = ref.sort_values(COORDS)
    data = data.sort_values(COORDS)

    snp = data[COORDS].drop_duplicates()
    log_.debug(snp.shape)
    n_snps = snp.shape[0]
    snp["snp_id"] = range(n_snps)
    data = data.merge(snp)

    bins, IX = bins_from_bed(
        bed=ref.iloc[:, :5], snp=snp, data=data, bin_size=bin_size, 
        pos_mode=pos_mode, sex=sex 
    )

    ref = ref.merge(snp[COORDS], 'right')
    log_.debug(ref.shape)

    P = data2probs(
        data,
        ref,
        IX,
        state_ids,
        cont_id,
        prior=prior,
        ancestral=ancestral
    )

    assert ref.shape[0] == P.alpha.shape[0] + P.alpha_hap.shape[0]
    del ref, snp

    pars = init_pars(state_ids, sex, F0, tau0, e0, c0, est_inbreeding, init_guess=init_guess)
    log_.info("done loading data")

    Z, G, pars, ll, emissions, hap_emissions, (alpha, beta, n), (ahap, bhap, nhap) = baum_welch(
        P, IX,
        pars,
        est_inbreeding=est_inbreeding,
        est_contamination=est_contamination,
        gt_mode=gt_mode,
        **kwargs
    )

    #pickle.dump((alpha, beta, n, emissions, pars), open("dump.pickle", "wb"))

    viterbi_path = viterbi(pars.alpha0, pars.trans, emissions)
    viterbi_path_hap = viterbi(pars.alpha0_hap, pars.trans_hap,  hap_emissions)


    #output sims
    df_pred = pred_sims(
        trans=pars.trans,
        emissions=emissions,
        beta=beta,
        alpha0=pars.alpha0,
        n=n,
        n_homo=len(state_ids),
        n_sims=n_post_replicates,
        est_inbreeding=est_inbreeding,
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
    if gt_mode:
        snp_df = pd.concat((D, pd.DataFrame(IX.SNP2BIN, columns=["bin"])), axis=1)
    else:
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

    CC = data.groupby(['lib']).agg(({"tref": sum, "talt": sum})).reset_index()
    CC['n_snps'] = CC.tref + CC.talt
    del CC['tref']
    del CC['talt']

    df_libs = df_libs.merge(CC)
    df_libs.sort_values("n_snps", ascending=False)

    df_pars = pd.DataFrame(pars.trans, columns=pars.gamma_names)
    df_pars["alpha0"] = pars.alpha0
    df_pars["state"] = pars.gamma_names
    df_pars["F"] = 0
    df_pars["tau"] = 0
    df_pars["ll"] = ll
    for i in range(len(pars.F)):
        df_pars.loc[i, "F"] = pars.F[i]
    for i in range(len(pars.tau)):
        df_pars.loc[i, "tau"] = pars.tau[i]

    df_rle = get_rle(df, state_ids, run_penalty)

    return df, snp_df, df_libs, df_pars, df_rle, df_pred
