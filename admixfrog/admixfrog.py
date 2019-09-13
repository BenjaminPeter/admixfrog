import numpy as np
import pandas as pd
import itertools
import yaml
from collections import Counter, defaultdict
from .io import load_read_data, load_ref, filter_ref
from .io import write_bin_table, write_pars_table, write_cont_table
from .io import write_snp_table, write_est_runs, write_sim_runs
from .utils import bins_from_bed, data2probs, init_pars, Pars, ParsHD
from .utils import guess_sex, parse_state_string
from .fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions
from .genotype_emissions import update_post_geno, update_Ftau, update_snp_prob
from .genotype_emissions import update_emissions
from .read_emissions import update_contamination
from .rle import get_rle
from .decode import pred_sims, resampling_pars
from .log import log_, setup_log
from .geno_io import read_geno_ref, read_geno

COORDS = ["chrom", "map", "pos"]

EST_DEFAULT = dict(
    [
        ("est_inbreeding", False),
        ("est_contamination", True),
        ("est_F", False),
        ("est_tau", False),
        ("est_error", False),
        ("freq_contamination", 1),
        ("freq_F", 1),
    ]
)


def baum_welch(P, IX, pars, est_options, max_iter=1000, ll_tol=1e-1, gt_mode=False):
    O = est_options
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
    # if gll_mode:
    SNP = np.zeros((IX.n_snps, n_states, n_gt))  # P(O, G | Z)
    PG = np.zeros((IX.n_snps, n_states, n_gt))  # P(G Z | O)
    # else:
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


    s_scaling = update_snp_prob(
        SNP, P, IX, cont, error, F, tau, O["est_inbreeding"], gt_mode
    )

    e_scaling = update_emissions(E, SNP, P, IX)  # P(O | Z)
    scaling = e_scaling + s_scaling

    for it in range(max_iter):
        alpha, beta, n = fwd_bwd_algorithm(alpha0, emissions, trans, gamma)
        alpha_hap, beta_hap, n_hap = fwd_bwd_algorithm(
            alpha0_hap, hap_emissions, trans_hap, hap_gamma
        )

        ll, old_ll = (
            np.sum([np.sum(np.log(n_i)) for n_i in itertools.chain(n, n_hap)])
            + scaling,
            ll,
        )
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
            trans, alpha, beta, gamma, emissions, n, est_inbreeding=O["est_inbreeding"]
        )
        trans_hap = update_transitions(
            trans_hap,
            alpha_hap,
            beta_hap,
            hap_gamma,
            hap_emissions,
            n_hap,
            est_inbreeding=False,
        )
        alpha0 = np.linalg.matrix_power(trans, 10000)[0]
        alpha0_hap = np.linalg.matrix_power(trans_hap, 10000)[0]

        if gamma_names is not None:
            log_.info("\t".join(gamma_names))
        log_.info("\t".join(["%.3f" % a for a in alpha0]))

        scaling = update_emission_stuff(
            it, E, P, PG, SNP, Z, IX, cont, error, F, tau, scaling, est_options, gt_mode
        )
        #breakpoint()

    update_post_geno(PG, SNP, Z, IX)
    pars = ParsHD(
        alpha0, alpha0_hap, trans, trans_hap, cont, error, F, tau, gamma_names, sex
    )
    return (
        Z,
        PG,
        pars,
        ll,
        emissions,
        hap_emissions,
        (alpha, beta, n),
        (alpha_hap, beta_hap, n_hap),
    )


def update_emission_stuff(
    it, E, P, PG, SNP, Z, IX, cont, error, F, tau, scaling, est_options, gt_mode
):
    O = est_options
    cond_cont = O["est_contamination"] and (it % O["freq_contamination"] == 0 or it < 3)
    cond_Ftau = O["est_F"] or O["est_tau"] and (it % O["freq_F"] == 0 or it < 3)

    if cond_Ftau or cond_cont:  # and gll_mode:
        update_post_geno(PG, SNP, Z, IX)

    if cond_cont and not gt_mode:
        delta = update_contamination(cont, error, P, PG, IX, est_options)
        if delta < 1e-5:  # when we converged, do not update contamination
            O["est_contamination"], cond_cont = False, False
            log_.info("stopping contamination updates")

    if cond_Ftau:
        delta = update_Ftau(F, tau, PG, P, IX, est_options)
        if delta < 1e-5:  # when we converged, do not update F
            O["est_F"], O["est_tau"] = False, False
            cond_Ftau, cond_F = False, False
            log_.info("stopping Ftau updates")
    if cond_Ftau or cond_cont:
        s_scaling = update_snp_prob(
            SNP,
            P,
            IX,
            cont,
            error,
            F,
            tau,
            est_inbreeding=O["est_inbreeding"],
            gt_mode=gt_mode,
        )

        e_scaling = update_emissions(E, SNP, P, IX)  # P(O | Z)
        log_.info("e-scaling: %s", e_scaling)
        log_.info("s-scaling: %s", s_scaling)
        scaling = e_scaling + s_scaling

    return scaling



def load_admixfrog_data(states,
                        state_dict=None,
                        target=None,
                        target_file=None, 
                        gt_mode=False,
                        filter=defaultdict(lambda: None),
                        ref_files=None,
                        geno_file=None,
                        ancestral=None,
                        cont_id=None,
                        split_lib=True,
                        pos_mode=False,
                        downsample=1,
                        guess_ploidy=True,
                        map_col='map',
                        autosomes_only=False):
    """
        we have the following possible input files
        1. target (standalone csv)
        2. ref (standalone csv)
        3. geno (target and/or ref)
    """

    breakpoint()

    "1. only geno file"
    if ref_files is None and target_file is None and geno_file and target:
        df = read_geno_ref(fname=geno_file, pops=state_dict, 
                           target_ind=target, guess_ploidy=guess_ploidy)
        df['lib'] = 'lib0'
        "2. standard ref and target from geno"
    elif ref_files and target_file is None and geno_file and target:
        raise NotImplementedError("")
        "3. standard input"
    elif ref_files and target_file and geno_file is None and target is None:
        if gt_mode:  # gt mode does not do read emissions, assumes genotypes are known
            data = load_read_data(target_file, high_cov_filter=filter.pop('filter_high_cov'))
        else:
            #breakpoint()
            data = load_read_data(target_file, split_lib, downsample, high_cov_filter=filter.pop('filter_high_cov'))

        ref = load_ref(ref_files, state_dict, cont_id, ancestral,
                       autosomes_only, map_col=map_col)
        ref = filter_ref(ref, states, **filter)
        if pos_mode:
            ref.reset_index('map', inplace=True)
            ref.map = ref.index.get_level_values('pos')
            ref.set_index('map', append=True, inplace=True)
        ref = ref.loc[~ref.index.duplicated()]
        df = ref.join(data, how='inner')


        "4. geno ref, standard target"
    elif ref_files is None and geno_file is None and target_file and target is None:
        raise NotImplementedError("")
    else:
        raise ValueError("ambiguous input")

    #get ids of unique snps
    snp_ids = df[~df.index.duplicated()].groupby(df.index.names).ngroup()
    snp_ids = snp_ids.rename('snp_id')
    snp_ids = pd.DataFrame(snp_ids)                              
    snp_ids.set_index('snp_id', append=True, inplace=True)       
    df = snp_ids.join(df)

    df.sort_index(inplace=True)

    #breakpoint()

    return df





def run_admixfrog(
    target_file,
    ref_files,
    geno_file=None,
    target=None,
    states=("AFR", "VIN", "DEN"),
    state_file = None,
    cont_id="AFR",
    split_lib=True,
    bin_size=1e4,
    prior=None,
    ancestral=None,
    sex=None,
    pos_mode=False,
    autosomes_only=False,
    map_col='map',
    downsample=1,
    n_post_replicates=100,
    gt_mode=False,
    keep_loc=True,
    output=defaultdict(lambda: True),
    outname='admixfrog',
    init=defaultdict(lambda: 1e-2),
    guess_ploidy=False,
    est=EST_DEFAULT,
    filter=defaultdict(lambda: None),
    **kwargs
):
    """admixture fragment inference
    this is typically run through the command-line interface. Type admixfrog --help for information
    on arguments
    """



    # numpy config
    np.set_printoptions(suppress=True, precision=4)
    np.seterr(divide="ignore", invalid="ignore")

    if (not est["est_contamination"] and init["c0"] == 0) or gt_mode:
        cont_id = None

    state_dict = parse_state_string(states + [ancestral, cont_id], state_file=state_file) 
    state_dict2 = parse_state_string(states, state_file=state_file) 
    states = list(set(state_dict2.values()))

    # by default, bin size is scaled by 10^6 - could be changed
    bin_size = bin_size if pos_mode else bin_size * 1e-6

    df = load_admixfrog_data(target_file = target_file,
                             ref_files=ref_files,
                             geno_file=geno_file,
                             target=target,
                             gt_mode = gt_mode,
                             states=states,
                             state_dict = state_dict,
                             ancestral=ancestral,
                             cont_id=cont_id,
                             split_lib=split_lib,
                             map_col=map_col,
                             pos_mode=pos_mode,
                             downsample=downsample,
                             guess_ploidy=guess_ploidy,
                             filter=filter,
                             autosomes_only=autosomes_only)
    log_.info("done loading data")


    # sexing stuff
    if sex is None:
        sex = guess_sex(df)


    bins, IX = bins_from_bed(df, bin_size=bin_size, sex=sex)
    log_.info("done creating bins")
    #breakpoint()

    
    P = data2probs(df, IX, states, cont_id, prior=prior, ancestral=ancestral)
    log_.info("done creating prior")

    pars = init_pars(states, sex, est_inbreeding=est["est_inbreeding"], **init)

    Z, G, pars, ll, emissions, hemissions, (_, beta, n), (_, bhap, nhap) = baum_welch(
        P, IX, pars, gt_mode=gt_mode, est_options=est, **kwargs
    )

    # output formating from here
    if output["output_pars"]:
        df_pars = write_pars_table(pars, outname=f'{outname}.pars.yaml')

    if output["output_rsim"]:
        df_pred = pred_sims(
            trans=pars.trans,
            emissions=emissions,
            beta=beta,
            alpha0=pars.alpha0,
            n=n,
            n_homo=len(states),
            n_sims=n_post_replicates,
            est_inbreeding=est["est_inbreeding"],
            keep_loc=keep_loc,
        )
        df_pred["chrom"] = [IX.diplo_chroms[i] for i in df_pred.chrom.values]
        df_pred["state"] = [states[i] for i in df_pred.state.values]

        if len(bhap) > 0:
            df_pred_hap = pred_sims(
                trans=pars.trans_hap,
                emissions=hemissions,
                beta=bhap,
                alpha0=pars.alpha0_hap,
                n=nhap,
                n_homo=len(states),
                n_sims=n_post_replicates,
                est_inbreeding=est["est_inbreeding"],
                keep_loc=keep_loc,
                decode=False,
            )
            df_pred_hap["chrom"] = [
                IX.haplo_chroms[i] for i in df_pred_hap.chrom.values
            ]
            df_pred_hap["state"] = [states[i] for i in df_pred_hap.state.values]
            df_pred = pd.concat((df_pred, df_pred_hap))

        write_sim_runs(df_pred, outname=f'{outname}.res.xz')
        resampling_pars(df_pred).to_csv(f'{outname}.res2.xz', compression='xz', index=True, float_format="%.6f")

    if output["output_bin"] or output["output_rle"]:
        viterbi_path = viterbi(pars.alpha0, pars.trans, emissions)
        viterbi_path_hap = viterbi(pars.alpha0_hap, pars.trans_hap, hemissions)
        V = np.array(pars.gamma_names)[np.hstack(viterbi_path)]
        viterbi_df = pd.Series(V, name="viterbi")

        df_bin = write_bin_table(Z, bins, viterbi_df, pars.gamma_names, IX, outname=f'{outname}.bin.xz')

    if output["output_snp"]:
        df_snp = write_snp_table(data=df, G=G, Z=Z, IX=IX, gt_mode=gt_mode, outname=f'{outname}.snp.xz')

    if output["output_cont"]:
        df_cont = write_cont_table(df, pars.cont, pars.error, outname=f'{outname}.cont.xz')

    if output["output_rle"]:
        df_rle = get_rle(df_bin, states, init["run_penalty"])
        write_est_runs(df_rle, outname=f'{outname}.rle.xz')

    return 
