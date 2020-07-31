import numpy as np
import pandas as pd
import itertools
import yaml
from collections import Counter, defaultdict
from .read_emissions2 import p_snps_given_gt
from .io import load_read_data, load_ref, filter_ref
from .io import write_bin_table, write_pars_table, write_cont_table
from .io import write_snp_table, write_est_runs, write_sim_runs, write_sfs
from .utils import bins_from_bed, sfs_from_bed, data2probs, init_pars_sfs, ParsSFS
from .utils import guess_sex, parse_state_string
from .fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions
from .genotype_emissions import update_post_geno,  update_snp_prob
from .genotype_emissions import update_emissions
from .read_emissions import update_contamination
from .rle import get_rle
from .decode import pred_sims, resampling_pars
from .log import log_, setup_log
from .geno_io import read_geno_ref, read_geno
from .admixfrog import load_admixfrog_data

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

def update_sfs_genotypes(G, P, IX, F, tau):
    G[:, 0] = F * (1-tau) + (1-F) * (1-tau) * (1-tau)
    G[:, 1] = (1-F) * 2 * (1-tau) * tau
    G[:, 2] = F * tau + (1-F) * tau * tau
    return 0

def update_snp(E, G, P, IX, cont, error, gt_mode):
    cflat = np.array([cont[lib] for lib in P.lib])
    eflat = np.array([error[lib] for lib in P.lib])

    E[:] = p_snps_given_gt(P, cflat, eflat, IX, gt_mode)
    E[:] = G[IX.SNP2BIN] * E
    return 0

def update_ftau(PG, IX, F, tau):
    for i in range(IX.n_sfs):
        g0, g1, g2 = np.mean(PG[IX.SNP2BIN==i],0)
        F[i] = 2 * g0 / (2 * g0 + g1 + 1e-8) - g1 / (2 * g2 + g1 + 1e-8)
        if np.isnan(F[i]):
            breakpoint()
        F[:] = np.minimum(1, np.maximum(0, F))
        tau[i] = g1/2 + g2

    


def em_sfs(P, IX, pars, est_options, max_iter=1000, 
               ll_tol=1e-1, gt_mode=False,
               do_hmm=True, **kwargs):
    O = est_options
    cont, error, F, tau, sex = pars
    ll = -np.inf
    n_gt = 3

    # create arrays for posterior, emissions
    G = np.zeros((IX.n_sfs, n_gt)) # Pr(G | F, tau)
    E = np.zeros((IX.n_snps, n_gt)) # P(O, G | tau, F)

    #update G: given F, tau, calculate Pr(G | F, tau) for each SFS class
    update_sfs_genotypes(G, P, IX, F, tau) 

    #update E: Pr(O, G | F, tau)
    update_snp(E, G, P, IX, cont, error, gt_mode)

    for it in range(max_iter):
        ll, old_ll = np.sum(np.log(np.sum(E,1))), ll
        if np.isnan(ll):
            breakpoint()
        assert not np.isnan(ll)

        tpl = (
            IX.n_reads / 1000,
            IX.n_obs / 1000,
            IX.n_snps / 1000,
            IX.n_bins / 1000,
            it,
            np.mean(tau),
            ll,
            ll - old_ll,
        )
        log_.info("[%dk|%dk|%dk|%dk]: iter:%d |tavg:%.3f\tLL:%.4f\tΔLL:%.4f" % tpl)
        if ll - old_ll < ll_tol:
            break

        if np.any(np.isnan(E)):
            raise ValueError("nan observed in emissions")

        E[:] = E / np.sum(E, 1)[:, np.newaxis]
        update_ftau(E, IX, F, tau)
        PG = E.reshape(E.shape[0], 1, E.shape[1])
        delta = update_contamination(cont, error, P, PG, IX, est_options)

        update_sfs_genotypes(G, P, IX, F, tau) 
        update_snp(E, G, P, IX, cont, error, gt_mode)


    pars = ParsSFS(
        cont, error, F, tau, sex
    )

    return G, PG, pars, ll






def run_admixslug(
    target_file,
    ref_files,
    geno_file=None,
    target=None,
    states=("AFR", "VIN", "DEN"),
    state_file = None,
    cont_id="AFR",
    split_lib=True,
    prior=None,
    map_col='map',
    ancestral=None,
    ancestral_prior = 0,
    sex=None,
    pos_mode=False,
    autosomes_only=False,
    downsample=1,
    gt_mode=False,
    output=defaultdict(lambda: True),
    outname='admixfrog',
    init=defaultdict(lambda: 1e-2),
    guess_ploidy=False,
    est=EST_DEFAULT,
    fake_contamination=0.,
    filter=defaultdict(lambda: None),
    **kwargs
):
    """contamination estimation using sfs 
    """

    from pprint import pprint
    pprint(kwargs)


    # numpy config
    np.set_printoptions(suppress=True, precision=4)
    np.seterr(divide="ignore", invalid="ignore")

    if (not est["est_contamination"] and init["c0"] == 0) or gt_mode:
        cont_id = None

    state_dict = parse_state_string(states + [ancestral, cont_id], state_file=state_file) 
    state_dict2 = parse_state_string(states, state_file=state_file) 
    states = list(dict(((x, None) for x in state_dict2.values())))


    df, sex, tot_n_snps = load_admixfrog_data(target_file = target_file,
                             ref_files=ref_files,
                             geno_file=geno_file,
                             target=target,
                             gt_mode = gt_mode,
                             states=states,
                             sex=sex,
                             state_dict = state_dict,
                             ancestral=ancestral,
                             cont_id=cont_id,
                             split_lib=split_lib,
                             map_col=map_col,
                             pos_mode=pos_mode,
                             downsample=downsample,
                             guess_ploidy=guess_ploidy,
                             fake_contamination=fake_contamination,
                             filter=filter,
                             autosomes_only=autosomes_only)
    log_.info("done loading data")

    sfs, IX = sfs_from_bed(df, states=states, ancestral=ancestral, sex=sex)
    log_.info("done creating sfs")

    P = data2probs(df, IX, states, cont_id, prior=prior, ancestral=ancestral,
                   ancestral_prior=ancestral_prior,
                   doalphabeta=False)
    log_.info("done creating prior")

    pars = init_pars_sfs(IX.n_sfs, sex=sex, **init)

    Z, G, pars, ll = em_sfs( P, IX, pars, gt_mode=gt_mode, 
        est_options=est, **kwargs
    )

    # output formating from here
    if output["output_pars"]:
        df_pars = write_pars_table(pars, outname=f'{outname}.pars.yaml')

    if output["output_snp"]:
        df_snp = write_snp_table(data=df, G=G, Z=Z, IX=IX, gt_mode=gt_mode, outname=f'{outname}.snp.xz')

    if output["output_cont"]:
        df_cont = write_cont_table(df, pars.cont, pars.error, tot_n_snps, outname=f'{outname}.cont.xz')

    if output["output_sfs"]:
        write_sfs(sfs, Z, pars.tau, pars.F, pars.cont, P, IX, outname=f'{outname}.sfs.xz')

    return 