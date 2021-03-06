import numpy as np
import pandas as pd
import yaml
from pprint import pprint
from collections import Counter, defaultdict
from copy import deepcopy
from .read_emissions2 import p_snps_given_gt
from .io import load_read_data, load_ref, filter_ref
from .io import write_bin_table, write_pars_table, write_cont_table
from .io import write_snp_table, write_est_runs, write_sim_runs, write_sfs
from .io import write_snp_table_slug, write_cont_table_slug
from .io import write_snp_table_slug2, write_sfs2, write_vcf
from .utils import data2probs, init_pars_sfs
from .utils import guess_sex, parse_state_string
from .fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions
from .genotype_emissions import update_post_geno,  update_snp_prob
from .genotype_emissions import update_emissions
from .read_emissions import update_contamination
from .rle import get_rle
from .decode import pred_sims, resampling_pars
from .log import log_, setup_log
from .geno_io import read_geno_ref, read_geno
from .slug_classes import SlugController
from .utils import make_slug_data, make_slug_reads_data
from .slug_em_reads import em, squarem
from .slug_emissions_reads import full_posterior_genotypes

EST_DEFAULT = dict(
    [
        ("est_contamination", True),
        ("est_F", True),
        ("est_tau", True),
        ("est_error", True),
        ("est_bias", True),
    ]
)


def run_admixslug(
    target_file,
    ref_files,
    geno_file=None,
    states=("AFR", "VIN", "DEN"),
    state_file = None,
    cont_id="AFR",
    split_lib=True,
    prior=None,
    ancestral=None,
    sex=None,
    autosomes_only=False,
    downsample=1,
    output=defaultdict(lambda: True),
    outname='admixfrog',
    init=defaultdict(lambda: 1e-2),
    est=EST_DEFAULT,
    fake_contamination=0.,
    bin_reads= True,
    deam_bin_size = 50000,
    len_bin_size = 1000,
    filter=defaultdict(lambda: None),
    ll_tol = 0.001,
    ptol = 1e-4,
    max_iter = 100,
    jk_resamples = 0,
    vcf_sample_name = "admixslug",
    **kwargs
):
    """contamination estimation using sfs 
    """

    pprint(kwargs)
    # numpy config
    np.set_printoptions(suppress=True, precision=4)
    np.seterr(divide="ignore", invalid="ignore")


    controller = SlugController(update_eb = est['est_error'] ,
                                update_ftau = est['est_tau'],
                                update_F = est['est_F'],
                                update_cont = est['est_contamination'],
                                update_bias = est['est_bias'],
                                n_iter=max_iter,
                                ll_tol = ll_tol,
                                param_tol = ptol,
                                n_resamples = jk_resamples
                                )




    if (not est["est_contamination"] and init["c0"] == 0):
        cont_id = None

    state_dict = parse_state_string(states + [ancestral, cont_id], state_file=state_file) 
    state_dict2 = parse_state_string(states, state_file=state_file) 
    states = list(dict(((x, None) for x in state_dict2.values())))


    df, sex, n_sites, ix = load_admixslug_data_native(target_file = target_file,
                             ref_files=ref_files,
                             states=states,
                             state_dict = state_dict,
                             ancestral=ancestral,
                             sex = sex,
                             cont_id=cont_id,
                             split_lib=split_lib,
                             downsample=downsample,
                             fake_contamination=fake_contamination,
                             filter=filter,
                             make_bins=bin_reads,
                             deam_bin_size = deam_bin_size,
                             len_bin_size = len_bin_size,
                             autosomes_only=autosomes_only)
    log_.info("done loading data")



    data, sfs = make_slug_reads_data(df, states=states, ancestral=ancestral, sex=sex,
                          cont_id=cont_id, flip=True)
    pars = init_pars_sfs(data.n_sfs, data.n_rgs, **init)
    pars0 = deepcopy(pars)

    #posterior_gt = em(pars, data, controller)
    pars = squarem(pars, data, controller)
    gt_ll, posterior_gt = full_posterior_genotypes(data, pars)

    if controller.n_resamples > 0:
        jk_pars_list = []

        for i in range(controller.n_resamples):
            jk_data = data.jackknife_sample(i, controller.n_resamples)
            jk_pars = squarem(pars0, jk_data, controller)
            jk_pars_list.append(jk_pars)
            print(f"done with jackknife sample {i+1} / {controller.n_resamples}")

        jk_table = np.vstack(tuple(p.pars for p in jk_pars_list))

        n = np.sum(~np.isnan(jk_table), 0)
        se = np.sqrt(((n-1)/(n)*np.nansum((jk_table - pars.pars) ** 2, 0)))
        se_pars = deepcopy(pars)
        se_pars._pars = se
    else:
        se_pars = deepcopy(pars)
        se_pars._pars[:] = np.nan




    # output formating from here
    if output["output_pars"]:
        df_pars = write_pars_table(pars, outname=f'{outname}.pars.yaml')


    if output["output_cont"]:
        df_cont = write_cont_table_slug(ix, data.rgs, pars.cont, 
                                   n_sites, 
                                        se= se_pars.cont,
                                        outname=f'{outname}.cont.xz')

    if output["output_sfs"]:
        write_sfs2(sfs, pars, data,
                   se_tau=se_pars.tau, se_F = se_pars.F, 
                   outname=f'{outname}.sfs.xz')

    if output["output_snp"]:
        df_snp = write_snp_table_slug2(df=df, data=data, posterior_gt = posterior_gt, outname=f'{outname}.snp.xz')

    if output["output_vcf"]:
        df_snp = write_vcf(df=df, 
                           data=data, 
                           posterior_gt = posterior_gt, 
                           genotype_ll = gt_ll,
                           sample_name = vcf_sample_name,
                           outname=f'{outname}.vcf')

    return 

def load_admixslug_data_native(states,
                        state_dict=None,
                        target_file=None, 
                        filter=defaultdict(lambda: None),
                        ref_files=None,
                        sex=None,
                        ancestral=None,
                        cont_id=None,
                        split_lib=True,
                        pos_mode=False,
                        downsample=1,
                        map_col='map',
                        fake_contamination=0,
                        make_bins=True,
                        deam_bin_size = 20_000,
                        len_bin_size = 5000,
                        autosomes_only=False):
     
    data, ix = load_read_data(target_file, split_lib, downsample, 
                              make_bins=make_bins,
                              deam_bin_size = deam_bin_size,
                              len_bin_size = len_bin_size,
                              high_cov_filter=filter.pop('filter_high_cov'))

    data = data[['tref', 'talt', 'rg']]

    ref = load_ref(ref_files, state_dict, cont_id, ancestral,
                   autosomes_only, map_col=map_col, large_ref=False)
    ref = filter_ref(ref, states, ancestral=ancestral, cont=cont_id, **filter)


    if pos_mode:
        ref.reset_index('map', inplace=True)
        ref.map = ref.index.get_level_values('pos')
        ref.set_index('map', append=True, inplace=True)



    #workaround a pandas join bug
    cats = ref.index.get_level_values(0).categories
    ref.index = ref.index.set_levels(ref.index.levels[0].astype(str), level=0)
    data.index = data.index.set_levels(data.index.levels[0].astype(str), level=0)


    df = ref.join(data, how='inner')
    df = make_snp_ids(df)

    # sexing stuff
    if sex is None:
        sex = guess_sex(ref, data)

    if fake_contamination and cont_id:
        add_fake_contamination(df, cont_id, prop_cont=fake_contamination)

    n_sites = df.shape[0]

    return df, sex, n_sites, ix

def make_snp_ids(df):
    """integer id for each SNP with available data"""
    snp_ids = df[~df.index.duplicated()].groupby(df.index.names).ngroup()
    snp_ids = snp_ids.rename('snp_id')
    snp_ids = pd.DataFrame(snp_ids)                              
    snp_ids.set_index('snp_id', append=True, inplace=True)       
    df = snp_ids.join(df)

    df.sort_index(inplace=True)
    return df

def add_fake_contamination(df, cont_id, prop_cont):
        cont_ref, cont_alt = f"{cont_id}_ref", f"{cont_id}_alt"
        mean_endo_cov = np.mean(df.tref + df.talt)
        """
            C = x / (e+x);
            Ce + Cx = x
            Ce = x - Cx
            Ce = x(1-C)
            Ce / ( 1- C) = x
        """

        target_cont_cov = prop_cont * mean_endo_cov / (1 - prop_cont)
        f_cont = df[cont_alt] / (df[cont_ref] + df[cont_alt])

        log_.debug(f"endogenous cov: {mean_endo_cov}")
        log_.debug(f"fake contamination cov: {target_cont_cov}")


        try:
            c_ref = np.random.poisson( (1 - f_cont) * target_cont_cov)
            c_alt = np.random.poisson( f_cont * target_cont_cov)
        except ValueError:
            breakpoint()
            raise ValueError()
        log_.debug(f"Added cont. reads with ref allele: {np.sum(c_ref)}")
        log_.debug(f"Added cont. reads with alt allele: {np.sum(c_alt)}")
        df.tref += c_ref
        df.talt += c_alt
        


