import numpy as np
import logging
import pandas as pd
import yaml
from pprint import pprint
from collections import Counter, defaultdict
from copy import deepcopy
from ..gll.read_emissions2 import p_snps_given_gt
from ..utils.input import load_read_data, load_ref, filter_ref
from ..utils.output import write_pars_table
from ..utils.output_slug import write_snp_table_slug, write_cont_table_slug
from ..utils.output_slug import write_sfs2, write_vcf
from ..utils.output_slug import write_f3_table, write_f4_table, write_f2_table
from ..utils.utils import data2probs
from ..utils.utils import guess_sex
from ..utils.states import States
from ..gll.genotype_emissions import update_post_geno, update_snp_prob
from ..gll.genotype_emissions import update_emissions
from ..gll.read_emissions import update_contamination
from ..utils.geno_io import read_geno_ref, read_geno
from .classes import SlugController, SlugReads, SlugPars
from .em import em, squarem
from .emissions import full_posterior_genotypes
from .fstats import calc_fstats, summarize_f3, summarize_f4, summarize_pi, summarize_f2

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
    state_file=None,
    cont_id="AFR",
    split_lib=True,
    prior=None,
    ancestral=None,
    sex=None,
    autosomes_only=False,
    downsample=1,
    output=defaultdict(lambda: True),
    outname="admixfrog",
    init=defaultdict(lambda: 1e-2),
    est=EST_DEFAULT,
    fake_contamination=0.0,
    deam_bin_size=50000,
    len_bin_size=1000,
    filter=defaultdict(lambda: None),
    ll_tol=0.001,
    ptol=1e-4,
    max_iter=100,
    jk_resamples=0,
    target="admixslug",
    sex_chroms=None,
    **kwargs,
):
    """contamination estimation using sfs"""

    # numpy config
    np.set_printoptions(suppress=True, precision=4)
    np.seterr(divide="ignore", invalid="ignore")

    controller = SlugController(
        update_eb=est["est_error"],
        update_ftau=est["est_tau"],
        update_F=est["est_F"],
        update_cont=est["est_contamination"],
        update_bias=est["est_bias"],
        n_iter=max_iter,
        ll_tol=ll_tol,
        param_tol=ptol,
        n_resamples=jk_resamples,
    )

    if not est["est_contamination"] and init["c0"] == 0:
        cont_id = None

    states = States.from_commandline(
        raw_states=[*states, ancestral], state_file=state_file, ancestral=ancestral, cont_id=cont_id
    )

    df, ix, sex, n_sites = load_admixslug_data_native(
        target_file=target_file,
        ref_files=ref_files,
        states=states,
        ancestral=ancestral,
        sex=sex,
        cont_id=cont_id,
        split_lib=split_lib,
        downsample=downsample,
        fake_contamination=fake_contamination,
        filter=filter,
        deam_bin_size=deam_bin_size,
        len_bin_size=len_bin_size,
        autosomes_only=autosomes_only,
    )
    logging.info("done loading data")

    data, sfs = SlugReads.load(
        df,
        states=states,
        ancestral=ancestral,
        sex=sex,
        cont_id=cont_id,
        flip=True,
        sex_chroms=sex_chroms,
    )

    pars = SlugPars.from_n(data.n_sfs, data.n_rgs, **init)
    pars0 = deepcopy(pars)

    pars = squarem(pars, data, controller)
    # pars = em(pars, data, controller)
    gt_ll, posterior_gt = full_posterior_genotypes(data, pars)

    if controller.n_resamples > 0 or output["output_fstats"]:
        jk_pars_list = []
        jk_sfs = list()

        for i in range(controller.n_resamples):
            jk_data = data.jackknife_sample(i, controller.n_resamples)
            jk_pars = squarem(pars0, jk_data, controller)
            # jk_pars = em(pars0, jk_data, controller)
            jk_pars_list.append(jk_pars)

            if output["output_jk_sfs"] or output["output_fstats"]:
                jk_sfs_row = write_sfs2(sfs, jk_pars, jk_data)
                jk_sfs_row["rep"] = i
                jk_sfs.append(jk_sfs_row)

            print(f"done with jackknife sample {i+1} / {controller.n_resamples}")

        jk_table = np.vstack(tuple(p.pars for p in jk_pars_list))
        if output["output_jk_sfs"] or output["output_fstats"]:
            jk_sfs = pd.concat(jk_sfs)

        n = np.sum(~np.isnan(jk_table), 0)
        se = np.sqrt(((n - 1) / (n) * np.nansum((jk_table - pars.pars) ** 2, 0)))
        se_pars = deepcopy(pars)
        se_pars._pars = se
    else:
        se_pars = deepcopy(pars)
        se_pars._pars[:] = np.nan

    # output formating from here
    if output["output_fstats"]:
        f2s, f3s, f4s, pis = calc_fstats(jk_sfs, states, name=target)
        df_f2 = write_f2_table(f2s, outname=f"{outname}.f2.jk.xz")
        df_f3 = write_f3_table(f3s, outname=f"{outname}.f3.jk.xz")
        df_f4 = write_f4_table(f4s, outname=f"{outname}.f4.jk.xz")

        pis.to_csv(f"{outname}.pi.jk.xz", float_format="%.6f", index=False)

        pi_summary = summarize_pi(pis)
        pi_summary.to_csv(f"{outname}.pi.xz", float_format="%.6f", index=False)

        f2_summary = summarize_f2(f2s)
        f2_summary.to_csv(f"{outname}.f2.xz", float_format="%.6f", index=False)
        f3_summary = summarize_f3(f3s)
        f3_summary.to_csv(f"{outname}.f3.xz", float_format="%.6f", index=False)
        f4_summary = summarize_f4(f4s)
        f4_summary.to_csv(f"{outname}.f4.xz", float_format="%.6f", index=False)

    if output["output_pars"]:
        df_pars = write_pars_table(pars, outname=f"{outname}.pars.yaml")

    if output["output_cont"]:
        df_cont = write_cont_table_slug(
            ix,
            data.rgs,
            pars.cont,
            n_sites,
            se=se_pars.cont,
            outname=f"{outname}.cont.xz",
        )
        ix.to_csv(f"{outname}.ix.xz", index=False)

    if output["output_sfs"]:
        write_sfs2(
            sfs,
            pars,
            data,
            se_tau=se_pars.tau,
            se_F=se_pars.F,
            outname=f"{outname}.sfs.xz",
        )

    if output["output_jk_sfs"]:
        jk_sfs.to_csv(f"{outname}.jksfs.xz", float_format="%5f", index=False)

    if output["output_snp"]:
        df_snp = write_snp_table_slug(
            df=df, data=data, posterior_gt=posterior_gt, outname=f"{outname}.snp.xz"
        )

    if output["output_vcf"]:
        df_snp = write_vcf(
            df=df,
            data=data,
            posterior_gt=posterior_gt,
            genotype_ll=gt_ll,
            sample_name=target,
            outname=f"{outname}.vcf",
        )

    return


def load_admixslug_data_native(
    states,
    target_file=None,
    filter=defaultdict(lambda: None),
    ref_files=None,
    sex=None,
    ancestral=None,
    cont_id=None,
    split_lib=True,
    pos_mode=False,
    downsample=1,
    map_col="map",
    fake_contamination=0,
    deam_bin_size=20_000,
    len_bin_size=5000,
    autosomes_only=False,
    sex_chroms=[],
):

    data, ix = load_read_data(
        target_file,
        split_lib,
        downsample,
        deam_bin_size=deam_bin_size,
        len_bin_size=len_bin_size,
        high_cov_filter=filter.pop("filter_high_cov"),
        make_bins=True,
    )

    ref = load_ref(
        ref_files,
        states,
        cont_id,
        ancestral,
        autosomes_only,
        map_col=map_col,
        large_ref=False,
    )
    ref = filter_ref(ref, states, ancestral=ancestral, cont=cont_id, **filter)

    if pos_mode:
        ref.reset_index("map", inplace=True)
        ref.map = ref.index.get_level_values("pos")
        ref.set_index("map", append=True, inplace=True)

    n_sites = ref.shape[0]

    # workaround a pandas join bug
    df = ref.join(data, how="inner")
    df = make_snp_ids(df)

    # sexing stuff
    if sex is None:
        sex = guess_sex(ref, data)

    # this might be a little bit problematic because it should ideally be before we
    # filter for existing categories
    if fake_contamination and cont_id:
        add_fake_contamination(df, cont_id, prop_cont=fake_contamination)

    return df, ix, sex, n_sites


def make_snp_ids(df):
    """integer id for each SNP with available data"""
    df.reset_index("rg", inplace=True)
    snp_ids = (
        df[~df.index.duplicated()].groupby(df.index.names, observed=False).ngroup()
    )
    snp_ids = snp_ids.rename("snp_id")
    snp_ids = pd.DataFrame(snp_ids)
    snp_ids.set_index("snp_id", append=True, inplace=True)
    df = snp_ids.join(df)
    df.set_index("rg", inplace=True, append=True)

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

    logging.debug(f"endogenous cov: {mean_endo_cov}")
    logging.debug(f"fake contamination cov: {target_cont_cov}")

    try:
        c_ref = np.random.poisson((1 - f_cont) * target_cont_cov)
        c_alt = np.random.poisson(f_cont * target_cont_cov)
    except ValueError:
        raise ValueError()
    logging.debug(f"Added cont. reads with ref allele: {np.sum(c_ref)}")
    logging.debug(f"Added cont. reads with alt allele: {np.sum(c_alt)}")
    df.tref += c_ref
    df.talt += c_alt
