import numpy as np
import pandas as pd
import itertools
import yaml
import logging
from pprint import pprint
from collections import Counter, defaultdict
from ..utils.input import load_read_data, load_ref, filter_ref
from ..utils.input import load_admixfrog_data_geno, load_admixfrog_data
from ..utils.output import write_pars_table
from ..utils.output_frog import write_bin_table, write_cont_table_frog
from ..utils.output_frog import write_snp_table, write_est_runs, write_sim_runs
from ..utils.output_slug import write_cont_table_slug
from ..utils.utils import bins_from_bed, data2probs, init_pars
from ..utils.states import States
from .fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions, fwd_bwd_algorithm2
from ..gll.genotype_emissions import update_post_geno, update_Ftau, update_snp_prob
from ..gll.genotype_emissions import update_emissions
from ..gll.read_emissions import update_contamination
from .rle import get_rle
from .decode import pred_sims, resampling_pars
from ..utils.classes import FrogOptions, FrogX
from ..utils.squarem import squarem
from copy import deepcopy

COORDS = ["chrom", "map", "pos"]

EST_DEFAULT = dict(
    [
        ("est_inbreeding", False),
        ("est_contamination", True),
        ("est_F", False),
        ("est_tau", False),
        ("est_error", False),
    ]
)


def bw_one_iter(pars, data, O, X):
    pars = deepcopy(pars) if O.copy_pars else pars

    H = X.H
    fwd_bwd_algorithm2(pars, X)
    fwd_bwd_algorithm2(pars.H, X.H)

    pars.ll, pars.prev_ll = (
        np.sum([np.sum(np.log(n_i)) for n_i in itertools.chain(X.n, H.n)]) + X.scaling,
        pars.ll,
    )

    assert np.allclose(np.sum(X.Z, 1), 1)
    assert not np.isnan(pars.ll)

    tpl = (
        data.n_reads / 1000,
        data.n_obs / 1000,
        data.n_snps / 1000,
        data.n_bins / 1000,
        np.mean(np.max(X.Z, 1) >= 0.95),
        pars.ll,
        pars.ll - pars.prev_ll,
    )
    logging.info("[%dk|%dk|%dk|%dk]: |p95:%.3f\tLL:%.4f\tΔLL:%.4f" % tpl)
    # if pars.ll - pars.prev_ll < O.ll_tol:
    #    return pars

    if np.any(np.isnan(X.Z)):
        raise ValueError("nan observed in state posterior")
    if np.any(np.isnan(X.E)):
        raise ValueError("nan observed in emissions")

    # update stuff
    if O.est_trans:
        pars.trans = update_transitions(pars, X, O)
        pars.trans_hap = update_transitions(pars.H, X.H, O)
        pars.alpha0 = np.linalg.matrix_power(pars.trans, 10000)[0]
        pars.alpha0_hap = np.linalg.matrix_power(pars.trans_hap, 10000)[0]

    logging.info("\t".join(data.states.state_names))
    logging.info("\t".join(["%.3f" % a for a in pars.alpha0]))

    update_full_emissions(data, X, pars, O)
    return pars


def update_full_emissions(P, X, pars, O):
    if O.est_contamination or O.est_F or O.est_tau or O.est_error:
        update_post_geno(X, P)

    if O.est_contamination and not O.gt_mode:
        update_contamination(pars.cont, pars.error, P, X.PG, O)

    if O.est_F or O.est_tau:
        update_Ftau(pars.F, pars.tau, X.PG, P, O)

    if O.est_contamination or O.est_F or O.est_tau or O.est_error:
        s_scaling = update_snp_prob(
            X.SNP,
            P,
            pars=pars,
            est_inbreeding=O.est_inbreeding,
            gt_mode=O.gt_mode,
            scale_probs=O.scale_probs,
        )

        e_scaling = update_emissions(
            X.E, X.SNP, P, scale_probs=O.scale_probs
        )  # P(O | Z)
        logging.info("e-scaling: %s", e_scaling)
        logging.info("s-scaling: %s", s_scaling)
        X.scaling = e_scaling + s_scaling

    return X.scaling


def run_admixfrog(
    target_file,
    ref_files,
    geno_file=None,
    target=None,
    states=("AFR", "VIN", "DEN"),
    state_file=None,
    homo_states=None,
    het_states=None,
    cont_id="AFR",
    split_lib=True,
    bin_size=1e4,
    snp_mode=False,
    prior=None,
    ancestral=None,
    ancestral_prior=0,
    sex=None,
    pos_mode=False,
    autosomes_only=False,
    map_col="map",
    downsample=1,
    n_post_replicates=100,
    gt_mode=False,
    keep_loc=True,
    output=defaultdict(lambda: True),
    outname="admixfrog",
    init=defaultdict(lambda: 1e-2),
    guess_ploidy=False,
    est=EST_DEFAULT,
    fake_contamination=0.0,
    filter=defaultdict(lambda: None),
    bin_reads=False,
    deam_bin_size=50000,
    len_bin_size=1000,
    **kwargs,
):
    """admixture fragment inference
    this is typically run through the command-line interface. Type admixfrog --help for information
    on arguments
    """
    pprint(kwargs)

    # numpy config for output
    np.set_printoptions(suppress=True, precision=4)
    np.seterr(divide="ignore", invalid="ignore")

    if (not est["est_contamination"] and init["c0"] == 0) or gt_mode:
        cont_id = None

    states = States.from_commandline(
        raw_states=states,
        state_file=state_file,
        ancestral=ancestral,
        cont_id=cont_id,
        homo_states=homo_states,
        het_states=het_states,
        est_inbreeding=est["est_inbreeding"],
    )

    # by default, bin size is scaled by 10^6 - could be changed
    bin_size = bin_size if pos_mode else bin_size * 1e-6

    df, ix, sex, tot_n_snps = load_admixfrog_data(
        target_file=target_file,
        ref_files=ref_files,
        geno_file=geno_file,
        target=target,
        gt_mode=gt_mode,
        states=states,
        sex=sex,
        ancestral=ancestral,
        cont_id=cont_id,
        split_lib=split_lib,
        map_col=map_col,
        pos_mode=pos_mode,
        downsample=downsample,
        guess_ploidy=guess_ploidy,
        fake_contamination=fake_contamination,
        filter=filter,
        autosomes_only=autosomes_only,
        bin_reads=bin_reads,
        deam_bin_size=deam_bin_size,
        len_bin_size=len_bin_size,
    )
    logging.info("done loading data")

    bins, IX = bins_from_bed(df, bin_size=bin_size, sex=sex, snp_mode=snp_mode)
    logging.info("done creating bins")

    P = data2probs(
        df,
        IX,
        states,
        cont_id,
        prior=prior,
        ancestral=ancestral,
        ancestral_prior=ancestral_prior,
    )
    logging.info("done creating prior")

    pars = init_pars(
        states,
        n_rgs=P.n_rgs,
        homo_ids=homo_states,
        het_ids=het_states,
        bin_size=bin_size,
        **init,
    )

    X = FrogX(P)
    O = FrogOptions(**est, gt_mode=gt_mode)

    # initialize emissions latent variable
    # scaling = update_full_emissions(P, X, pars, O)
    s_scaling = update_snp_prob(
        SNP=X.SNP,
        P=P,
        pars=pars,
        est_inbreeding=O.est_inbreeding,
        gt_mode=O.gt_mode,
        scale_probs=O.scale_probs,
    )  # return value here is SNP

    e_scaling = update_emissions(X.E, X.SNP, P, scale_probs=O.scale_probs)  # P(O | Z)
    logging.info("e-scaling: %s", e_scaling)
    logging.info("s-scaling: %s", s_scaling)
    X.scaling = e_scaling + s_scaling

    # run accelerted EM
    pars = squarem(pars, data=P, latents=X, controller=O, updater=bw_one_iter)

    # final update of latent variables
    update_post_geno(X, P)

    # output formating from here
    if output["output_pars"]:
        df_pars = write_pars_table(pars, outname=f"{outname}.pars.yaml")

    if output["output_cont"]:
        if ix is None:
            df_cont = write_cont_table_frog(
                df,
                P.rgs,
                pars.cont,
                pars.error,
                tot_n_snps,
                outname=f"{outname}.cont.xz",
            )
        else:
            rgs = dict((l, i) for i, l in enumerate(P.rgs))

            df_cont = write_cont_table_slug(
                ix,
                rgs,
                pars.cont,
                tot_n_snps,
                se=None,
                outname=f"{outname}.cont.xz",
            )
            ix.to_csv(f"{outname}.ix.xz", index=False)

    if output["output_bin"] or output["output_rle"]:
        viterbi_path = viterbi(pars.alpha0, pars.trans, X.emissions)
        viterbi_path_hap = viterbi(pars.alpha0_hap, pars.trans_hap, X.H.emissions)
        V = np.array([*states.state_names])[np.hstack(viterbi_path)]
        viterbi_df = pd.Series(V, name="viterbi")

        df_bin = write_bin_table(
            X.Z, bins, viterbi_df, [*states.state_names], P, outname=f"{outname}.bin.xz"
        )

    if output["output_rle"]:
        df_rle = get_rle(df_bin, states, init["run_penalty"])
        write_est_runs(df_rle, outname=f"{outname}.rle.xz")

    if output["output_snp"]:
        df_snp = write_snp_table(
            data=df, G=X.PG, Z=X.Z, P=P, gt_mode=gt_mode, outname=f"{outname}.snp.xz"
        )

    if output["output_rsim"]:
        df_pred = pred_sims(
            trans=pars.trans,
            emissions=X.emissions,
            beta=X.beta,
            alpha0=pars.alpha0,
            n=X.n,
            states=P.states,
            n_sims=n_post_replicates,
            decode=not O.est_inbreeding,
            keep_loc=keep_loc,
        )
        df_pred["chrom"] = [P.diplo_chroms[i] for i in df_pred.chrom.values]
        if est["est_inbreeding"]:
            l = [*states.state_names]
            df_pred["state"] = [l[i] for i in df_pred.state.values]
        else:
            df_pred["state"] = [states[i] for i in df_pred.state.values]

        if len(X.H.beta) > 0:
            df_pred_hap = pred_sims(
                trans=pars.trans_hap,
                emissions=X.H.emissions,
                beta=X.H.beta,
                alpha0=pars.alpha0_hap,
                n=X.H.n,
                states=states,
                n_sims=n_post_replicates,
                keep_loc=keep_loc,
                decode=False,
            )
            df_pred_hap["chrom"] = [P.haplo_chroms[i] for i in df_pred_hap.chrom.values]
            df_pred_hap["state"] = [states[i] for i in df_pred_hap.state.values]
            df_pred = pd.concat((df_pred, df_pred_hap))

        if keep_loc:
            df_pred = df_pred[["state", "chrom", "start", "end", "len", "it"]]
        else:
            df_pred = df_pred[["state", "chrom", "len", "it"]]

        write_sim_runs(df_pred, outname=f"{outname}.res.xz")
        resampling_pars(df_pred).to_csv(
            f"{outname}.res2.xz", index=True, float_format="%.6f"
        )

    return
