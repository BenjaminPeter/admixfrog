import numpy as np
import pandas as pd
import itertools
import yaml
import logging
from collections import Counter, defaultdict
from ..utils.input import load_read_data, load_ref, filter_ref
from ..utils.output import write_pars_table
from ..utils.output_frog import write_bin_table, write_cont_table_frog
from ..utils.output_frog import write_snp_table, write_est_runs, write_sim_runs
from ..utils.output_slug import write_cont_table_slug
from ..utils.utils import bins_from_bed, data2probs, init_pars, Pars
from ..utils.utils import guess_sex
from ..utils.states import States
from .fwd_bwd import fwd_bwd_algorithm, viterbi, update_transitions
from ..gll.genotype_emissions import update_post_geno, update_Ftau, update_snp_prob
from ..gll.genotype_emissions import update_emissions
from ..gll.read_emissions import update_contamination
from .rle import get_rle
from .decode import pred_sims, resampling_pars
from ..utils.geno_io import read_geno_ref, read_geno

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


def baum_welch(
    P,
    IX,
    pars,
    est_options,
    max_iter=1000,
    ll_tol=1e-1,
    gt_mode=False,
    scale_probs=True,
):
    O = est_options
    alpha0, alpha0_hap, trans, trans_hap, cont, error, F, tau, sex = pars
    gll_mode = not gt_mode
    ll = -np.inf
    n_states = P.S.n_states
    n_hap_states = P.S.n_hap
    n_gt = 3

    # create arrays for posterior, emissions
    Z = np.zeros((IX.n_bins, n_states))  # P(Z | O)
    E = np.ones((IX.n_bins, n_states))  # P(O | Z)

    SNP = np.zeros(
        (IX.n_snps, n_states, n_gt)
    )  # P(O, G | Z), scaled such that  the max is 1
    PG = np.zeros((IX.n_snps, n_states, n_gt))  # P(G Z | O)

    # pointers to the same data, but split by chromosome
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
        SNP,
        P,
        IX,
        cont,
        error,
        F,
        tau,
        O["est_inbreeding"],
        gt_mode,
        scale_probs=scale_probs,
    )  # return value here is SNP

    e_scaling = update_emissions(E, SNP, P, IX, scale_probs=scale_probs)  # P(O | Z)
    logging.info("e-scaling: %s", e_scaling)
    logging.info("s-scaling: %s", s_scaling)
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
        logging.info("[%dk|%dk|%dk|%dk]: iter:%d |p95:%.3f\tLL:%.4f\tΔLL:%.4f" % tpl)
        if ll - old_ll < ll_tol:
            break

        if np.any(np.isnan(Z)):
            raise ValueError("nan observed in state posterior")
        if np.any(np.isnan(E)):
            raise ValueError("nan observed in emissions")

        # update stuff
        if O["est_trans"]:
            trans = update_transitions(
                trans,
                alpha,
                beta,
                gamma,
                emissions,
                n,
                est_inbreeding=O["est_inbreeding"],
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
            logging.debug(trans)
            logging.debug(trans_hap)
            alpha0 = np.linalg.matrix_power(trans, 10000)[0]
            alpha0_hap = np.linalg.matrix_power(trans_hap, 10000)[0]

        logging.info("\t".join(P.S.state_names))
        logging.info("\t".join(["%.3f" % a for a in alpha0]))

        scaling = update_emission_stuff(
            it,
            E,
            P,
            PG,
            SNP,
            Z,
            IX,
            cont,
            error,
            F,
            tau,
            scaling,
            est_options,
            gt_mode,
            scale_probs=scale_probs,
        )

    update_post_geno(PG, SNP, Z, IX)
    pars = Pars(alpha0, alpha0_hap, trans, trans_hap, cont, error, F, tau, sex)

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
    it,
    E,
    P,
    PG,
    SNP,
    Z,
    IX,
    cont,
    error,
    F,
    tau,
    scaling,
    est_options,
    gt_mode,
    scale_probs,
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
            logging.info("stopping contamination updates")

    if cond_Ftau:
        delta = update_Ftau(F, tau, PG, P, IX, est_options)
        if delta < 1e-5:  # when we converged, do not update F
            O["est_F"], O["est_tau"] = False, False
            cond_Ftau, cond_F = False, False
            logging.info("stopping Ftau updates")
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
            scale_probs=scale_probs,
        )

        e_scaling = update_emissions(E, SNP, P, IX, scale_probs=scale_probs)  # P(O | Z)
        logging.info("e-scaling: %s", e_scaling)
        logging.info("s-scaling: %s", s_scaling)
        scaling = e_scaling + s_scaling

    return scaling


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
    df = filter_ref(df, states, ancestral=ancestral, **filter)
    df["rg"] = "rg0"
    if pos_mode:
        df.reset_index("map", inplace=True)
        df["map"] = df.index.get_level_values("pos")
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
    ancestral=None,
    cont_id=None,
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
):
    """
    we have the following possible input files
    1. only a geno file (for target and reference)
    2. custom ref/target (from standalone csv)
    3. geno ref (and target from csv)
    """
    tot_n_snps = 0

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
        ref = load_ref(
            ref_files, states, cont_id, ancestral, autosomes_only, map_col=map_col
        )
        ref = filter_ref(ref, states, ancestral=ancestral, **filter)

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

        if fake_contamination and cont_id:
            """filter for SNP with fake ref data"""
            cont_ref, cont_alt = f"{cont_id}_ref", f"{cont_id}_alt"
            ref = ref[ref[cont_ref] + ref[cont_alt] > 0]

        tot_n_snps = ref.shape[0]

        if pos_mode:
            ref.reset_index("map", inplace=True)
            ref["map"] = ref.index.get_level_values("pos")
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

    if fake_contamination and cont_id:
        cont_ref, cont_alt = f"{cont_id}_ref", f"{cont_id}_alt"
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

    return df, ix, sex, tot_n_snps


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
    haplo_chroms=None,
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

    from pprint import pprint

    pprint(kwargs)

    # numpy config
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

    bins, IX = bins_from_bed(
        df, bin_size=bin_size, sex=sex, snp_mode=snp_mode, haplo_chroms=haplo_chroms
    )
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
        homo_ids=homo_states,
        het_ids=het_states,
        sex=sex,
        bin_size=bin_size,
        **init,
    )

    Z, G, pars, ll, emissions, hemissions, (_, beta, n), (_, bhap, nhap) = baum_welch(
        P, IX, pars, gt_mode=gt_mode, est_options=est, **kwargs
    )

    # output formating from here
    if output["output_pars"]:
        df_pars = write_pars_table(pars, outname=f"{outname}.pars.yaml")

    if output["output_cont"]:
        if ix is None:
            df_cont = write_cont_table_frog(
                df, pars.cont, pars.error, tot_n_snps, outname=f"{outname}.cont.xz"
            )
        else:
            rgs = dict((l, i) for i, l in enumerate(IX.rgs))
            cont = [pars.cont[k] for k in rgs.keys()]

            df_cont = write_cont_table_slug(
                ix,
                rgs,
                cont,
                tot_n_snps,
                se=None,
                outname=f"{outname}.cont.xz",
            )
            ix.to_csv(f"{outname}.ix.xz", index=False)

    if output["output_bin"] or output["output_rle"]:
        viterbi_path = viterbi(pars.alpha0, pars.trans, emissions)
        viterbi_path_hap = viterbi(pars.alpha0_hap, pars.trans_hap, hemissions)
        V = np.array([*states.state_names])[np.hstack(viterbi_path)]
        viterbi_df = pd.Series(V, name="viterbi")

        df_bin = write_bin_table(
            Z, bins, viterbi_df, [*states.state_names], IX, outname=f"{outname}.bin.xz"
        )

    if output["output_rle"]:
        df_rle = get_rle(df_bin, states, init["run_penalty"])
        write_est_runs(df_rle, outname=f"{outname}.rle.xz")

    if output["output_snp"]:
        df_snp = write_snp_table(
            data=df, G=G, Z=Z, IX=IX, gt_mode=gt_mode, outname=f"{outname}.snp.xz"
        )

    if output["output_rsim"]:
        df_pred = pred_sims(
            trans=pars.trans,
            emissions=emissions,
            beta=beta,
            alpha0=pars.alpha0,
            n=n,
            states=states,
            n_sims=n_post_replicates,
            decode=not est["est_inbreeding"],
            keep_loc=keep_loc,
        )
        df_pred["chrom"] = [IX.diplo_chroms[i] for i in df_pred.chrom.values]
        if est["est_inbreeding"]:
            l = [*states.state_names]
            df_pred["state"] = [l[i] for i in df_pred.state.values]
        else:
            df_pred["state"] = [states[i] for i in df_pred.state.values]

        if len(bhap) > 0:
            df_pred_hap = pred_sims(
                trans=pars.trans_hap,
                emissions=hemissions,
                beta=bhap,
                alpha0=pars.alpha0_hap,
                n=nhap,
                states=states,
                n_sims=n_post_replicates,
                keep_loc=keep_loc,
                decode=False,
            )
            df_pred_hap["chrom"] = [
                IX.haplo_chroms[i] for i in df_pred_hap.chrom.values
            ]
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
