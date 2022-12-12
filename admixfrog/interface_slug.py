"""basic command line interfaces
"""
import argparse
import os
from pprint import pprint, pformat
import admixfrog
import numpy as np
import logging
from os.path import isfile

from .slug import run_admixslug
from .utils.bam import process_bam, process_bam2
from .utils.vcf import (
    load_pop_file,
    vcf_to_ref,
    load_random_read_samples,
    vcf_to_sample,
)
from .utils.log import setup_log
from .options import INFILE_OPTIONS, REFFILE_OPTIONS, GENO_OPTIONS
from .options import ALGORITHM_OPTIONS_SLUG
from .options import add_target_file_options, add_rle_options
from .options import add_filter_options
from .options import add_geno_options, add_ref_options
from .options import add_pop_options, add_base_options_slug
from .options import add_output_options_slug, add_estimation_options_slug


def run_sfs():
    parser = argparse.ArgumentParser(
        description="Infer sfs and contamination from low-coverage and contaminated genomes"
        #        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + admixfrog.__version__
    )

    parser.add_argument(
        "--target-file",
        "--infile",
        "--in",
        help="""Sample input file (csv). Contains individual specific data, obtained from
                        a bam file.

                    - Fields are chrom, pos, map, lib, tref, talt"
                        - chrom: chromosome
                        - pos : physical position (int)
                        - map : rec position (float)
                        - lib : read group. Any string, same string assumes same contamination 
                        - tref : number of reference reads observed
                        - talt: number of alt reads observed""",
    )
    parser.add_argument(
        "--ref",
        "--ref-file",
        default=None,
        action="append",
        dest="ref_files",
        help="""refernce input file (csv). 
                    - Fields are chrom, pos, ref, alt, map, X_alt, X_ref
                        - chrom: chromosome
                        - pos : physical position (int)
                        - ref : refrence allele
                        - alt : alternative allele
                        - map : rec position (float)
                        - X_alt, X_ref : alt/ref alleles from any number of sources / contaminant populations.
                        these are used later in --cont-id and --state-id flags
                        """,
    )

    parser.add_argument(
        "--force-ref", "--force-vcf", default=False, action="store_true"
    )

    parser.add_argument(
        "--seed", help="random number generator seed for resampling", default=None
    )
    # removed for now, as it is default and requierd
    # parser.add_argument(
    #    "--bin-reads",
    #    default=False,
    #    action="store_true",
    #    help="""Input file has info for binning reads. If false,
    #                    reads are grouped by the `lib` column. Otherwise, it
    #                    uses a deam and length column to bin itself""",
    # )

    add_target_file_options(parser)
    add_estimation_options_slug(parser)
    add_base_options_slug(parser)
    add_output_options_slug(parser)
    add_pop_options(parser)
    add_filter_options(parser)

    from . import __version__

    args = parser.parse_args()
    V = vars(args)

    seed = V.pop("seed")
    if seed is not None:
        seed = int(seed)
        logging.info(f"numpy random seed:{seed}")
    np.random.seed(seed)

    # this reorganizes options into sensible groups
    output_options = dict()
    filter_options = dict()
    init_pars = dict()
    est_options = dict()
    target_pars = dict()
    reffile_pars = dict()
    algo_pars = dict()
    geno_pars = dict()

    for k in list(V.keys()):
        if k.startswith("output_"):
            output_options[k] = V.pop(k)
        elif k.startswith("filter_"):
            filter_options[k] = V.pop(k)
        elif k.startswith("est_") or k.startswith("freq_"):
            est_options[k] = V.pop(k)
        elif k in ["init_guess", "F0", "e0", "c0", "tau0", "run_penalty"]:
            init_pars[k] = V.pop(k)
        elif k in INFILE_OPTIONS:
            target_pars[k] = V.pop(k)
        elif k in REFFILE_OPTIONS:
            reffile_pars[k] = V.pop(k)
        elif k in ALGORITHM_OPTIONS_SLUG:
            algo_pars[k] = V.pop(k)
        elif k in GENO_OPTIONS:
            geno_pars[k] = V.pop(k)

    del V["random_read_samples"]

    # all the stuff that will get passed to admixfrog.run_admixfrog
    V["output"] = output_options
    V["filter"] = filter_options
    V["est"] = est_options
    V["init"] = init_pars

    # V["algorithm"] = algo_pars
    # V["geno"] = geno_pars
    # V['ref'] = reffile_pars
    setup_log(filename=f"{args.outname}.log")

    logging.info(
        "running admixslug %s with the following arguments:\n%s ",
        __version__,
        pformat(V),
    )

    # 1. find ref-file the following precedent is set:
    # . 1.1. directly specified with --ref-files
    # . 1.2. create from --vcf-ref option
    # . 1.3. in geno-format, specified as --gfile

    if V["ref_files"] is None:
        if reffile_pars["vcf_ref"] is None and geno_pars["geno_file"] is not None:
            logging.info("reading geno file \n %s", pformat(geno_pars))
        elif reffile_pars["vcf_ref"] is not None:
            logging.info(
                "creating admixslug reference file from vcf\n %s", pformat(reffile_pars)
            )
            V["ref_files"] = [V["outname"] + ".ref.xz"]
            if isfile(V["ref_files"][0]) and not reffile_pars["force_ref"]:
                raise ValueError(
                    """ref-file exists. Use this or set --force-ref to 
                    regenerate the file"""
                )
            logging.info("creating ref from vcf file")

            pop2sample = load_pop_file(reffile_pars["state_file"], V["states"])
            random_read_samples = load_random_read_samples(reffile_pars["state_file"])
            logger.debug(pformat(random_read_samples))
            vcf_to_ref(
                V["ref_files"][0],
                reffile_pars["vcf_file"],
                reffile_pars["rec_file"],
                pop2sample,
                random_read_samples,
                reffile_pars["pos_id"],
                reffile_pars["map_id"],
                reffile_pars["rec_rate"],
                reffile_pars["chroms"],
            )
        else:
            raise ValueError(
                "no ref defined (set either --ref-files, --vcf-ref or --gfile)"
            )

    # 2. find genome to-be decoded
    #  2.1  directly specified with --target-file
    #  2.2  specified with --gfile and --target
    #  2.3  specified with --bam
    #  2.4  specified with --vcfgt and --target
    if target_pars["target_file"] is None:
        if target_pars["bamfile"] is not None:
            target_pars["target_file"] = V["outname"] + ".in.xz"
            if (
                isfile(target_pars["target_file"])
                and not target_pars["force_target_file"]
            ):
                raise ValueError(
                    """target_file exists. Use this or set --force-target-file to 
                                 regenerate"""
                )
            logging.info("creating input from bam file")
            process_bam2(
                outfile=target_pars["target_file"],
                bamfile=target_pars["bamfile"],
                ref=V["ref_files"][0],
                deam_cutoff=target_pars["deam_cutoff"],
                length_bin_size=target_pars["length_bin_size"],
                chroms=reffile_pars["chroms"],
            )
        elif target_pars["vcfgt"] is not None and target_pars["target"] is not None:
            target_pars["target_file"] = V["outname"] + ".in.xz"
            if (
                isfile(target_pars["target_file"])
                and not target_pars["force_target_file"]
            ):
                raise ValueError(
                    """target_file exists. Use this or set --force-target_file to 
                                 regenerate"""
                )
            logging.info("creating input from bam file")
            process_bam(
                outfile=target_pars["target_file"],
                bamfile=target_pars["bamfile"],
                ref=V["ref_files"][0],
                deam_cutoff=target_pars["deam_cutoff"],
                length_bin_size=target_pars["length_bin_size"],
            )
        elif geno_pars["geno_file"] is not None and target_pars["target"] is not None:
            pass
        else:
            raise ValueError(
                "no sample defined (set either --target-file --bam --vcfgt or --gfile must be set"
            )

    if target_pars["target"] is None:
        target_pars["target"] = os.path.basename(target_pars["target_file"]).split("_")[
            0
        ]

    from . import __version__

    logging.info(f"admixslug {__version__}")

    # run stuff
    run_admixslug(
        **V,
        **algo_pars,
        **geno_pars,
        target=target_pars["target"],
        target_file=target_pars["target_file"],
    )


def profile():
    import cProfile
    import pstats

    cProfile.runctx("run_sfs()", globals(), locals(), filename="profile.txt")
    pstats.Stats("profile.txt").strip_dirs().sort_stats(1).print_stats(200)
    import cProfile
