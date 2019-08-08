"""basic command line interfaces
"""
import argparse
import random
from pprint import pprint, pformat
import admixfrog

from .admixfrog import run_admixfrog
from .bam import process_bam
from .rle import get_rle
from .vcf import load_pop_file, vcf_to_ref, load_random_read_samples, vcf_to_sample
from os.path import isfile
from .log import log_, setup_log
from .options import POP_OPTIONS, INFILE_OPTIONS, REFFILE_OPTIONS, ALGORITHM_OPTIONS, GENO_OPTIONS
from .options import add_output_options, add_target_file_options, add_rle_options
from .options import add_geno_options, add_ref_options, add_estimation_options
from .options import add_base_options, add_pop_options


def bam():
    """create input file from bam/vcf"""
    logger = setup_log()
    parser = argparse.ArgumentParser(description="Parse bam file for admixfrog")
    parser.add_argument(
        "--outfile", "--out", required=True, help="output file name (xz-zipped)"
    )
    parser.add_argument(
        "--ref",
        "--ref-file",
        help="""reference input file (csv). 
                    - Fields are chrom, pos, ref, alt, map, X_alt, X_ref
                        - chrom: chromosome
                        - pos : physical position (int)
                        - ref : refrence allele
                        - alt : alternative allele
                        - map : rec position (float)
                        - X_alt, X_ref : alt/ref alleles from any number of sources / contaminant populations.
                        these are used later in --cont-id and --states flags
                        """,
        required=True,
    )
    parser.add_argument(
        "--random-read-sample",
        default=False,
        action="store_true",
        help="""At each position, just use a single read (bam-mode), or assume
        that the genotype reflects a single randomly sampled read
        """,
    )
    parser.add_argument(
        "--chroms",
        "--chromosome-files",
        default="1-22,X",
        help="""The chromosomes to be used in vcf-mode.
        """,
    )
    add_target_file_options(parser)
    args = vars(parser.parse_args())

    logger.info(pformat(args))
    force_target_file = args.pop("force_target_file")
    if isfile(args["outfile"]) and not force_target_file:
        raise ValueError("""target file exists. set --force-target-file to regenerate""")

    if "vcfgt" in args and args["vcfgt"] is not None:
        vcf_to_sample(
            outfile=args["outfile"],
            vcf_file=args["vcfgt"],
            ref_file=args["ref"],
            chroms=args["chroms"],
            random_read=args["random_read_sample"],
            sample_id=args["target"],
        )
    else:
        del args["vcfgt"]
        del args["target"]
        del args["chroms"]
        process_bam(**args)


def do_rle():
    logger = setup_log()
    parser = argparse.ArgumentParser(description="Do RLE encoding")
    parser.add_argument(
        "--outfile", "--out", required=True, help="output file name (xz-zipped)"
    )
    parser.add_argument(
        "--target-file", "--in", required=True, help="input file name *.bin.xz"
    )
    add_rle_options(parser)
    args = parser.parse_args()
    logger.info(pformat(args))
    import pandas as pd

    dtype_ = dict(chrom="category")
    data = pd.read_csv(args.target_file, dtype=dtype_)
    data.chrom.cat.reorder_categories(pd.unique(data.chrom), inplace=True)
    states = list(data.columns)[6:]
    homo = [s for s in states if sum(s in ss for ss in states) > 1]

    rle = get_rle(data, homo, args.run_penalty)
    rle.to_csv(args.outfile, float_format="%.6f", index=False, compression="xz")


def do_ref():
    """subprogram to create reference file
    """
    logger = setup_log()
    parser = argparse.ArgumentParser(description="create reference file from vcf")
    parser.add_argument(
        "--outfile", "--out", required=True, help="output file name (xz-zipped)"
    )

    add_pop_options(parser)
    add_ref_options(parser)
    args = parser.parse_args()
    logger.info(pformat(args))

    pop2sample = load_pop_file(args.state_file, args.states)
    random_read_samples = load_random_read_samples(args.state_file)
    logger.debug(pformat(random_read_samples))
    vcf_to_ref(
        args.outfile,
        args.vcf_ref,
        args.rec_file,
        pop2sample,
        random_read_samples,
        args.pos_id,
        args.map_id,
        rec_rate=args.rec_rate,
        chroms=args.chroms,
    )


def run():
    logger = setup_log()

    parser = argparse.ArgumentParser(
        description="Infer admixture frogments from low-coverage and contaminated genomes"
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
        "--map-col",
        "--map-column",
        default='map',
        help="""Name of reference-file column that contains the 
        recombination map. (default: 'map')
                        """,
    )
    parser.add_argument(
        "--filter-delta",
        type=float,
        help="""only use sites with allele frequency difference bigger than DELTA (default off)""",
    )
    parser.add_argument(
        "--filter-pos",
        type=int,
        help="""greedily prune sites to be at least POS positions apart""",
    )
    parser.add_argument(
        "--filter-map",
        type=float,
        help="""greedily prune sites to be at least MAP recombination distance apart""",
    )
    parser.add_argument(
        "--male",
        dest="sex",
        action="store_const",
        const="m",
        default=None,
        help="Assumes haploid X chromosome. Default is guess from coverage. currently broken",
    )
    parser.add_argument(
        "--female",
        dest="sex",
        action="store_const",
        const="f",
        default=None,
        help="Assumes diploid X chromosome. Default is guess from coverage",
    )

    parser.add_argument("--seed", help="random number generator seed for resampling",
                        default=None)

    add_target_file_options(parser)
    add_geno_options(parser)
    add_estimation_options(parser)
    add_base_options(parser)
    add_ref_options(parser)
    add_rle_options(parser)
    add_output_options(parser)
    add_pop_options(parser)

    from . import __version__

    args = parser.parse_args()
    V = vars(args)

    random.seed(V.pop("seed"))

    #this reorganizes options into sensible groups
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
        elif k in ALGORITHM_OPTIONS:
            algo_pars[k] = V.pop(k)
        elif k in GENO_OPTIONS:
            geno_pars[k] = V.pop(k)

    #all the stuff that will get passed to admixfrog.run_admixfrog
    V["output"] = output_options
    V["filter"] = filter_options
    V["est"] = est_options
    V["init"] = init_pars

    #V["algorithm"] = algo_pars
    #V["geno"] = geno_pars
    #V['ref'] = reffile_pars

    log_.info(
        "running admixfrog %s with the following arguments:\n%s ",
        __version__,
        pformat(V),
    )

    #1. find ref-file the following precedent is set:
    #. 1.1. directly specified with --ref-files
    #. 1.2. create from --vcf-ref option
    #. 1.3. in geno-format, specified as --gfile


    if V["ref_files"] is None:
        if reffile_pars['vcf_ref'] is None and geno_pars['geno_file'] is not None:
            log_.info(
                "reading geno file \n %s", 
                pformat(geno_pars),
            )
        elif reffile_pars['vcf_ref'] is not None:
            log_.info(
                "creating admixfrog reference file from vcf\n %s", 
                pformat(reffile_pars),
            )
            V["ref_files"] = [V["outname"] + ".ref.xz"]
            if isfile(V["ref_files"][0]) and not reffile_pars["force_ref"]:
                raise ValueError(
                    """ref-file exists. Use this or set --force-ref to 
                    regenerate the file"""
                )
            log_.info("creating ref from vcf file")

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
            raise ValueError("no ref defined (set either --ref-files, --vcf-ref or --gfile)")

    #2. find genome to-be decoded
    #  2.1  directly specified with --target-file
    #  2.2  specified with --gfile and --target
    #  2.3  specified with --bam
    #  2.4  specified with --vcfgt and --target
    if target_pars["target_file"] is None:
        if target_pars["bamfile"] is not None:
            target_pars["target_file"] = V["outname"] + ".in.xz"
            if isfile(target_pars["target_file"]) and not target_pars["force_target_file"]:
                raise ValueError(
                    """target_file exists. Use this or set --force-target-file to 
                                 regenerate"""
            )
            log_.info("creating input from bam file")
            process_bam(
                outfile=target_pars["target_file"],
                bamfile=target_pars["bamfile"],
                ref=V["ref_files"][0],
                deam_cutoff=target_pars["deam_cutoff"],
                length_bin_size=target_pars["length_bin_size"],
            )
        elif target_pars["vcfgt"] is not None and target_pars['target'] is not None:
            target_pars["target_file"] = V["outname"] + ".in.xz"
            if isfile(target_pars["target_file"]) and not target_pars["force_target_file"]:
                raise ValueError(
                    """target_file exists. Use this or set --force-target_file to 
                                 regenerate"""
            )
            log_.info("creating input from bam file")
            process_bam(
                outfile=target_pars["target_file"],
                bamfile=target_pars["bamfile"],
                ref=V["ref_files"][0],
                deam_cutoff=target_pars["deam_cutoff"],
                length_bin_size=target_pars["length_bin_size"],
            )
        elif geno_pars["geno_file"] is not None and target_pars['target'] is not None:
            pass
        else:
            raise ValueError("no sample defined (set either --target-file --bam --vcfgt or --gfile must be set")


    from . import __version__

    log_.info("admixfrog %s", __version__)
    
    #run stuff
    run_admixfrog(**V, **algo_pars, **geno_pars, 
                  target=target_pars['target'],
                  target_file=target_pars['target_file'])



def profile():
    import cProfile
    import pstats

    cProfile.runctx("run()", globals(), locals(), filename="profile.txt")
    pstats.Stats("profile.txt").strip_dirs().sort_stats(1).print_stats(200)
