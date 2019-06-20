"""basic command line interfaces
"""
import argparse
from .admixfrog import run_admixfrog
from pprint import pprint, pformat
from .bam import process_bam
from .rle import get_rle
from .vcf import load_pop_file, vcf_to_ref, load_random_read_samples, vcf_to_sample
from os.path import isfile
import admixfrog
from .log import log_, setup_log
import pdb

# set up a list of option groups to keep things organized

#population comparison options
POP_OPTIONS = ["cont_id", 
               "infile", 
               "ref_files", 
               "sex", 
               "states", 
               "ancestral"]


#generate infile from vcf, bam or geno
INFILE_OPTIONS = [
    "bamfile",
    "deam_cutoff",
    "minmapq",
    "force_infile",
    "length_bin_size",
    "tstv",
    "alleles",
    "vcfgt",
    "sample_id",
]

#generate reference from vcf or geno
REFFILE_OPTIONS = [
    "vcf_file",
    "pop_file",
    "rec_file",
    "rec_rate",
    "pos_id",
    "map_id",
    "chroms",
    "force_ref",
]

#algorithm control options
ALGORITHM_OPTIONS = [
    "bin_size",
    "autosomes_only",
    "downsample",
    "gt_mode",
    "ll_tol",
    "max_iter",
    "n_post_replicates",
    "pos_mode",
    "prior",
    "split_lib",
]

def add_output_parse_group(parser):
    g = parser.add_argument_group(
        "output name and files to be generated",
        """By default, all files are generated. However, if any of 
                                  the --no-* options are used to disable specific files
                                  """,
    )
    g.add_argument(
        "--out",
        "-o",
        default="admixfrog",
        help="""Output file path (without extensions)""",
    )
    g.add_argument(
        "--no-rle",
        action="store_false",
        default=True,
        dest="output_rle",
        help="""Disabble Estimating  runs and writeing to file with extension .rle.xz""",
    )
    g.add_argument(
        "--no-snp",
        action="store_false",
        default=True,
        dest="output_snp",
        help="""Disable writing posterior genotype likelihood to file with extension .snp.xz""",
    )
    g.add_argument(
        "--no-bin",
        action="store_false",
        default=True,
        dest="output_bin",
        help="""Disable writing posterior states to file with extension .bin.xz""",
    )

    g.add_argument(
        "--no-cont",
        action="store_false",
        default=True,
        dest="output_cont",
        help="""Disable writing contamination estimates to file with extension .bin.xz""",
    )
    g.add_argument(
        "--no-rsim",
        action="store_false",
        default=True,
        dest="output_rsim",
        help="""Disable writing posterior simulations of runs 
            to file with extension .res.xz""",
    )
    g.add_argument(
        "--no-pars",
        action="store_false",
        default=True,
        dest="output_pars",
        help="""Disable writing parameters 
            to file with extension .pars.yaml""",
    )


def add_infile_parse_group(parser):
    g = parser.add_argument_group("bam parsing")
    g.add_argument(
        "--bamfile",
        "--bam",
        help="""Bam File to process. Choose this or infile. 
                   The resulting input file will be writen in {out}.in.xz,
                   so it doesn't need to be regenerated.
                   If the input file exists, an error is generated unless
                   --force-infile is set
                   """,
    )
    g.add_argument("--force-infile", "--force-bam", default=False, action="store_true")
    # g.add_argument("--bedfile", "--bed",
    #               help="Bed file with anc/der allele to restrict to")
    g.add_argument(
        "--deam-cutoff",
        type=int,
        default=3,
        help="""reads with deamination in positions < deam-cutoff are
                   considered separately""",
    )
    g.add_argument(
        "--minmapq",
        type=int,
        default=25,
        help="""reads with mapq < MINMAPQ are removed""",
    )
    g.add_argument(
        "--length-bin-size",
        type=int,
        default=None,
        help="""if set, reads are binned by length for contamination estimation""",
    )


def add_rle_parse_group(parser):
    g = parser.add_argument_group("call introgressed fragments")
    g.add_argument(
        "--run-penalty",
        type=float,
        default=0.2,
        help="""penalty for runs. Lower value means runs are called more
        stringently (default 0.2)""",
    )


def add_geno_parse_group(parser):
    g = parser.add_argument_group("""geno (Eigenstrat/Admixtools/Reich) format
                                  parser options""")
    g.add_argument(
        "--gfile",
        help="""geno file name (without extension, expects .snp/.ind/.geno
        files). Only reads binary format for now""",
    )
    g.add_argument(
        "--gtarget",
        help="""
        id of target individual from geno file
        """
    )
    g.add_argument(
        "--gpops",
        help="""
        target populations. These should be either entries in the third column
        of the .ind file, or a combination thereof. If two populations YRI and
        SAN are defined, a syntax of the form 
        AFR=YRI+SAN can be used to define a reference `AFR` that is a
        combination of the wo.
        """,
        nargs="*"
    )
    g.add_argument(
        "--guess-ploidy",
        action='store_true',
        default=False,
        help="""guess ploidy of individuals (use if e.g. random read sample
        inds are present)
        """
    )


def add_ref_parse_group(parser):
    g = parser.add_argument_group("creating reference file")
    g.add_argument(
        "--vcf-file",
        "--vcf",
        help="""VCF File to process. Choose this or reffile. 
                   The resulting ref file will be writen as {out}.ref.xz,
                   so it doesn't need to be regenerated.
                   If the input file exists, an error is generated unless
                   --force-ref is set
                   """,
    )
    g.add_argument(
        "--pop-file", default=None, help="""Population assignments (yaml format)"""
    )
    g.add_argument(
        "--rec-file",
        "--rec",
        help="""Recombination rate file. Modelled after 
        https://www.well.ox.ac.uk/~anjali/AAmap/
        If file is split by chromosome, use {CHROM} as 
        wildcards where the chromosome id will be included
        """,
    )
    g.add_argument(
        "--rec-rate",
        default=1e-8,
        type=float,
        help="""Constant recombination rate (per generation per base-pair)""",
    )
    g.add_argument(
        "--pos-id",
        default="Physical_Pos",
        help="""column name for position (default: Physical_Pos)
        """,
    )
    g.add_argument(
        "--map-id",
        default="AA_Map",
        help="""column name for genetic map (default: AA_Map)
        """,
    )
    parser.add_argument(
        "--chroms",
        "--chromosome-files",
        default="1-22,X",
        help="""The chromosomes to be used in vcf-mode.
        """,
    )
    g.add_argument("--force-ref", "--force-vcf", default=False, action="store_true")


def add_estimation_parse_group(P):
    parser = P.add_argument_group("""options that control estimation of model
                                  parameters""")
    parser.add_argument(
        "--dont-est-contamination",
        action="store_false",
        dest="est_contamination",
        default=True,
        help="""Don't estimate contamination (default do)""",
    )
    parser.add_argument(
        "--est-error",
        action="store_true",
        default=False,
        help="""estimate sequencing error per rg""",
    )
    parser.add_argument(
        "--freq-contamination",
        "--fc",
        type=int,
        default=1,
        help="""update frequency for contamination/error (default 1)""",
    )
    parser.add_argument(
        "--est-F",
        "-f",
        action="store_true",
        default=False,
        help="""Estimate F (distance from ref, default False)""",
    )
    parser.add_argument(
        "--est-tau",
        "-tau",
        action="store_true",
        default=False,
        help="""Estimate tau (population structure in references)""",
    )
    parser.add_argument(
        "--freq-F",
        "--f",
        type=int,
        default=1,
        help="""update frequency for F (default 1)""",
    )
    parser.add_argument(
        "--est-inbreeding",
        "-I",
        default=False,
        action="store_true",
        help="""allow  haploid (i.e. inbreed) stretches. Experimental""",
    )
    parser.add_argument(
        "--F0",
        nargs="*",
        type=float,
        default=0.5,
        help="initial F (should be in [0;1]) (default 0)",
    )
    parser.add_argument(
        "--tau0",
        nargs="*",
        type=float,
        default=0,
        # help="initial tau (should be in [0;1]) (default 1), at most 1 per source",
        help="initial log-tau (default 0), at most 1 per source",
    )
    parser.add_argument(
        "--e0", "-e", type=float, default=1e-2, help="initial error rate"
    )
    parser.add_argument(
        "--c0", "-c", type=float, default=1e-2, help="initial contamination rate"
    )


def add_base_parse_group(P):
    parser = P.add_argument_group("options that control the algorithm behavior")
    parser.add_argument(
        "--gt-mode",
        "--gt",
        help="""Assume genotypes are known.
        """,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-b",
        "--bin-size",
        type=float,
        default=10000,
        help="""Size of bins. By default, this is given in 1e-8 cM, so that the unit is
        approximately the same for runs on physical / map positions""",
    )
    parser.add_argument(
        "--prior",
        "-p",
        type=float,
        default=None,
        help="""Prior of reference allele frequencies. If None (default, recommended), this is 
        estimated from the data
        
        This number is added to both the
        ref and alt allele count for each reference, to reflect the uncertainty in allele
        frequencies from a sample. If references are stationary with size 2N, this is
        approximately  [\sum_i^{2N}(1/i) 2N]^{-1}.
          """,
    )
    parser.add_argument(
        "-P",
        "--pos-mode",
        default=False,
        action="store_true",
        help="""Instad of recombination distances, use physical distances for binning""",
    )
    parser.add_argument(
        "--max-iter",
        "-m",
        type=int,
        default=1000,
        help="""maximum number of iterations""",
    )
    parser.add_argument(
        "--ll-tol", type=float, default=1e-2, help="""stop EM when DeltaLL < ll-tol"""
    )
    parser.add_argument(
        "--dont-split-lib",
        action="store_false",
        dest="split_lib",
        default=True,
        help="""estimate one global contamination parameter (default: one per read
        group)""",
    )
    parser.add_argument(
        "--autosomes-only",
        action="store_true",
        default=False,
        help="Only run autosomes",
    )
    parser.add_argument(
        "--downsample",
        type=float,
        default=1.0,
        help="downsample coverage to a proportion of reads",
    )
    parser.add_argument(
        "--init-guess",
        nargs="*",
        help="""init transition so that one state is favored. should be a 
        state in --state-ids """,
    )


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
        "--vcfgt",
        "--vcf-gt",
        "--vcf-infile",
        help="""VCF input file. To generate input format for admixfrog in genotype mode, use this.
        """,
    )
    parser.add_argument(
        "--sample-id",
        help="""sample id for vcf mode.
        """,
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
    add_infile_parse_group(parser)
    args = vars(parser.parse_args())

    logger.info(pformat(args))
    force_infile = args.pop("force_infile")
    if isfile(args["outfile"]) and not force_infile:
        raise ValueError("""infile exists. set --force-infile to regenerate""")

    if "vcfgt" in args and args["vcfgt"] is not None:
        vcf_to_sample(
            outfile=args["outfile"],
            vcf_file=args["vcfgt"],
            ref_file=args["ref"],
            chroms=args["chroms"],
            random_read=args["random_read_sample"],
            sample_id=args["sample_id"],
        )
    else:
        del args["vcfgt"]
        del args["sample_id"]
        del args["chroms"]
        process_bam(**args)


def do_rle():
    logger = setup_log()
    parser = argparse.ArgumentParser(description="Do RLE encoding")
    parser.add_argument(
        "--outfile", "--out", required=True, help="output file name (xz-zipped)"
    )
    parser.add_argument(
        "--infile", "--in", required=True, help="input file name *.bin.xz"
    )
    add_rle_parse_group(parser)
    args = parser.parse_args()
    logger.info(pformat(args))
    import pandas as pd

    dtype_ = dict(chrom="category")
    data = pd.read_csv(args.infile, dtype=dtype_)
    data.chrom.cat.reorder_categories(pd.unique(data.chrom), inplace=True)
    states = list(data.columns)[6:]
    homo = [s for s in states if sum(s in ss for ss in states) > 1]

    rle = get_rle(data, homo, args.run_penalty)
    rle.to_csv(args.outfile, float_format="%.6f", index=False, compression="xz")


def do_ref():
    logger = setup_log()
    parser = argparse.ArgumentParser(description="create reference file from vcf")
    parser.add_argument(
        "--outfile", "--out", required=True, help="output file name (xz-zipped)"
    )
    parser.add_argument(
        "--states",
        nargs="*",
        help="""Populations to be parsed. Need to be either present in the pop-file,
        or, as single sample, in the vcf-file. If unset, all clusters defined in the 
        pop-file will be used
        """,
    )

    add_ref_parse_group(parser)
    args = parser.parse_args()
    logger.info(pformat(args))

    pop2sample = load_pop_file(args.pop_file, args.states)
    random_read_samples = load_random_read_samples(args.pop_file)
    logger.debug(pformat(random_read_samples))
    vcf_to_ref(
        args.outfile,
        args.vcf_file,
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
        "--states",
        "--state-ids",
        nargs="*",
        default=["AFR", "VIN", "DEN"],
        help="""the allowed sources. The target will be made of a mix of all homozygous
        and heterozygous combinations of states. More than 4 or 5 sources have not been
        tested and are not recommended. Must be present in the ref file
        """,
    )


    parser.add_argument(
        "--cont-id",
        "--cont",
        default="AFR",
        help="""the source of contamination. Must be specified in ref file""",
    )

    parser.add_argument(
        "--ancestral",
        "-a",
        type=str,
        default=None,
        help="""Outgroup population with the ancestral allele. By default, assume
        ancestral allele is unknown
        """,
    )
    parser.add_argument(
        "--n-post-replicates",
        type=int,
        default=100,
        help="""Number of replicates that are sampled from posterior. Useful for
        parameter estimation and bootstrapping
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

    add_infile_parse_group(parser)
    add_estimation_parse_group(parser)
    add_base_parse_group(parser)
    add_ref_parse_group(parser)
    add_rle_parse_group(parser)
    add_output_parse_group(parser)

    from . import __version__

    args = parser.parse_args()
    V = vars(args)

    output_options = dict()
    filter_options = dict()
    init_pars = dict()
    est_options = dict()
    infile_pars = dict()
    reffile_pars = dict()
    algo_pars = dict()

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
            infile_pars[k] = V.pop(k)
        elif k in REFFILE_OPTIONS:
            reffile_pars[k] = V.pop(k)
        elif k in ALGORITHM_OPTIONS:
            algo_pars[k] = V.pop(k)
    V["output"] = output_options
    V["filter"] = filter_options
    V["est"] = est_options
    V["init"] = init_pars
    V["algorithm"] = algo_pars

    log_.info(
        "running admixfrog %s with the following arguments:\n%s ",
        __version__,
        pformat(V),
    )

    if reffile_pars["vcf_file"] is not None:
        if V["ref_files"] is not None:
            raise ValueError("cant specify ref and vcf input")
        V["ref_files"] = [V["out"] + ".ref.xz"]
        if isfile(V["ref_files"][0]) and not reffile_pars["force_ref"]:
            raise ValueError(
                """ref-file exists. Use this or set --force-ref to 
                regenerate the file"""
            )
        log_.info("creating ref from vcf file")

        pop2sample = load_pop_file(reffile_pars["pop_file"], V["states"])
        random_read_samples = load_random_read_samples(reffile_pars["pop_file"])
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

    if infile_pars["bamfile"] is not None:
        if V["infile"] is not None:
            raise ValueError("cant specify csv and bam input")
        V["infile"] = V["out"] + ".in.xz"
        if isfile(V["infile"]) and not infile_pars["force_infile"]:
            raise ValueError(
                """infile exists. Use this or set --force-infile to 
                             regenerate"""
            )
        log_.info("creating input from bam file")
        process_bam(
            outfile=V["infile"],
            bamfile=infile_pars["bamfile"],
            ref=V["ref_files"][0],
            # bedfile=V.pop('bedfile'),
            deam_cutoff=infile_pars["deam_cutoff"],
            length_bin_size=infile_pars["length_bin_size"],
        )

    if "infile" not in V or V["infile"] is None:
        raise ValueError("require infile, specify --infile or --bam or --vcfgt")
    if "ref_files" not in V or V["ref_files"] is None:
        raise ValueError("require ref, specify with --ref or create using --vcf")

    out = V.pop("out")

    from . import __version__

    log_.info("admixfrog %s", __version__)
    del V["algorithm"]
    bins, snps, cont, pars, rle, res = run_admixfrog(**V, **algo_pars)

    res.to_csv("%s.res.xz" % out, float_format="%.6f", index=False, compression="xz")
    rle.to_csv("%s.rle.xz" % out, float_format="%.6f", index=False, compression="xz")
    bins.to_csv("%s.bin.xz" % out, float_format="%.6f", index=False, compression="xz")
    cont.to_csv("%s.cont.xz" % out, float_format="%.6f", index=False, compression="xz")
    with open("%s.pars.yaml" % out, "wt") as f:
        f.write(pars)
    snps.to_csv("%s.snp.xz" % out, float_format="%.6f", index=False, compression="xz")


def profile():
    import cProfile
    import pstats

    cProfile.runctx("run()", globals(), locals(), filename="profile.txt")
    pstats.Stats("profile.txt").strip_dirs().sort_stats(1).print_stats(200)
