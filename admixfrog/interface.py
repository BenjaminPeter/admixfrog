import argparse
from .admixfrog import run_admixfrog
from pprint import pprint, pformat
from .bam import process_bam
from .rle import get_rle
from os.path import isfile
import admixfrog
import logging


def add_bam_parse_group(parser):
    g = parser.add_argument_group("bam parsing")
    g.add_argument(
        "--bamfile",
        "--bam",
        help="""Bam File to process. Choose this or infile. 
                   The resulting input file will be writen in {out}.in.xz,
                   so it doesn't need to be regenerated.
                   If the input file exists, an error is generated unless
                   --force-bam is set
                   """,
    )
    g.add_argument("--force-bam", default=False, action="store_true")
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
        "--length-bin-size",
        type=int,
        default=None,
        help="""if set, reads are binned by length for contamination estimation""",
    )


def add_rle_parse_group(parser):
    g = parser.add_argument_group("bam parsing")
    g.add_argument(
        "--run-penalty",
        type=float,
        default=0.5,
        help="""penalty for runs. Lower value means runs are called more
        stringently (default 0.5)""",
    )


def bam():
    parser = argparse.ArgumentParser(description="Parse bam file for admixfrog")
    parser.add_argument(
        "--outfile", "--out", required=True, help="output file name (xz-zipped)"
    )
    parser.add_argument(
        "--ref",
        "--ref-file",
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
    add_bam_parse_group(parser)
    args = vars(parser.parse_args())
    logging.info(pformat(args))
    force_bam = args.pop("force_bam")
    if isfile(args["outfile"]) and not force_bam:
        raise ValueError("""infile exists. set --force-bam to regenerate""")

    process_bam(**args)


def do_rle():
    parser = argparse.ArgumentParser(description="Do RLE encoding")
    parser.add_argument(
        "--outfile", "--out", required=True, help="output file name (xz-zipped)"
    )
    parser.add_argument(
        "--infile", "--in", required=True, help="input file name *.bin.xz"
    )
    add_rle_parse_group(parser)
    args = parser.parse_args()
    logging.info(pformat(args))
    import pandas as pd

    data = pd.read_csv(args.infile)
    states = list(data.columns)[8:]
    homo = [s for s in states if sum(s in ss for ss in states) > 1]

    rle = get_rle(data, homo, args.run_penalty)
    rle.to_csv(args.outfile, float_format="%.6f", index=False, compression="xz")


def run():
    parser = argparse.ArgumentParser(
        description="Infer admixture frogments from low-coverage and contaminated genomes"
        #        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + admixfrog.__version__
    )

    parser.add_argument(
        "--infile", "--in",
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
        "--ref-file",
        "--ref",
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
        "--out", "-o", required=True, help="""Output file path (without extensions)"""
    )
    parser.add_argument(
        "--state-ids",
        "--states",
        nargs="*",
        default=["AFR", "VIN", "DEN"],
        help="""the allowed sources. The target will be made of a mix of all homozygous
        and heterozygous combinations of states. More than 4 or 5 sources have not been
        tested and are not recommended. Must be present in the ref file, with a few
        additional ones:
        - REF : always reference allele
        - NRE : always non-ref allele
        - SFS : allele frequencies are drawn from a  Beta(prior, prior) distribution
        - UNIF : allele frequencies are drawn from a uniform / Beta(1, 1) distribution
        - HALF : allele frequencies are drawn from a  Beta(0.5, 0.5) distribution
        - ZERO : allele frequencies are drawn from a  Beta(0, 0) distribution
        """,
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
        "-P",
        "--pos-mode",
        default=False,
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
        "--prior",
        "-p",
        type=float,
        default=None,
        help="""Prior of reference allele frequencies. IF None (default), this is 
        estimated from the data
        
        This number is added to both the
        ref and alt allele count for each reference, to reflect the uncertainty in allele
        frequencies from a sample. If references are stationary with size 2N, this is
        approximately  [\sum_i^{2N}(1/i) 2N]^{-1}.
          """,
    )
    parser.add_argument(
        "--dont-est-contamination",
        action="store_false",
        dest="est_contamination",
        default=True,
        help="""Don't estimate contamination (default do)""",
    )
    parser.add_argument(
        "--freq-contamination",
        "--fc",
        type=int,
        default=1,
        help="""update frequency for contamination (default 1)""",
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
        "--cont-id",
        "--cont",
        default="AFR",
        help="""the source of contamination. Must be specified in ref file""",
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
        "--male",
        dest="sex",
        action="store_const",
        const="m",
        default=None,
        help="Assumes haploid X chromosome. Default is guess from coverage",
    )
    parser.add_argument(
        "--female",
        dest="sex",
        action="store_const",
        const="f",
        default=None,
        help="Assumes diploid X chromosome. Default is guess from coverage",
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
        help="downsample coveragem keep only DS",
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
        #help="initial tau (should be in [0;1]) (default 1), at most 1 per source",
        help="initial log-tau (default 0), at most 1 per source",
    )
    parser.add_argument(
        "--e0", "-e", type=float, default=1e-2, help="initial error rate"
    )
    parser.add_argument(
        "--c0", "-c", type=float, default=1e-2, help="initial contamination rate"
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
        "--est-inbreeding",
        "-I",
        default=False,
        action="store_true",
        help="""allow  haploid (i.e. inbreed) stretches""",
    )

    add_bam_parse_group(parser)
    add_rle_parse_group(parser)

    from . import __version__
    logging.info("running admixfrog version %s", __version__)
    args = parser.parse_args()
    V = vars(args)
    logging.info(pformat(V))
    force_bam = V.pop("force_bam")

    if V["infile"] is not None and V["bamfile"] is not None:
        raise ValueError("cant specify csv and bam input")
    # elif V["bamfile"] is not None and V["bedfile"] is None:
    #    raise ValueError("require bed file to create input from bam")
    # if V["bamfile"] is not None and V["bedfile"] is not None:
    if V["bamfile"] is not None:
        V["infile"] = V["out"] + ".in.xz"
        if isfile(V["infile"]) and not force_bam:
            raise ValueError(
                """infile exists. Use this or set --force-bam to 
                             regenerate"""
            )
        print("creating input from bam file")
        process_bam(
            outfile=V["infile"],
            bamfile=V.pop("bamfile"),
            ref=V["ref_file"],
            # bedfile=V.pop('bedfile'),
            deam_cutoff=V.pop("deam_cutoff"),
            length_bin_size=V.pop("length_bin_size"),
        )
    else:
        del V["bamfile"]
        del V["deam_cutoff"]
        del V["length_bin_size"]

    out = V.pop("out")

    bins, snps, cont, pars, rle, res = run_admixfrog(**vars(args))
    # bins.to_csv(f"{out}.bin.xz", float_format="%.6f", index=False, compression="xz")
    # cont.to_csv(f"{out}.cont.xz", float_format="%.6f", index=False, compression="xz")
    # pars.to_csv(f"{out}.pars.xz", float_format="%.6f", index=False, compression="xz")
    # snps.to_csv(f"{out}.snp.xz", float_format="%.6f", index=False, compression="xz")
    res.to_csv("%s.res.xz" % out, float_format="%.6f", index=False, compression="xz")
    rle.to_csv("%s.rle.xz" % out, float_format="%.6f", index=False, compression="xz")
    bins.to_csv("%s.bin.xz" % out, float_format="%.6f", index=False, compression="xz")
    cont.to_csv("%s.cont.xz" % out, float_format="%.6f", index=False, compression="xz")
    pars.to_csv("%s.pars.xz" % out, float_format="%.6f", index=False, compression="xz")
    snps.to_csv("%s.snp.xz" % out, float_format="%.6f", index=False, compression="xz")


def profile():
    import cProfile
    import pstats

    cProfile.runctx("run()", globals(), locals(), filename="profile.txt")
    pstats.Stats("profile.txt").strip_dirs().sort_stats(1).print_stats(200)
