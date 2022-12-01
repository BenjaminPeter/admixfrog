import argparse
from pprint import pprint, pformat
import numpy as np
import logging
from os.path import isfile

from .utils.log import setup_log
from .utils.bam import process_bam, process_bam2
from .utils.vcf import (
    load_pop_file,
    vcf_to_ref,
    load_random_read_samples,
    vcf_to_sample,
)
from .options import add_target_file_options
from .options import add_pop_options
from .options import add_ref_options


def do_ref():
    """subprogram to create reference file"""
    setup_log()
    parser = argparse.ArgumentParser(description="create reference file from vcf")
    parser.add_argument(
        "--outfile", "--out", required=True, help="output file name (xz-zipped)"
    )

    add_pop_options(parser, states_only=True)
    add_ref_options(parser)
    args = parser.parse_args()
    logging.info(pformat(args))

    pop2sample = load_pop_file(args.state_file, args.states)
    random_read_samples = load_random_read_samples(
        args.state_file, args.random_read_samples
    )
    logging.debug(pformat(random_read_samples))
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


def bam():
    """create input file from bam/vcf"""
    setup_log()
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
    parser.add_argument(
        "--max-reads",
        "--max-reads-per-snp",
        type=int,
        default="100",
        help="""Maximum number of reads to retain per SNP. High numbers may lead
            to underflow errors in admixfrog
        """,
    )
    add_target_file_options(parser)
    args = vars(parser.parse_args())

    logging.info(pformat(args))
    force_target_file = args.pop("force_target_file")
    if isfile(args["outfile"]) and not force_target_file:
        raise ValueError(
            """target file exists. set --force-target-file to regenerate"""
        )

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
        #del args["chroms"]
        process_bam(**args)


def bam2():
    """create input file from bam/vcf"""
    setup_log()
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
        "--pseudo_haploid",
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

    logging.info(pformat(args))
    force_target_file = args.pop("force_target_file")
    if isfile(args["outfile"]) and not force_target_file:
        raise ValueError(
            """target file exists. set --force-target-file to regenerate"""
        )

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
        process_bam2(**args)
