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
from .utils.geno_io import read_geno
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

    parser.add_argument(
        "--geno",
        help="""geno format input file"""
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
    elif "geno" in args and args["geno"] is not None:
        geno_to_sample(
            outfile=args["outfile"],
            geno=args["geno"],
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

def geno_to_sample(
        outfile, geno, ref_file, random_read, sample_id):
    raise NotImplemented("Not finished")
    df = read_geno(geno, target_ind = sample_id)


    ref = pd.read_csv(ref_file)
    ref.chrom = ref.chrom.astype(str)

    with lzma.open(outfile, "wt") as infile:
        infile.write(f"chrom,pos,tref,talt\n")

        for chrom in chroms:

            ref_local = ref[ref.chrom == chrom]
            logging.info("sites on chrom %s: %s", chrom, len(ref_local))
            with VariantFile(vcf_file.format(CHROM=chrom)) as vcf:
                vcf.subset_samples([sample_id])
                for i, row in enumerate(vcf.fetch(chrom)):
                    ALT_INDEX = 1
                    if i % 1000 == 0 and i > 0:
                        logging.info(f"{i} i sites processed")
                    if row.pos in ref_local.pos.values:
                        if len(row.alleles) == 1:
                            ALT_INDEX = -1
                        elif len(row.alleles) != 2:
                            ref_row = ref_local[ref_local.pos == row.pos]
                            ref_alt = ref_row.alt.values
                            ALT_INDEX = np.where(ref_alt == row.alleles)[0].item()
                    else:
                        continue

                    snp_str = f"{row.chrom},{row.pos}"

                    sample_data = row.samples
                    for s in sample_data:
                        if random_read:
                            allele = sample_data[s]["GT"][0]
                            if allele == 0:
                                allele_str = "1,0"
                            elif allele == ALT_INDEX:
                                allele_str = "0,1"
                            else:
                                allele_str = "0,0"
                        else:
                            alleles = Counter(sample_data[s]["GT"])
                            allele_str = f"{alleles[0]},{alleles[ALT_INDEX]}"

                    if allele_str != "0,0":
                        infile.write(f"{snp_str},{allele_str}\n")
            logging.info(f"done processing {chrom}")


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
        process_bam2(**args)
