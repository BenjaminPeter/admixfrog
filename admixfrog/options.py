"""option groups
"""
# set up a list of option groups to keep things organized
#population comparison options
POP_OPTIONS = ["cont_id", 
               "ref_files", 
               "sex", 
               "states", 
               "state_file",
               "random_read_samples"
               "ancestral"]


#generate target_file from vcf, bam or geno
INFILE_OPTIONS = [
    "bamfile",
    "deam_cutoff",
    "minmapq",
    "force_target_file",
    "length_bin_size",
    "report_alleles",
    "tstv",
    "alleles",
    "vcfgt",
    "target",
    "target_file"
]

#generate reference from vcf or geno
REFFILE_OPTIONS = [
    "vcf_ref",
    "pop_file",
    "rec_file",
    "rec_rate",
    "pos_id",
    "map_id",
    'default_map',
    "chroms",
    "force_ref",
]

#algorithm control options
ALGORITHM_OPTIONS = [
    "bin_size",
    "autosomes_only",
    "downsample",
    "fake_contamination",
    "gt_mode",
    "ll_tol",
    "max_iter",
    "n_post_replicates",
    "pos_mode",
    "prior",
    "split_lib",
]

#geno format options
GENO_OPTIONS=[
    'geno_file',
    'guess_ploidy'
]
def add_pop_options(parser):
    parser.add_argument(
        "--states",
        "--state-ids",
        nargs="*",
        default=["AFR", "NEA", "DEN"],
        help="""the allowed sources. The target will be made of a mix of all homozygous
        and heterozygous combinations of states. More than 4 or 5 sources have not been
        tested and are not recommended. Must be present in the ref file
        """,
    )
    parser.add_argument(
        "--state-file", 
        "--pop-file",
        default=None, help="""Population assignments (yaml format)"""
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
        "--random-read-samples",
        "--pseudo-haploid",
        nargs="*",
        default=[],
        help="""Set a sample as a pseudo-haploid random-read sample for the reference. This means when creating a reference,
        only one allele is taken.
        """
    )

def add_output_options(parser):
    g = parser.add_argument_group(
        "output name and files to be generated",
        """By default, all files are generated. However, if any of 
                                  the --no-* options are used to disable specific files
                                  """,
    )
    g.add_argument(
        "--outname", 
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


def add_target_file_options(parser):
    g = parser.add_argument_group("bam parsing")
    g.add_argument(
        "--bamfile",
        "--bam",
        help="""Bam File to process. Choose this or target_file. 
                   The resulting input file will be writen in {out}.in.xz,
                   so it doesn't need to be regenerated.
                   If the input file exists, an error is generated unless
                   --force-target-file is set
                   """,
    )
    g.add_argument("--force-target-file", "--force-bam", '--force-infile', 
                   default=False, action="store_true")
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
    g.add_argument(
        "--report-alleles",
        dest="report_alleles",
        action="store_true",
        default=False,
        help="""whether contamination/error rates should be conditioned on alleles present at locus""",
    )
    parser.add_argument(
        "--vcfgt",
        "--vcf-gt",
        "--vcf-target_file",
        help="""VCF input file. To generate input format for admixfrog in genotype mode, use this.
        """,
    )
    parser.add_argument(
        "--target",
        "--sample-id",
        help="""sample id if target is read from vcf or geno file. No effect for bam-file
        """,
    )


def add_rle_options(parser):
    g = parser.add_argument_group("call introgressed fragments")
    g.add_argument(
        "--run-penalty",
        type=float,
        default=0.2,
        help="""penalty for runs. Lower value means runs are called more
        stringently (default 0.2)""",
    )
    g.add_argument(
        "--n-post-replicates",
        type=int,
        default=100,
        help="""Number of replicates that are sampled from posterior. Useful for
        parameter estimation and bootstrapping
        """
    )


def add_geno_options(parser):
    g = parser.add_argument_group("""geno (Eigenstrat/Admixtools/Reich) format
                                  parser options""")
    g.add_argument(
        "--geno-file", 
        "--gfile",
        help="""geno file name (without extension, expects .snp/.ind/.geno
        files). Only reads binary format for now""",
    )
    g.add_argument(
        "--guess-ploidy",
        action='store_true',
        default=True,
        help="""guess ploidy of individuals (use if e.g. random read sample
        inds are present)
        """
    )


def add_ref_options(parser):
    g = parser.add_argument_group("creating reference file")
    g.add_argument(
        "--vcf-ref",
        "--vcf",
        help="""VCF File to process. Choose this or reffile. 
                   The resulting ref file will be writen as {out}.ref.xz,
                   so it doesn't need to be regenerated.
                   If the input file exists, an error is generated unless
                   --force-ref is set
                   """,
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
        default=["AA_Map", "deCODE", "YRI_LD", "CEU_LD"],
        nargs="*",
        help="""column name for genetic map (default: AA_Map)
        """,
    )
    parser.add_argument(
        "--default-map",
        default="AA_Map",
        help="""default recombination map column"""
    )
    parser.add_argument(
        "--chroms",
        "--chromosome-files",
        default="1-22,X",
        help="""The chromosomes to be used in vcf-mode.
        """,
    )
    g.add_argument("--force-ref", "--force-vcf", default=False, action="store_true")


def add_estimation_options(P):
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


def add_base_options(P):
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
        approximately  [\\sum_i^{2N}(1/i) 2N]^{-1}.
          """
    )
    parser.add_argument(
        "-P",
        "--pos-mode",
        "--posmode",
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
        "--fake-contamination",
        type=float,
        default=0.0,
        help="Adds fake-contamination from the contamination panel",
    )
    parser.add_argument(
        "--init-guess",
        nargs="*",
        help="""init transition so that one state is favored. should be a 
        state in --state-ids """,
    )



