# Admixfrog
Admixfrog is a HMM to infer ancestry frogments (fragments) from low-coverage, contaminated data. 

Briefly, we try to fit the allele frequency at each genomic position in a _target_ by
comparing it with a number of _sources_. In the motivating example, the target would be a
modern human, and the sources would be modern humans (AFR), Neandertals (VIN) or
Denisovans (DEN).

We fit a hidden Markov Model across the genome, with the hidden states being all possible
combinations of ancestry between one or two sources.

## Installation
Requires `python3`
Install dependencies:
```
pip install cython numba numpy scipy pandas setuptools --upgrade
```

Install `admixfrog` (from source directory):
```
pip install .
```

For the pipeline to work, also needs `pgdirect`:
[https://vcs.eva.mpg.de/benjamin_peter/pgdirect](https://vcs.eva.mpg.de/benjamin_peter/pgdirect)

## Output
There are four output files:
 - `*.cont.xz` : contamination estimates for each read group
 - `*.bin.xz` : posterior decoding for each bin along the genome
 - `*.snp.xz` : posterior genotype likelihoods for each SNP, taking contamination into
   acccount
- `*.pars.xz` : parameter estimates



## Documentation
The documentation is under construction, this is a dump of the help file for now.

For a description of the algorithm, see [docs/admixfrog.pdf](docs/admixfrog.pdf)


## Contact
Benjamin Peter [benjamin_peter@eva.mpg.de](benjamin_peter@eva.mpg.de)

```python
    def run():
	parser = argparse.ArgumentParser(
	    description="Infer admixture frogments from low-coverage and contaminated genomes"
	    #        formatter_class=argparse.RawTextHelpFormatter
	)

	parser.add_argument(
	    "-v", "--version", action="version", version="%(prog)s " + admixfrog.__version__
	)

	parser.add_argument(
	    "--infile",
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
	    default=1e-5,
	    help="""Prior of reference allele frequencies. This number is added to both the
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
	    default=1,
	    help="initial tau (should be in [0;1]) (default 1), at most 1 per source",
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

	from . import __version__

	logging.info("admixfrog %s", __version__)
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
```
