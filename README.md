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

```python
    def run():
        parser = argparse.ArgumentParser(
            description="Infer admixture frogments from low-coverage and contaminated genomes"
        )
        parser.add_argument('infile', 
            help="""Sample input file (csv). Contains individual specific data, obtained from
                            a bam file.

                        - Fields are chrom, pos, map, lib, tref, talt"
                            - chrom: chromosome
                            - pos : physical position (int)
                            - map : rec position (float)
                            - lib : read group. Any string, same string assumes same contamination 
                            - tref : number of reference reads observed
                            - talt: number of alt reads observed"""
        )
        parser.add_argument('ref_file', 
            help="""refernce input file (csv). 
                        - Fields are chrom, pos, ref, alt, map, X_alt, X_ref
                            - chrom: chromosome
                            - pos : physical position (int)
                            - ref : refrence allele
                            - alt : alternative allele
                            - map : rec position (float)
                            - X_alt, X_ref : alt/ref alleles from any number of sources / contaminant populations.
                            these are used later in --cont-id and --state-id flags
                            """
        )
        parser.add_argument('--out', "-o",
                            required=True,
            help="""Output file path (without extensions)"""
        )
        parser.add_argument(
            "--state-ids",
            nargs="*",
            default = ["AFR", "VIN", "DEN"],
            help="""the allowed sources. The target will be made of a mix of all homozygous
            and heterozygous combinations of states. More than 4 or 5 sources have not been
            tested and are not recommended. Must be present in the ref file, with a few
            additional ones:
            - REF : always reference allele
            - ALT : always alternate allele
            - UNIF : allele frequencies are drawn from a uniform / Beta(1, 1) distribution
            - SFS : allele frequencies are drawn from a neutral SFS / Beta(0.5, 0.5) distribution
            """
        )
        parser.add_argument(
            "-b", "--bin-size",
            type = float,
            default = 10000,
            help="""Size of bins. By default, this is given in 1e-8 cM, so that the unit is
            approximately the same for runs on physical / map positions"""
        )
        parser.add_argument(
            "-P", "--pos-mode",
            default = False,
            help="""Instad of recombination distances, use physical distances for binning"""
        )
        parser.add_argument(
            "--max-iter", "-m",
            type = int,
            default = 1000,
            help="""maximum number of iterations"""
        )
        parser.add_argument(
            "--ll-tol",
            type = float,
            default = 1e-2,
            help="""stop EM when DeltaLL < ll-tol"""
        )
        parser.add_argument(
            "--prior", "-p",
            type = float,
            default = 0.5,
            help="""Prior of reference allele frequencies. This number is added to both the
            ref and alt allele count for each reference, to reflect the uncertainty in allele
            frequencies from a sample. IMPORTANT: low values (<1e-4) may lead to
            underflow issues in the betabinomial calculations and should be avoided"""
        )
        parser.add_argument(
            "--no-est-contamination", 
            action="store_false",
            dest="est_contamination",
            default = True,
            help="""Don't estimate contamination (default do)"""
        )
        parser.add_argument(
            "--est-tau", 
            action="store_true",
            default = False,
            help="""Estimate tau (distance from ref, default False)"""
        )
        parser.add_argument(
            "--cont-id",
            default = "AFR",
            help="""the source of contamination. Must be specified in ref file"""
        )
        parser.add_argument(
            "--dont-split-lib",
            action = 'store_false',
            dest = 'split_lib',
            default = True,
            help="""estimate one global contamination parameter (default: one per read
            group)"""
        )

        parser.add_argument(
            "--male",
            dest="sex",
            action="store_const",
            const='m',
            default=None,
            help="Assumes haploid X chromosome. Default is guess from coverage"
        )
        parser.add_argument(
            "--female",
            dest="sex",
            action="store_const",
            const='f',
            default=None,
            help="Assumes diploid X chromosome. Default is guess from coverage"
        )

        args = parser.parse_args()
        pprint(vars(args))
        V = vars(args)
        out = V['out']
        del V['out']

        bins, snps, cont, pars, ll= run_hmm_bb(**vars(args))
        bins.to_csv("%s.bin.xz" % out, float_format="%.6f", index=False, compression="xz")
        snps.to_csv("%s.snp.xz" % out, float_format="%.6f", index=False, compression="xz")
        cont.to_csv("%s.cont.xz" % out, float_format="%.6f", index=False, compression="xz")
        pars.to_csv("%s.pars.xz" % out, float_format="%.6f", index=False, compression="xz")

```

## Contact
Benjamin Peter [benjamin_peter@eva.mpg.de](benjamin_peter@eva.mpg.de)
