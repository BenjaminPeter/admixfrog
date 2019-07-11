# Admixfrog
Admixfrog is a HMM to infer ancestry frogments (fragments) from low-coverage, contaminated data. 

Briefly, we try to fit the allele frequency at each genomic position in a _target_ by
comparing it with a number of _sources_. In the motivating example, the target would be a
modern human, and the sources would be modern humans (AFR), Neandertals (NEA) or
Denisovans (DEN).

We fit a hidden Markov Model across the genome, with the hidden states being all possible
combinations of ancestry between one or two sources.

## Installation
Requires `python3.6+`
Install dependencies:
```
pip install cython scipy  --upgrade
```
Install `admixfrog` (from github):
```
pip install git+https://github.com/benjaminpeter/admixfrog
```

Install `admixfrog` (from source directory):
```
pip install .
```

## Data
Admixfrog requires (binary-only) eigenstrat data, vcf and bam-files. Supplementary files are typically in yaml-format. The bam-file is used for the _target_, individual, if genotypes are unknown. If genotypes are known, they can be specified in either the eigenstrat or vcf format. In addition, a set of references are required. These too are specified in the reference.

## Quickstart
To get things started, consider an analysis where we would like to learn to local Human, Neandertal and Denisovan ancestry of the Oase1 specimen:

```admixfrog --gfile data/oase --target Oase1_d --states NEA=Vindija.DG+Altai.DG YRI=Yoruba.DG  Denisova.DG  --cont YRI --out quickstart```
this will do the following:

1. `--gfile data/oase`: read the file data/oase.geno|snp|ind (eigenstrat-format)
2. `--target Oase1_d`: declare that we would use the sample named `Oase1_d` as the target
3. with `--states NEA=Vindija.DG+Altai.DG YRI=Yoruba.DG  Denisova.DG` we declare the three sources: a) combine the Vindija and Altai populations from the file (third column in the `.ind`) file into a population named NEA, b) use the population Yoruba.DG, but rename it to YRI and c) Denisova.DG is the third possible source
4. `--cont YRI` designates YRI as a proxy for the contaminant. If there is no contamination, estimating it can be disabled using the `--c0 0 --dont-est-contamination` flags.
5. `--out quickstart`: a prefix for all output files

## Running the program
For most analyses, it is often useful to generate the reference-file and target-file before running the main analysis. This is because parsing these files is quite time-consuming, and is not needed for replicate analyses. However, this is not required, and the program will perform all steps automatically if required


Thus, we might run these three commands:
```
mkdir res/
admixfrog-ref --out res/ref_example.xz --vcf-ref data/oase.vcf.gz \
    --state-file data/pops.yaml \
    --rec-file data/maps_chr.9 \
    --states AFR NEA=Altai_snpAD.DG \
    --map-id AA_Map deCODE COMBINED_LD \
    --default-map AA_Map \
    --chroms 9                        
```
To create the target file, we might run
```
admixfrog-bam --bam data/oase_chr9.bam --ref ref_example.xz  --out oase_example.in.xz 
```
and finally, the analysis can be run using
```
admixfrog --infile oase_example.in.xz --ref ref_example.xz --out example1 -b 10000 \
    --states AFR NEA --contamination AFR
```

Thee most useful command is `admixfrog --help` that will give an up-to-date summary of all the parameters.



there are a few optional parameters, the most important are
 - `-b` the bin size (in 10^6cM), when using a recombination map, or in bp when running without (using `-P`)
 - `--ancestral`: a taxon in that specifies the ancestral allele (must be in the
   reference file)
 - `--states`: the potential admixture sources. (Must be in the reference)
 - `--contamination`: the source of contamination. (Must also be in the reference)

For other parameters, see below or type `admixfrog --help`

There are also utilities to create the input file (from a bam file ) and the reference file (from a vcf file) 
from standard formats. These can be called using `admixfrog-bam` or
`admixfrog-ref`, respectively. Their arguments are also accepted by the main 
`admixfrog` program. However, as parsing and creating these files takes typically much
longer than running admixfrog, I recommend generating them first.


The input file is optionally generated from a bam-file:
```
admixfrog --bamfile {x}.bam --ref {y}.ref.xz --out {z} -b 10000 --ancestral PAN --states AFR NEA DEN
```
but this takes quite long for high-coverage genomes.

### Creating the Reference File:
The input file for `admixfrog` can be created from an (indexed) vcf-file using the 
`admixfrog-ref` subprogram:
```bash
    admixfrog-ref --vcf x_{CHROM}.vcf.gz --out x.ref.xz  \
        --states AFR VIN=Vindija33.19 DEN=Denisova \
        --pop-file data.yaml \
        --rec-file rec.{CHROM}
```

The options are:
    - `--vcf` : an indexed file in vcf format. Non-biallelic variants are
      skipped,but everything else is used. Hence, filtering should be done on this file. Use the wildcard `{CHROM}` if files are split by chromosome
    - `--out` : the name of the output file
    - `--states` : the names of the states, which will be used as sources of
      admixture, contamination and ancestral alleles. By convention I use
      all-caps, 3-4 letter abbreviations. There are three possibilities:

        1. a population define in the `pop file`
        2. a sample name from the vcf file. This will create a single-sample
           reference with the same name as the sample
        3. a string of the form `NEA=Altai,Vindija33.19`. This will create a 
           reference named NEA from the samples `Altai` and `Vindija33.19`

    - `--pop-file`: A `yaml`-format file that defines
        which  samples are in which population, and which samples are
        (pseudo)-haploid

    - `--rec-file` A file specifying the recombination map. I use the file  from here: [https://www.well.ox.ac.uk/~anjali/AAmap/](https://www.well.ox.ac.uk/~anjali/AAmap/)

#### File Format Specification
The reference file has the following columns:
- `chrom` is the chromosome (or contig) id
- `pos` is the physical position of this chromosome
- `ref`, `alt` are the two alleles present at this locus
- `map`,  is the genetic position (in cM)
- a number of pairs of `{ID}_alt`, `{ID}_ref` that give the number of 
  non-reference and reference alleles observed for reference `{ID}`,
respectively.

```
    chrom,pos,ref,alt,map,AFK_alt,AFR_alt,ALT_alt,CHA_alt,DEN_alt,EAS_alt,EUR_alt,NEA_alt,PAN_alt,UST_alt,VIN_alt,AFK_ref,AFR_ref,ALT_ref,CHA_ref,DEN_ref,EAS_ref,EUR_ref,NEA_ref,PAN_ref,UST_ref,VIN_ref
    1,570094,G,A,0,20,0,0,0,2,0,0,0,2,0,0,394,2,2,2,0,10,36,6,0,2,2
    1,714019,A,G,0,168,11,2,2,2,0,0,6,2,0,2,246,27,0,0,0,54,118,0,0,2,0
    1,724289,C,A,0,0,0,1,0,0,0,0,1,0,0,0,414,80,1,2,2,94,148,5,2,2,2
    1,724290,A,C,0,0,0,1,0,0,0,0,1,0,0,0,414,80,1,2,2,94,148,5,2,2,2
    1,725389,C,T,0,5,0,2,2,1,0,0,6,2,0,2,409,0,0,0,1,0,0,0,0,2,0
```

#### Population file format
I use `yaml`-formatted files to define populations, as they are an easily
readable data storage format. The format specification is as follows:
The `sampleset`-section defines sources. For example, below we make a source
panel containing the two Neandertals (AltaiNeandertal and Vindija33.19), and a 
source named `EUR` containing three individuals from the SGDP data set. Finally,
I create a panel named `ANC` which contains the aligned chimp (`panTro4`)
sequence.

In addition, I designate two samples (`panTro4` and `Denisova11`) as
pseudo-haploid by listing them under `pseudo_haploid`. For the outgroup
`panTro4`, this is because we do not care about within-chimp variation, and for
Denisova 11, because it is a low-coverage genome and we cannot get confident
genotype calls.

                         
```yaml
    sampleset:               
        NEA:                 
            - AltaiNeandertal
            - Vindija33.19   
        EUR:             
            - "B_Crete-1"    
            - "B_Crete-2"    
            - "B_French-3"   
        ANC:                 
            - panTro4        
                         
    pseudo_haploid:     
        - Denisova11
        - panTro4            
```


### Creating the Input File:
The input file for `admixfrog` can be created from a bam-file using the 
`admixfrog-bam` subprogram:

```
    admixfrog-bam --bam {x}.bam --ref {y}.ref.xz --deam-cutoff 3 --length-bin-size 35  --out {x}.in.xz
```
This will create a file named `{x}.in.xz` in admixfrog input format from
`{x}.bam`. The site will be ascertained on the sites in `{y}.ref.xz`. Reads with
a deamination (C-\>T) in strand direction in the first 3 bases will be considered
separately for purposes of contamination estimations. Reads will also be binned
in bins of size 35bp for contamination estimation.


##### File Format
the infile has 5 mandatory columns, called `chrom`, `pos`, `tref` and `talt`.  `lib` is optional.

The columns are

    - `chrom` is the chromosome (or contig) id
    - `pos` is the physical position of this chromosome
    - `lib` is a library/read group id. Reads are split by `lib` for contamination
      estimates
    - `tref`, `talt` are the number of refernce and non-reference reads observed for
      this position.

```
    chrom,pos,lib,tref,talt
    1,570094,L5733_0_deam,1,0
    1,570094,R9873_0_deam,1,0
    1,714019,R9880_2_nodeam,0,1
    1,724289,L5736_0_nodeam,1,0
    1,724289,L5736_1_nodeam,1,0
    1,724289,L5734_0_nodeam,1,0
    1,724290,L5736_0_nodeam,1,0
    1,724290,L5736_1_nodeam,1,0
    1,724290,L5734_0_nodeam,1,0
```



#### visualization
a simple viz is 
```R
    library(tidyverse)
    a = read_csv("admixfrog/5000/AFR_VIN_DEN/Papuan_archaicadmixture.bin.xz")
    a %>% gather(k, v, -chrom:-n_snps) %>% 
        filter(k!="AFR", v>.1) %>%
        ggplot(aes(x=map, y=v, fill=k)) + geom_col() + 
        facet_wrap(~chrom, ncol=1, strip='l')
```

## Output
There are currently six output files. All of them are compressed with LZMA.
 - `*.cont.xz` : contamination estimates for each read group
 - `*.bin.xz` : posterior decoding for each bin along the genome
 - `*.snp.xz` : posterior genotype likelihoods for each SNP, taking contamination into
   acccount
 - `*.pars.xz` : parameter estimates
 - `*.rle.xz` : called runs of ancesstry
 - `*.res.xz` : simulated runs of ancestry


## Documentation
Full documentation is not yeat available, this is a dump of the help file for now.
Changes are that `admixfrog --help` will give more up-to-date info

For the detailed description of the algorithm, see [docs/admixfrog.pdf](docs/admixfrog.pdf)


## Contact
Benjamin Peter [benjamin_peter@eva.mpg.de](benjamin_peter@eva.mpg.de)

```python
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
            "--gt-mode", "--gt",
            help="""Assume genotypes are known.
            """,
            action="store_true",
            default=False
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
            help="""Prior of reference allele frequencies. If None (default, recommended), this is 
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
            help="downsample coverage to a proportion of reads"
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
            "--ld-weighting", "--ld", 
            default=False,
            action="store_true",
            help="""downweight SNP in the same bins to counter ancient LD. Very experimental, optimization appears to be broken"""
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
            help="""allow  haploid (i.e. inbreed) stretches. Experimental""",
        )

        add_bam_parse_group(parser)
        add_rle_parse_group(parser)

        from . import __version__
        log_.info("running admixfrog version %s", __version__)
        args = parser.parse_args()
        V = vars(args)
        log_.info(pformat(V))
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

        log_.info("admixfrog %s", __version__)
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
