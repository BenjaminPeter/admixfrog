# Admixfrog
Admixfrog is a HMM to infer ancestry frogments (fragments) from low-coverage, contaminated data. 

Briefly, we try to fit the allele frequency at each genomic position in a _target_ by
comparing it with a number of _sources_. In the motivating example, the target would be a
modern human, and the sources would be modern humans (AFR), Neandertals (NEA) or
Denisovans (DEN).

We fit a hidden Markov Model across the genome, with the hidden states being all possible
combinations of ancestry between one or two sources.

## Installation
Requires `python3.8+`
Install dependencies:
```
pip install cython scipy  --upgrade
```
Install `admixfrog` (from github):
```
pip install git+https://github.com/benjaminpeter/admixfrog@0.7.3
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
 - `*.pars.yaml` : parameter estimates
 - `*.rle.xz` : called runs of ancesstry
 - `*.res.xz` : simulated runs of ancestry

### Contamination estimates (admixfrog.cont.xz)
The contamination and error estimates are in an  xz-compressed csv format and
will look like this:
```
lib,cont,error,rg,len_bin,deam,n_snps           
SR_nodeam,0.356971,0.010000,SR_nodeam,0,NA,467  
SR_deam,0.000554,0.010000,SR_deam,0,NA,76       
```

Eac row represents a subset of reads for which error and contamination rates estimated
independently. The columns are

- `lib` : a unique string used to group reads. This can be any value, but the
  program tries to split the string according to the format `{rg}_{len_bin}_{deam}`.
If present in this way, the corresponding columns will be filled
- `rg` : Read group
- `len_bin` : Length-bin
- `deam` : whether reads have a terminal deamination
- `cont` : contamination estimate
- `error` : sequencing error estimate
- `n_snps `: how many reads are in this class


### Posterior decoding (admixfrog.bin.xz)
The posterior decoding is in xz-compressed csv format and will look like this
```
chrom,map,pos,id,haploid,viterbi,n_snps,AFK,ARC,AFKARC         
9,200000.000000,281845,0,False,AFK,1,0.541009,0.208888,0.250103
9,300000.000000,300000,1,False,AFK,0,0.540493,0.205282,0.254225
9,400000.000000,400000,2,False,AFK,0,0.539910,0.200351,0.259739
```
Each row represents a bin used in the HMM-algorithm, and the columns are

 - `chrom`: chromosome of bin
 - `map` : map (genetic) coordinate of lower bin boundary
 - `pos` : physical coordinate of lower bin boundary
 - `id` : id of bin (unique number, starting from 0, ordered along chromosome)
 - `haploid` : flag set to True if bin is haploid
 - `viterbi` : Viterbi (Maximum-likelihood) decoding of bin state
 - `n_snps` : number of observed SNP present in bin

the remaining columns (`AFK`, `ARC`, `AFKARC` in the example) give the posterior probability
for the bin being in a given state. The number of columns will vary according to the references
used, and their values sum up to 1. In the example, there are two homozygous  states (`AFK`, `ARC`)
and a heterozygous state `AFKARC`, designated by a concatenation of the two
strings.

### Posterior genotype likelihood (admixfrog.snp.xz)
Results by SNP. xz-compressed csv format.

```
chrom,pos,map,tref,talt,G0,G1,G2,p,bin                         
9,281845,281845,1,0,-0.409389,-2.114613,-2.795105,0.005443,0   
9,635998,635998,0,1,-1.744065,-0.723350,-0.713062,0.288156,4   
9,660473,660473,1,0,-0.401219,-2.487784,-3.318218,0.002107,4   
9,1004958,1004958,0,1,-3.356600,-0.726711,-0.530971,0.388274,8 
9,1463080,1463080,1,0,-0.361344,-2.485787,-4.174832,0.001701,12
```

Each row is a SNP

 - `chrom`: chromosome SNP is on
 - `pos`: physical position of SNP
 - `map`: genetic position of SNP
 - `tref`: number of reference reads at SNP
 - `talt`: number of alt reads at SNP
 - `G0,G1,G2`: log10-likelihood of SNP state 0, 1, 2
 - `p`: estimated allele frequency of derived allele
 - `bin`: bin-id this SNP is in

### Posterior samples (admixfrog.res.xz)
Samples of the posterior given the learned parameters and data are given in xz-compressed csv format and will look like this
```
len,start,end,state,it,chrom        
1,0,1,AFK,0,9                       
1,0,1,AFK,0,9                       
16,6,22,AFK,0,9                     
8,26,34,AFK,0,9                     
7,36,43,AFK,0,9                     
5,1,6,ARC,0,9   
25,1,26,ARC,0,9 
2,34,36,ARC,0,9 
1,43,44,ARC,0,9 

```
Each row represents a segment in the same state, and the columns are:
- `len` : Length (in bins) of segment
- `start` :  Start(id) of segment
- `end` :  End(id) of segment
- `state` :  State of segment
- `it` :  iteration / sample number of posterior sample
- `chrom` :  chromosome sampled

For example, the above snipped designates the 0th iteration of chromosome 9, 
the first bin is homozyogus for the `AFK` state, then two segments, one 5 bins
long, one 25 bins long, start in the `ARC` state.

### Estimated introgressed fragments (admixfrog.rle.xz)

Called introgressed tracts. Calls are done in two formats: 
1. `state` refers to calls where tracts are continued regardless whether they
   are homozygous or heterozygous
2. `het` and `homo` designate runs that are strictly heterozygous or homozygous

```
chrom,start,end,score,target,type,map,pos,id,map_end,pos_end,id_end,len,map_len,pos_len,nscore
9,154,156,0.145364,AFKARC,het,15600000.000000,15600000,154,15800000.000000,15800000,156,2,200000.000000,200000,0.072682
9,1216,1404,28.729700,AFK,state,121800000.000000,121800000,1216,140600000.000000,140600000,1404,188,18800000.000000,18800000,0.152818
9,1187,1193,0.223771,AFK,state,118900000.000000,118900000,1187,119500000.000000,119500000,1193,6,600000.000000,600000,0.037295
9,250,919,78.011711,AFK,state,25200000.000000,25200000,250,92100000.000000,92100000,919,669,66900000.000000,66900000,0.116609
```
Each row represents a segment in the same state, and the columns are:
- `chrom` : Chromosome the segment is on
- `score`, `nscore` :  Numerical score giving certainty of fragment,
  unnormalized or normalized by bin size
- `target` :  iteration / sample number of posterior sample
- `map_start` `map_end`, `map_len` :  start, end and length in genetic  map
- `pos_start` `pos_end`, `pos_len` :  start, end and length in physical map
- `start`, `end`, `len` : start, end and length in Bin id
- `type` :  type of segment call (zygosity vs simple state)
- `target` :  target state for the segment


### Other parameters (admixfrog.pars.yaml)
In yaml format

- `gamma_names`: names of states. All other parameters are given in this order
- `F`, `tau`: estimates of drift parameters per homozygous state
- `alpha0, alpha0_hap`: stationary probabilities for diploid and haploid states,
  respectively
- `trans`, `trans_hap`: diploid and haploid tranition probability
- `error` : error estimates
- `cont`: contamination estimates
- `sex` : assumed sex of individual



## Documentation
Full documentation is not yeat available, this is a dump of the help file for now.
Changes are that `admixfrog --help` will give more up-to-date info

For the detailed description of the algorithm, see [docs/admixfrog.pdf](docs/admixfrog.pdf)




## Contact
Benjamin Peter [benjamin_peter@eva.mpg.de](benjamin_peter@eva.mpg.de)

```
usage: admixfrog [-h] [-v] [--target-file TARGET_FILE] [--ref REF_FILES]
                 [--filter-delta FILTER_DELTA] [--filter-pos FILTER_POS]
                 [--filter-map FILTER_MAP] [--male] [--female]
                 [--bamfile BAMFILE] [--force-target-file]
                 [--deam-cutoff DEAM_CUTOFF] [--minmapq MINMAPQ]
                 [--length-bin-size LENGTH_BIN_SIZE] [--vcfgt VCFGT]
                 [--target TARGET] [--geno-file GENO_FILE] [--guess-ploidy]
                 [--dont-est-contamination] [--est-error]
                 [--freq-contamination FREQ_CONTAMINATION] [--est-F]
                 [--est-tau] [--freq-F FREQ_F] [--est-inbreeding]
                 [--F0 [F0 [F0 ...]]] [--tau0 [TAU0 [TAU0 ...]]] [--e0 E0]
                 [--c0 C0] [--gt-mode] [-b BIN_SIZE] [--prior PRIOR] [-P]
                 [--max-iter MAX_ITER] [--ll-tol LL_TOL] [--dont-split-lib]
                 [--autosomes-only] [--downsample DOWNSAMPLE]
                 [--init-guess [INIT_GUESS [INIT_GUESS ...]]]
                 [--vcf-ref VCF_REF] [--rec-file REC_FILE]
                 [--rec-rate REC_RATE] [--pos-id POS_ID] [--map-id MAP_ID]
                 [--chroms CHROMS] [--force-ref] [--run-penalty RUN_PENALTY]
                 [--n-post-replicates N_POST_REPLICATES] [--outname OUTNAME]
                 [--no-rle] [--no-snp] [--no-bin] [--no-cont] [--no-rsim]
                 [--no-pars] [--states [STATES [STATES ...]]]
                 [--state-file STATE_FILE] [--cont-id CONT_ID]
                 [--ancestral ANCESTRAL]

Infer admixture frogments from low-coverage and contaminated genomes

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --target-file TARGET_FILE, --infile TARGET_FILE, --in TARGET_FILE
                        Sample input file (csv). Contains individual specific
                        data, obtained from a bam file. - Fields are chrom,
                        pos, map, lib, tref, talt" - chrom: chromosome - pos :
                        physical position (int) - map : rec position (float) -
                        lib : read group. Any string, same string assumes same
                        contamination - tref : number of reference reads
                        observed - talt: number of alt reads observed
  --ref REF_FILES, --ref-file REF_FILES
                        refernce input file (csv). - Fields are chrom, pos,
                        ref, alt, map, X_alt, X_ref - chrom: chromosome - pos
                        : physical position (int) - ref : refrence allele -
                        alt : alternative allele - map : rec position (float)
                        - X_alt, X_ref : alt/ref alleles from any number of
                        sources / contaminant populations. these are used
                        later in --cont-id and --state-id flags
  --filter-delta FILTER_DELTA
                        only use sites with allele frequency difference bigger
                        than DELTA (default off)
  --filter-pos FILTER_POS
                        greedily prune sites to be at least POS positions
                        apart
  --filter-map FILTER_MAP
                        greedily prune sites to be at least MAP recombination
                        distance apart
  --male                Assumes haploid X chromosome. Default is guess from
                        coverage. currently broken
  --female              Assumes diploid X chromosome. Default is guess from
                        coverage
  --vcfgt VCFGT, --vcf-gt VCFGT, --vcf-target_file VCFGT
                        VCF input file. To generate input format for admixfrog
                        in genotype mode, use this.
  --target TARGET, --sample-id TARGET
                        sample id if target is read from vcf or geno file. No
                        effect for bam-file
  --chroms CHROMS, --chromosome-files CHROMS
                        The chromosomes to be used in vcf-mode.
  --states [STATES [STATES ...]], --state-ids [STATES [STATES ...]]
                        the allowed sources. The target will be made of a mix
                        of all homozygous and heterozygous combinations of
                        states. More than 4 or 5 sources have not been tested
                        and are not recommended. Must be present in the ref
                        file
  --state-file STATE_FILE, --pop-file STATE_FILE
                        Population assignments (yaml format)
  --cont-id CONT_ID, --cont CONT_ID
                        the source of contamination. Must be specified in ref
                        file
  --ancestral ANCESTRAL, -a ANCESTRAL
                        Outgroup population with the ancestral allele. By
                        default, assume ancestral allele is unknown

bam parsing:
  --bamfile BAMFILE, --bam BAMFILE
                        Bam File to process. Choose this or target_file. The
                        resulting input file will be writen in {out}.in.xz, so
                        it doesn't need to be regenerated. If the input file
                        exists, an error is generated unless --force-target-
                        file is set
  --force-target-file, --force-bam, --force-infile
  --deam-cutoff DEAM_CUTOFF
                        reads with deamination in positions < deam-cutoff are
                        considered separately
  --minmapq MINMAPQ     reads with mapq < MINMAPQ are removed
  --length-bin-size LENGTH_BIN_SIZE
                        if set, reads are binned by length for contamination
                        estimation

geno (Eigenstrat/Admixtools/Reich) format
                                  parser options:
  --geno-file GENO_FILE, --gfile GENO_FILE
                        geno file name (without extension, expects
                        .snp/.ind/.geno files). Only reads binary format for
                        now
  --guess-ploidy        guess ploidy of individuals (use if e.g. random read
                        sample inds are present)

options that control estimation of model
                                  parameters:
  --dont-est-contamination
                        Don't estimate contamination (default do)
  --est-error           estimate sequencing error per rg
  --freq-contamination FREQ_CONTAMINATION, --fc FREQ_CONTAMINATION
                        update frequency for contamination/error (default 1)
  --est-F, -f           Estimate F (distance from ref, default False)
  --est-tau, -tau       Estimate tau (population structure in references)
  --freq-F FREQ_F, --f FREQ_F
                        update frequency for F (default 1)
  --est-inbreeding, -I  allow haploid (i.e. inbreed) stretches. Experimental
  --F0 [F0 [F0 ...]]    initial F (should be in [0;1]) (default 0)
  --tau0 [TAU0 [TAU0 ...]]
                        initial log-tau (default 0), at most 1 per source
  --e0 E0, -e E0        initial error rate
  --c0 C0, -c C0        initial contamination rate

options that control the algorithm behavior:
  --gt-mode, --gt       Assume genotypes are known.
  -b BIN_SIZE, --bin-size BIN_SIZE
                        Size of bins. By default, this is given in 1e-8 cM, so
                        that the unit is approximately the same for runs on
                        physical / map positions
  --prior PRIOR, -p PRIOR
                        Prior of reference allele frequencies. If None
                        (default, recommended), this is estimated from the
                        data This number is added to both the ref and alt
                        allele count for each reference, to reflect the
                        uncertainty in allele frequencies from a sample. If
                        references are stationary with size 2N, this is
                        approximately [\sum_i^{2N}(1/i) 2N]^{-1}.
  -P, --pos-mode        Instad of recombination distances, use physical
                        distances for binning
  --max-iter MAX_ITER, -m MAX_ITER
                        maximum number of iterations
  --ll-tol LL_TOL       stop EM when DeltaLL < ll-tol
  --dont-split-lib      estimate one global contamination parameter (default:
                        one per read group)
  --autosomes-only      Only run autosomes
  --downsample DOWNSAMPLE
                        downsample coverage to a proportion of reads
  --init-guess [INIT_GUESS [INIT_GUESS ...]]
                        init transition so that one state is favored. should
                        be a state in --state-ids

creating reference file:
  --vcf-ref VCF_REF, --vcf VCF_REF
                        VCF File to process. Choose this or reffile. The
                        resulting ref file will be writen as {out}.ref.xz, so
                        it doesn't need to be regenerated. If the input file
                        exists, an error is generated unless --force-ref is
                        set
  --rec-file REC_FILE, --rec REC_FILE
                        Recombination rate file. Modelled after
                        https://www.well.ox.ac.uk/~anjali/AAmap/ If file is
                        split by chromosome, use {CHROM} as wildcards where
                        the chromosome id will be included
  --rec-rate REC_RATE   Constant recombination rate (per generation per base-
                        pair)
  --pos-id POS_ID       column name for position (default: Physical_Pos)
  --map-id MAP_ID       column name for genetic map (default: AA_Map)
  --force-ref, --force-vcf

call introgressed fragments:
  --run-penalty RUN_PENALTY
                        penalty for runs. Lower value means runs are called
                        more stringently (default 0.2)
  --n-post-replicates N_POST_REPLICATES
                        Number of replicates that are sampled from posterior.
                        Useful for parameter estimation and bootstrapping

output name and files to be generated:
  By default, all files are generated. However, if any of the --no-\* options
  are used to disable specific files

  --outname OUTNAME, --out OUTNAME, -o OUTNAME
                        Output file path (without extensions)
  --no-rle              Disabble Estimating runs and writeing to file with
                        extension .rle.xz
  --no-snp              Disable writing posterior genotype likelihood to file
                        with extension .snp.xz
  --no-bin              Disable writing posterior states to file with
                        extension .bin.xz
  --no-cont             Disable writing contamination estimates to file with
                        extension .bin.xz
  --no-rsim             Disable writing posterior simulations of runs to file
                        with extension .res.xz
  --no-pars             Disable writing parameters to file with extension
                        .pars.yaml
```

# Admixslug

Admixslug is a genotype likelihood method for contaminated low-coverage data from Neandertals. It
works by computing a conditional site-frequency spectrum. It uses similar input file formats 
as admixfrog and is therefore for now in the same repository. Outputs of admixslug include 
a contamination estimates and f2, f3, f4 statistics. 

Documentation is still under construction, but a typical command would be


## Input files
Like admixfrog, admixslug requires two input files;

  - a reference file with information from high-quality samples, including the contaminant panel
  - a sample file that stores read information for a sample in compact format

```
admixfrog-bam2 --ref ref/ref_bigsteffi.csv.xz --bamfile bams/bigsteffi/Broion.bam  --out samples2/Broion_bigsteffi.in.xz  --length-bin-size 1 
```
The reference file is created exactly the same way as in admixfrog. The bamfile
contains the reads to be analyzed, and the `--out` flag designates where the
input file will be stored. see `admixfrog-bam2 --help` for details.


## Quickstart
The following command runs admixslug on a single sample stored in
`samples2/Brion_bigsteffi.in.xz` using the sites from `ref/ref_bigsteffi.csv.xz` 
and saving the output files in `admixslug/jk10/ALT_VIN_CHA_DEN/Broion_bigsteffi` 

```
admixslug --infile samples2/Broion_bigsteffi.in.xz \
        --ref ref/ref_bigsteffi.csv.xz \
        -o admixslug/jk10/ALT_VIN_CHA_DEN/Broion_bigsteffi  
        --states ALT VIN CHA DEN  
        --cont-id EUR  
        --ancestral PAN  
        --ll-tol 0.01  
        --ptol 0.001   
        --max-iter 100 
        --filter-pos 50 
        --filter-ancestral  
        --len-bin 2000 
        --jk-resamples 10
        --output-jk-sfs
        --output-fstats
```

The remaining arguments are

  - `--states` : The reference samples or populations to condition the SFS on
  - `--cont-id` : The putative contaminant panel
  - `--ancestral` The ancestral state (these three need to be defined in the
    reference file)
  - `--ll-tol, -ptol`: Convergence criteria in terms of log-likelihood and
    changes in parameter values, respectively
  - `--max-iter` : The maximum number of iterations
  - `--filter-pos` : filter position to be at least x bases apart
  - `--filter-ancestral` : only retain sites with ancestral allele info
  - `--len-bin k `: attempt to bin reads into bins with around k sites. Higher
    numbers of k will result in fewer length-bins for contamination estimation,
    and lower numbers will result in many uncertain estimates
  - `--jk-resamples` : the nubmer of jackknife resamples for standard error
    estimation
  - `--output-jk-sfs`: write a SFS file for each JK resample
  - `--output-fstats`: write all F-stats involving target


## Output
The main outputs are 

#### Contamination file
This file, named {out}.cont.xz contains contamination info for different read groups
and different categories (sequences with no deamination, sequences with terminal 
deamination etc.). 

```
rg            cont      n_exact  n_sites  se_cont   l_cont    h_cont
library1_255_0  0.972531  53135    669098   0.006565  0.959664  0.985398
library2_255_0  0.889277  36158    669098   0.009955  0.869765  0.908789
library3_255_0  0.947004  76550    669098   0.009063  0.929241  0.964766
```

#### SFS file
This file, named {out}.sfs.xz contains info on the estimated SFS

```
sex_chrom  VIN_anc  VIN_der  CHA_anc  CHA_der  ALT_anc  ALT_der  DEN_anc  DEN_der  PAN_anc  PAN_der  F         tau       n_snps  n_reads  n_endo        read_ratio  cont_est  psi        se_tau    l_tau     h_tau     se_F      l_F       h_F
autosome   0        2        0        2        0        2        0        2        1        0        0.289165  0.970069  47982   80951    9278.798637   0.617151    0.885378  0.571462   0.024662  0.921731  1.000000  0.074776  0.142604  0.435727
autosome   0        2        0        2        0        2        2        0        1        0        0.232955  0.893832  35125   57736    6577.306616   0.344759    0.886080  0.274166   0.068626  0.759325  1.000000  0.080282  0.075602  0.390308
autosome   2        0        2        0        2        0        0        2        1        0        0.367088  0.010731  42742   69756    7865.868251   0.169577    0.887237  0.189765   0.010181  0.000000  0.030685  0.239400  0.000000  0.836312
autosome   1        1        0        2        1        1        2        0        1        0        0.423531  0.519948  665     1068     128.319719    0.183521    0.879850  0.137579   0.093208  0.337260  0.702636  0.237326  0.000000  0.888690
```

Admixslug outputs another file named {out}.jksfs.xz which contains the same 
information but for each JK resample. 

#### vcf-file
This file, named {out}.vcf contains a vcf file with i) random read samples, ii)
genotype likelihoods and iii) genotype probabilities for all sites with coverage

#### snp-file
This file, named {out}.snp.xz contains similar info as the VCF file, but more
easily readable in R

```
chrom  pos        map         ref  alt  tref  talt  G0          G1          G2          p         random_read  sfs
1      834832     0.000000    G    C    1     1     -2.209128   -1.459541   -0.018131   0.976466  1            0
1      839495     0.000000    G    T    1     0     -1.454317   -0.826969   -0.088351   0.890396  1            1
1      846864     0.000000    G    C    3     0     -0.004820   -2.028945   -2.774279   0.006359  0            2
1      851204     0.000000    G    C    4     1     -1.674939   -0.910296   -0.067567   0.917391  1            1
1      853267     0.000000    G    T    1     0     -2.084439   -1.408942   -0.021013   0.972267  1            0
```

#### pi-file

```
sex_chrom  pop1              pop2              is_between  pi        sd        sterr
autosome   ALT               ALT               False       0.046972  0.002154  0.046409
autosome   ALT               CHA               True        0.079208  0.001168  0.034179
autosome   ALT               DEN               True        0.308555  0.001155  0.033984
autosome   ALT               PAN               True        0.305752  0.001218  0.034899
```
#### f-files
Currently there are sepeate output files for f2, f3 and f4 statistics, names 
{out}.f2.xz, {out}.f3.xz and {out}.f4.xz. 
These files contain the names of the individuals the statistics is calculated for, 
type of data (autosomal vs. sex chromosome), value of the statistics and the uncertainity. 

Admixslug outputs additional files named {out}.f2.jk.xz, {out}.f3.jk.xz and {out}.f4.jk.xz. 
These files contain the names of the individuals the statistics is calculated for, type of data 
(autosomal vs. sex chromosome) and the value of the statistic for each JK resample. 


Full command here:
```
    usage: admixslug [-h] [-v] [--target-file TARGET_FILE] [--ref REF_FILES] [--seed SEED]
    [--sex-chroms SEX_CHROMS] [--bamfile BAMFILE] [--force-target-file] [--deam-cutoff DEAM_CUTOFF]
    [--minmapq MINMAPQ] [--min-length MIN_LENGTH][--length-bin-size LENGTH_BIN_SIZE] [--report-alleles]
    [--vcfgt VCFGT] [--target TARGET] [--dont-est-contamination] [--dont-est-error] [--est-bias] [--dont-est-F]
    [--est-tau] [--F0 [F0 ...]] [--tau0 [TAU0 ...]] [--e0 E0] [--b0 B0] [--c0 C0] [--max-iter MAX_ITER]
    [--ll-tol LL_TOL] [--ptol PTOL] [--dont-split-lib] [--autosomes-only] [--downsample DOWNSAMPLE]
    [--fake-contamination FAKE_CONTAMINATION] [--deam-bin-size DEAM_BIN_SIZE] [--len-bin-size LEN_BIN_SIZE]
    [--jk-resamples JK_RESAMPLES] [--male] [--female] [--chroms CHROMS] [--outname OUTNAME] [--no-snp] [--no-cont]
    [--no-pars] [--no-sfs] [--output-vcf][--output-jk-sfs] [--output-fstats] [--states [STATES ...]]
    [--het-states [HET_STATES ...]] [--homo-states [HOMO_STATES ...]] [--state-file STATE_FILE]
    [--random-read-samples [RANDOM_READ_SAMPLES ...]] [--cont-id CONT_ID][--ancestral ANCESTRAL]
    [--filter-delta FILTER_DELTA] [--filter-pos FILTER_POS] [--filter-map FILTER_MAP]
    [--filter-high-cov FILTER_HIGH_COV] [--filter-ancestral]


Infer sfs and contamination from low-coverage and contaminated genomes

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --target-file TARGET_FILE, --infile TARGET_FILE, --in TARGET_FILE Sample input file (csv).
                        Contains individual specific data, obtained from a bam file.
                        - Fields are chrom, pos, map, lib, tref, talt"
                        - chrom: chromosome - pos : physical position (int)
                        - map : rec position (float)
                        - lib : read group.
                        Any string, same string assumes same contamination
                        - tref : number of reference reads observed - talt: number of
                        alt reads observed
  --ref REF_FILES, --ref-file REF_FILES refernce input file (csv).
                        - Fields are chrom, pos, ref, alt, map, X_alt, X_ref
                        - chrom: chromosome
                        - pos : physical position (int)
                        - ref : refrence allele
                        - alt : alternative allele
                        - map : rec position (float)
                        - X_alt, X_ref : alt/ref alleles
                        from any number of sources / contaminant populations. these are used later in --cont-id and --state-id flags
  --seed SEED           random number generator seed for resampling
  --sex-chroms SEX_CHROMS
                        The chromosomes to be used as sex chromosomes.
                        If not set, - chromsomes starting wth any of [XYZW] are sex chromosomes
  --vcfgt VCFGT, --vcf-gt VCFGT, --vcf-target_file VCFGT
                        VCF input file. To generate input format for admixfrog in genotype mode, use this.
  --target TARGET, --name TARGET, --sample-id TARGET
                        sample name if target is read from vcf or geno file. written in output of f-stats
  --no-sfs              Disable output of sfs
  --output-vcf          Enable output of vcf
  --output-jk-sfs       write a SFS file for each JK resample
  --output-fstats       write all F-stats involving target
  --states [STATES ...], --state-ids [STATES ...]
                        the allowed sources. The target will be made of a mix of all homozygous and heterozygous
                        combinations of states. More than 4 or 5 sources have not been tested and are not recommended.
                        Must be present in the ref file
  --het-states [HET_STATES ...]
                        Exact het states to be given. If missing or empty, will use all possible het states
  --homo-states [HOMO_STATES ...]
                        Which homozygous states to include. If missing or empty, use all homozygous states
  --state-file STATE_FILE, --pop-file STATE_FILE
                        Population assignments (yaml format). Doesn't currently support het/homo states
  --random-read-samples [RANDOM_READ_SAMPLES ...], --pseudo-haploid [RANDOM_READ_SAMPLES ...]
                        Set a sample as a pseudo-haploid random-read sample for the reference.
                        This means when creating a reference, only one allele is taken.
  --cont-id CONT_ID, --cont CONT_ID
                        the source of contamination. Must be specified in ref file
  --ancestral ANCESTRAL, -a ANCESTRAL
                        Outgroup population with the ancestral allele. By default, assume ancestral allele is unknown
  --filter-delta FILTER_DELTA
                        only use sites with allele frequency difference bigger than DELTA (default off)
  --filter-pos FILTER_POS
                        greedily prune sites to be at least POS positions apart
  --filter-map FILTER_MAP
                        greedily prune sites to be at least MAP recombination distance apart
  --filter-high-cov FILTER_HIGH_COV, --filter-highcov FILTER_HIGH_COV
                        remove SNP with highest coverage (default 0.001, i.e. 0.1% of SNP are removed)
  --filter-ancestral    remove sites with no ancestral allele information

    
bam parsing:
  --bamfile BAMFILE, --bam BAMFILE
                        Bam File to process. Choose this or target_file. The resulting input file will
                        be writen in {out}.in.xz, so it doesn't need to be regenerated. If the input file
                        exists, an error is generated unless --force-target-file is set
  --force-target-file, --force-bam, --force-infile
  --deam-cutoff DEAM_CUTOFF
                        reads with deamination in positions < deam-cutoff are considered separately
  --minmapq MINMAPQ     reads with mapq < MINMAPQ are removed
  --min-length MIN_LENGTH
                        reads with length < MIN_LENGTH are removed
  --length-bin-size LENGTH_BIN_SIZE
                        if set, reads are binned by length for contamination estimation
  --report-alleles      whether contamination/error rates should be conditioned on alleles present at locus

options that control estimation of model
                                  parameters:
  --dont-est-contamination
                        Don't estimate contamination (default do)
  --dont-est-error      estimate sequencing error per rg
  --est-bias            estimate reference bias independent from error
  --dont-est-F          Estimate F (distance from ref, default False)
  --est-tau, -tau       Estimate tau (population structure in references)
  --F0 [F0 ...]         initial F (should be in [0;1]) (default 0)
  --tau0 [TAU0 ...]     initial log-tau (default 0), at most 1 per source
  --e0 E0, -e E0        initial error rate
  --b0 B0, -b B0        initial ref bias rate
  --c0 C0, -c C0        initial contamination rate

options that control the algorithm behavior:
  --max-iter MAX_ITER, -m MAX_ITER
                        maximum number of iterations
  --ll-tol LL_TOL       stop EM when DeltaLL < ll-tol
  --ptol PTOL           stop EM when parameters change by less than ptol
  --dont-split-lib      estimate one global contamination parameter (default: one per read group)
  --autosomes-only      Only run autosomes
  --downsample DOWNSAMPLE
                        downsample coverage to a proportion of reads
  --fake-contamination FAKE_CONTAMINATION
                        Adds fake-contamination from the contamination panel
  --deam-bin-size DEAM_BIN_SIZE, --deam-bin DEAM_BIN_SIZE
                        bin size for deamination
  --len-bin-size LEN_BIN_SIZE, --len-bin LEN_BIN_SIZE
                        bin size for read length
  --jk-resamples JK_RESAMPLES, --n-resamples JK_RESAMPLES
                        number of resamples for Jackknife standard error estimation
  --male                Assumes haploid X chromosome. Default is guess from coverage. currently broken
  --female              Assumes diploid X chromosome. Default is guess from coverage
  --chroms CHROMS, --chromosome-files CHROMS
                        The chromosomes to be used in vcf-mode.


output name and files to be generated:
  By default, all files are generated. However, if any of the --no-* options are used to disable specific files

  --outname OUTNAME, --out OUTNAME, -o OUTNAME
                        Output file path (without extensions)
  --no-snp              Disable writing posterior genotype likelihood to file with extension .snp.xz
  --no-cont             Disable writing contamination estimates to file with extension .bin.xz
  --no-pars             Disable writing parameters to file with extension .pars.yaml

```


