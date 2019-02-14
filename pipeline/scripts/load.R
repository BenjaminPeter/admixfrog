library(admixr)
library(glue)
library(reshape2)
library(tidyverse)
library(vcfR)
load_bed <- function(bedfile){
    bed <- read_csv(bedfile, 
                    col_types=cols(chrom=col_character(),
                                   pos0=col_skip(),
                                   pos=col_integer(),
                                   ref=col_character(),
                                   alt=col_character(),
                                   .default = col_double()))
}

flip_geno <- function(genofile, bed){
        anames = names(genofile)[-1:-6]
        bs <- left_join(bed, genofile, by=c("chrom", "pos")) %>%
                dplyr::filter(!is.na(ref.y)) %>% 
                select(-id, -gen)
        bs_flip = bs$ref.x == bs$alt.y & bs$ref.y == bs$alt.x
        bs_same = bs$ref.x == bs$ref.y & bs$alt.y == bs$alt.x
        genofile <- bs[,anames]
        snp <- select(bs, chrom, pos, ref=ref.x, alt=alt.x)
        genofile[bs_flip,] = 1 - genofile[bs_flip,]
        genofile <- genofile[bs_flip | bs_same,]
        snp <- snp[bs_flip | bs_same,]
        genofile <- bind_cols(snp, genofile)
}

load_phase_geno <- function(fname, bed){
    e <- eigenstrat(fname)
    snp <- read_snp(e)
    
    geno <- read_geno(e); geno[geno==9] <- NA
    geno <- 1 - geno / 2
    geno <- as.tibble(data.frame(snp, geno)      )   
        load_geno() %>%
        flip_geno(bed)
}

load_vcf_ad <- function(fname){
    v = read.vcfR(fname)
    G = extract_gt_tidy(v, "GT", alleles=F, gt_column_prepend="") %>%
	separate(GT, sep="/", into=c("A1", "A2"), convert=T) %>%
	mutate(n_alleles = 2 - is.na(A1) - is.na(A2),
	       n_alt = replace_na(A1 == 1,0) + replace_na(A2 == 1, 0),
	       n_ref = replace_na(A1 == 0,0) + replace_na(A2 == 0, 0)) %>%
	select(Key, Indiv, n_ref, n_alt)
    G = v@fix %>% as_tibble %>% 
	select(CHROM, POS, REF, ALT) %>% 
	mutate(Key=1:nrow(v@fix), POS=as.integer(POS)) %>% 
	filter(str_length(ALT) == 1) %>%
	left_join(G) %>%
	select(-Key)
    names(G) <- tolower(names(G))
    G
}
