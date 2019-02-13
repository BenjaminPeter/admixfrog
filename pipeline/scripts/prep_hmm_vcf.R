suppressPackageStartupMessages({
library(admixr)
library(glue)
library(reshape2)
library(tidyverse)
source("scripts/load.R")
})

VIN_ = snakemake@params$VIN
DEN_ = snakemake@params$DEN
CHA_ = snakemake@params$CHA
ALT_ = snakemake@params$ALT
EUR_ = snakemake@config$sampleset$sgdpeur
PAN_ = snakemake@params$PAN
AFR_ = snakemake@config$sampleset$sgdpafr
pops = list(VIN=VIN_, 
	    DEN=DEN_, 
	    CHA=CHA_, 
	    ALT=ALT_, 
	    EUR=EUR_, 
	    PAN=PAN_, 
	    AFR=AFR_)
map_ = snakemake@params$recmap


bed <- load_bed(snakemake@input$bed) %>%
    select(chrom, pos, ref, alt, map=snakemake@params$recmap)
vcfs <- lapply(snakemake@input$vcfs, load_vcf_ad) %>% bind_rows

V = lapply(pops, function(x) vcfs %>% 
	   filter(indiv %in% x) %>% 
	   group_by(chrom, pos, ref, alt) %>% 
	   summarize(n_ref=sum(n_ref), n_alt=sum(n_alt))) %>%
	   bind_rows(.id="panel")

alphas = dcast(V, chrom+pos ~ panel, value.var="n_alt") %>% as_tibble
ref_alleles = dcast(V, chrom+pos ~ panel, value.var="n_ref") %>% as_tibble
    
bed %>% left_join(alphas) %>% 
    inner_join(ref_alleles, by=c("chrom", "pos"), suffix=c("_alt","_ref")) %>%
    write_csv(snakemake@output$ref)
