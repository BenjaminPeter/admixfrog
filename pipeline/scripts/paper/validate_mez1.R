library(tidyverse)
source("scripts/paper/long_frag_generic.R")


paper_d8 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="error2CAFR/5000/NEA_DEN/mez1_hcneaden"

    TARGET = 'DEN'

    runs = read_csv(sprintf("rle/%s.rle0.4.xz", base), col_types=cols(chrom='c'))
    snp = read_csv(sprintf("admixfrog/%s.snp.xz", base), col_types=cols(chrom=col_character()))
    #runs = read_csv(snakemake@input$rle_d8)
    #snp = read_csv(snakemake@input$snp_d8, col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_fixed_strict=F,
                     filter_ambiguous=F, 
                     plot_coverage=T,
                     plot_est = F,
                     ignore_bins=T,
                     min_cov = 1,
                     large_top = 2.5,
                     min_freq = 0.1,
                     one_snp_per_bin=T,
                     p_read_name='Mez1',
                     ext=c(1e4, 1e4), 
                     pops=c("ALT", "VIN", "DEN"), base_pop='DEN') 
    ggsave("figures/paper/longest/mez1_run.png", P, width=7.2, height=2.5)
}

paper_d8 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="error2CAFR/5000/NEA_DEN/mez2_hcneaden"

    TARGET = 'DEN'

    runs = read_csv(sprintf("rle/%s.rle0.4.xz", base), col_types=cols(chrom='c'))
    snp = read_csv(sprintf("admixfrog/%s.snp.xz", base), col_types=cols(chrom=col_character()))
    #runs = read_csv(snakemake@input$rle_d8)
    #snp = read_csv(snakemake@input$snp_d8, col_types=cols(chrom=col_character()))
    run = runs %>% filter(chrom==7, target==TARGET, type=='state') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_fixed_strict=F,
                     filter_ambiguous=F, 
                     plot_coverage=T,
                     plot_est = F,
                     ignore_bins=T,
                     min_cov = 1,
                     large_top = 2.5,
                     min_freq = 0.1,
                     one_snp_per_bin=F,
                     p_read_name='Mez1',
                     ext=c(1e5, 1e5), 
                     pops=c("NEA", "DEN"), base_pop='DEN') 
    ggsave("figures/paper/longest/mez1_run.png", P, width=7.2, height=2.5)
}


paper_d8()
