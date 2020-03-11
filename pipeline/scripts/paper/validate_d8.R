library(tidyverse)
source("scripts/paper/long_frag_generic.R")


paper_d8 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/denisova8_hcneaden"

    TARGET = 'NEA'

    #runs = read_csv(sprintf("%s.rle.xz", base))
    #snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    runs = read_csv(snakemake@input$rle_d8)
    snp = read_csv(snakemake@input$snp_d8, col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom == 'X') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_fixed_strict=F,
                     filter_ambiguous=T, 
                     plot_coverage=T,
                     plot_est = F,
                     ignore_bins=T,
                     min_cov = 1,
                     large_top = 2.5,
                     min_freq = 0.1,
                     one_snp_per_bin=F,
                     p_read_name='Denisova 8',
                     ext=c(2e6, 2e6), 
                     pops=c("NEA", "DEN", "AFR"), base_pop='AFR') 
    ggsave("figures/paper/longest/d8_run.png", P, width=7.2, height=2.5)
}


paper_d8()
