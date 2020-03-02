library(tidyverse)
source("scripts/paper/long_frag_generic.R")




paper_d2_chr11 <- function(){
    ref = read_csv(snakemake@input$ref, col_types=cols(chrom=col_character()))
    base="error2CAFR/5000/NEA_DEN/denisova2_hcneaden"


    TARGET = 'NEA'

    #snp = read_csv(sprintf("admixfrog/%s.snp.xz", base), col_types=cols(chrom=col_character()))
    #runs = read_csv(sprintf("rle/%s.rle0.5.xz", base))
    snp = read_csv(snakemake@input$snp_d2, col_types=cols(chrom=col_character()))
    runs = read_csv(snakemake@input$rle_d2)
    run = runs %>% filter(target==TARGET, type=='state' ) %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     filter_fixed_strict=T,
                     plot_coverage=T,
                     plot_est = T,
                     ignore_bins=T,
                     min_cov = 1,
                     large_top = 2.5,
                     min_freq = 0.1,
                     one_snp_per_bin=F,
                     p_read_name='Denisova 2',
                     ext=c(2e6, 2e6), 
                     filter_ambiguous=T, 
                     pops=c("NEA", "DEN", "AFR"), base_pop='AFR')
    ggsave("figures/paper/longest/d2_run11.png", P, width=7.2, height=2.5)
}

paper_d2_chr11()
