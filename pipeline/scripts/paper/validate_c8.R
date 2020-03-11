library(tidyverse)
source("scripts/paper/long_frag_generic.R")

paper_c8_c6 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="error2CEUR/5000/NEA_DEN/chagyrskaya08_hcneaden"

    TARGET = 'DEN'

    snp = read_csv(sprintf("admixfrog/%s.snp.xz", base), col_types=cols(chrom=col_character()))
    runs = read_csv(sprintf("rle/%s.rle0.4.xz", base))
    run = runs %>% filter(chrom==6, pos>29e6, pos < 34e6, target=='DEN', map_len > .2)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_ambiguous=T, 
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 3,
                     ignore_bins=F,
                     large_top = 2,
                     one_snp_per_bin=F,
                     ext=c(2e5, 2e5), 
                     pops=c("NEA", "DEN"), base_pop='NEA',
                     p_read_name='Chagyrskaya08') +
        theme(axis.text.x = element_blank())
    ggsave("figures/paper/longest/c8_run6.png", P, width=7.2, height=1.75)
}

paper_c8_c11 <- function(){
    ref = read_csv(snakemake@input$ref, col_types=cols(chrom=col_character()))
    base="error2CEUR/5000/NEA_DEN/chagyrskaya08_hcneaden"

    TARGET = 'DEN'

    snp = read_csv(snakemake@input$snp_c8, col_types=cols(chrom=col_character()))
    runs = read_csv(snakemake@input$rle_c8, col_types=cols(chrom='c'))
    #snp = read_csv(sprintf("admixfrog/%s.snp.xz", base), col_types=cols(chrom=col_character()))
    #runs = read_csv(sprintf("rle/%s.rle0.2.xz", base), col_types=cols(chrom='c'))
    run = runs %>% filter(chrom==11, target=='DEN', map_len > .2)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_ambiguous=T, 
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 3,
                     ignore_bins=F,
                     large_top = 2,
                     one_snp_per_bin=F,
                     ext=c(2e5, 2e5), 
                     pops=c("VIN", "ALT", "DEN"), base_pop='VIN',
                     p_read_name='Chag8') +
        theme(axis.text.x = element_blank())
    ggsave("figures/paper/longest/c8_run11.png", P, width=7.2, height=1.75)
}

paper_c8_c11()
