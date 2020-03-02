library(tidyverse)
source("scripts/paper/long_frag_generic.R")


paper_altai <- function(){
    ref = read_csv(snakemake@input$ref, col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/altai_hcneaden"

    TARGET = 'DEN'

    #runs = read_csv(sprintf("%s.rle.xz", base))
    #snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    runs = read_csv(snakemake@input$rle_altai)
    snp = read_csv(snakemake@input$snp_altai, col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom!='X' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     filter_ambiguous=T, 
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 1,
                     ignore_bins=T,
                     large_top = 2,
                     min_freq = 0.0,
                     one_snp_per_bin=F,
                     ext=c(5e5, 3e5), 
                     p_read_name = 'Denisova 5',
                     pops=c("VIN", "DEN"), base_pop='VIN') +
        theme(axis.text.x = element_blank())
    ggsave("figures/paper/longest/d5_run.png", P, width=7, height=1.75)
}

paper_altai()
