library(tidyverse)
source("scripts/paper/long_frag_generic.R")

paper_d3_old <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="error2CEUR/5000/A=ALT+D12_D=DEN+D11/denisova3_hcneaden"

    TARGET = 'A'

    runs = read_csv(sprintf("rle/%s.rle0.2.xz", base))
    run = runs %>% filter(target==TARGET, map_len>0.2, chrom==6, pos>30e6, pos < 60e6)  %>%
        arrange(pos)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 3,
                     ignore_bins=T,
                     large_top = 2,
                     one_snp_per_bin=F,
                     ext=c(1e6, 1e6), filter_ambiguous=F, 
                     pops=c("ALT", "DEN"), base_pop='DEN',
                     p_read_name='Denisova3') +
        theme(axis.text.x = element_blank())
    ggsave("figures/paper/longest/d3_run6.png", P, width=7.2, height=2)
}

paper_d3 <- function(){
    ref = read_csv(snakemake@input$ref, col_types=cols(chrom=col_character()))

    TARGET = 'ALT'

    snp = read_csv(snakemake@input$snp_d3, col_types=cols(chrom=col_character()))
    runs = read_csv(snakemake@input$rle_d3)
    run = runs %>% filter(target==TARGET, map_len>0.2, chrom==6, pos>30e6, pos < 60e6)  %>%
        arrange(pos)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 3,
                     ignore_bins=T,
                     large_top = 2,
                     one_snp_per_bin=F,
                     ext=c(1e6, 1e6), filter_ambiguous=F, 
                     pops=c("ALT", "DEN"), base_pop='DEN',
                     p_read_name='Denisova3') +
        theme(axis.text.x = element_blank())
    ggsave(snakemake@output$frag_d3, P, width=7.2, height=1.75)
}

paper_d3()
