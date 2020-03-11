library(tidyverse)
source("scripts/paper/long_frag_generic.R")




paper_d2_chr11 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="error2CAFR/5000/NEA_DEN/denisova2_hcneaden"


    TARGET = 'NEA'

    snp = read_csv(sprintf("admixfrog/%s.snp.xz", base), col_types=cols(chrom=col_character()))
    runs = read_csv(sprintf("rle/%s.rle0.5.xz", base))
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

paper_d8 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/denisova8_hcneaden"

    TARGET = 'NEA'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom == 'X') %>% arrange(-len) %>% head(1) %>% tail(1)

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
                     p_read_name='Denisova 8',
                     ext=c(2e6, 2e6), 
                     pops=c("NEA", "DEN", "AFR"), base_pop='AFR') 
    ggsave("figures/paper/longest/d8_run.png", P, width=7.2, height=2.5)
}

paper_altai <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/altai_hcneaden"

    TARGET = 'DEN'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
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

paper_d2_chr11()
paper_d8()
paper_altai()
