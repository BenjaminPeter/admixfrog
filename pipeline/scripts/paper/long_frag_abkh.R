library(tidyverse)
source("scripts/plotting/lib.R")
source("scripts/paper/long_frag_generic.R")
library(yaml)

#' run is a data frame with chrom, pos and pos_end
#' ext: how many bp to add left/right
#'




test_abkh1 <- function(){
    ref = read_csv("ref/ref_archaicadmixture.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/gtmode/5000/AFR_NEA_DEN/S_Abkhasian-1_archaicadmixture"

    TARGET = 'DEN'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character())) %>%
        mutate(p=talt/(talt+tref))
    run = runs %>% filter(target==TARGET, type=='state', chrom!='6' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_ambiguous=F, 
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 1,
                     ignore_bins=F,
                     large_top = 1.5,
                     min_freq = 0.0,
                     one_snp_per_bin=F,
                     ext=c(3e4, 5e4), 
                     p_read_name = 'S_Abkhasian-1',
                     pops=c("AFR", "NEA", "DEN"), base_pop='AFR')
    ggsave("figures/paper/longest/abkh_run.png", P, width=8, height=2.2)
}
