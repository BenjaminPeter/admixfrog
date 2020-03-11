require(tidyverse)
get_ABBA <- function(p1, p2, p3, p4)
    ABBA = (1 - p1) * p2 * p3 * (1-p4)
get_BABA <- function(p1, p2, p3, p4)
    BABA = p1 * (1 - p2) * p3 * (1-p4)

nABBA <- function(...)sum(get_ABBA(...))
nBABA <- function(...)sum(get_BABA(...))

D = function(ABBA, BABA){
    ABBA = sum(ABBA)
    BABA = sum(BABA)
    return(sum(ABBA - BABA) / (ABBA + BABA ))
}

B = read_csv("tables/paper/bvals.snp.csv.xz", col_types=cols(chrom='c'))

BBIN_SIZE = 50
N_BINS =  2 
P_BIN_SIZE = 5e6
M_BIN_SIZE = 5

data = B %>% mutate(NEA=NEA_alt / (NEA_alt + NEA_ref),
                    DEN=DEN_alt / (DEN_alt + DEN_ref), 
                    PAN=PAN_alt / (PAN_alt + PAN_ref), 
                    AFR=AFR_alt / (AFR_alt + AFR_ref)) %>% 
    dplyr::select(chrom, pos, map, NEA, DEN, PAN, AFR, bval) %>% filter(!is.na(NEA+DEN+PAN+AFR)) %>%
    mutate(bbin = cut(bval, quantile(bval, 0:N_BINS/N_BINS), include.lowest=T)) %>%
    mutate(nbin = as.integer(bbin)) %>%
    mutate(pbin = floor(pos / P_BIN_SIZE)) %>%
    mutate(mbin = floor(map / M_BIN_SIZE)) 


require(bootstrap)
df0 = data %>% 
    group_by(bbin, mbin, chrom) %>% 
    summarize(ABBA = nABBA(NEA, DEN, AFR, PAN), 
              BABA=nBABA(NEA, DEN, AFR, PAN)) 
x =  df0 %>%
    group_by(bbin) %>% 
    summarize(m=D(ABBA, BABA), s=jackknife(1:n(), function(i)D(ABBA[i], BABA[i]))$jack.se)
P = x %>% ggplot(aes(x=bbin, y=m, ymax=m+2*s, ymin=m-2*s)) + geom_point() + 
    geom_errorbar() + 
    geom_hline(yintercept=0) +
    ylab("D(NEA, DEN | AFR, PAN)")

