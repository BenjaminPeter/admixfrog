require(tidyverse)
require(scales)
a = read_csv("stats/frags2panel/gtmode_error/state_NEA_0.1/AFR_NEA_DEN/10000/simons_paper5b_archaicadmixture.frags")
bins = read_csv("stats/frags2panel/gtmode_error/state_NEA_0.1/AFR_NEA_DEN/10000/simons_paper5b_archaicadmixture.frags")

data0 = a %>% 
    filter(map_len > 0.1) 

data = data0 %>%
    rowwise %>%
    do(idx=seq(.$start, .$end)) %>%
    bind_cols(data0, .) %>%
    unnest

D = data %>% 
    group_by(chrom, idx, age) %>% 
    tally %>%
    arrange(-age)

P = D %>% 
    ggplot(aes(x=idx, y=n, color=age, fill=age))+ 
    geom_col(position='stack') + 
    scale_color_viridis_c(limits = c(0, 30000), oob=squish, aesthetic='both') +
    facet_wrap(~chrom, ncol=1, strip='l', scale='free_x')

D2 = D %>% mutate(present=age==0) %>% 
    group_by(chrom, idx, present) %>% 
    summarize(n=sum(n))
cuts = c(-1, 1, 15000, 25000, 50000)
D3 = D %>% mutate(age=cut(age, cuts)) %>% 
    group_by(chrom, idx, age) %>% 
    summarize(n=sum(n)) %>%
    ungroup


ss = data %>% select(sample, age) %>% distinct %>% mutate(age=cut(age, cuts))  %>% 
    group_by(age) %>% tally %>%mutate(age = as.factor(as.numeric(age))) %>%
    rename(tot=n)

D4 = D3 %>% mutate(age = as.factor(as.numeric(age))) %>% filter(age %in% c(1,4)) %>%
    left_join(ss) %>%
    mutate(n=ifelse(age==1, -n, n))

P = D4 %>% filter(chrom==9) %>% ggplot(aes(x=idx, y=n/tot, color=age, fill=age)) + geom_col(position='stack')
ggsave("figures/chrom9.png", P, widt=11, height=7)


ids = D4 %>% filter(chrom==9 ) %>% arrange(n) %>% filter(n < - 100)
snps  = read_csv("admixfrog/error/10000/AFR_NEA_DEN/Yana1_archaicadmixture.snp.xz")
V = ids %>% left_join(snps, by=c('chrom', idx='bin'))
V %>% select(chrom, pos) %>% as.data.frame %>% arrange(-pos) %>% write_csv("chr9_high_freq_snps.csv")
