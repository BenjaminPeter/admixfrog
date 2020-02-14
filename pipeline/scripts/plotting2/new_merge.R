library(tidyverse)
library(scales)

STATES = "AFR_NEA_DEN"
STATE = "NEA"
FRAG_CUTOFF = 0.1
BIN_SIZE = 5000
PANEL1='simons'
PANEL2='paper5c'
MODE1='gtmode'
MODE2='error'




a = read_csv("stats/frags2panel/gtmode_error/state_NEA_0.1/AFR_NEA_DEN/5000/simons_paper5c_archaicadmixture.frags")
meta = read_csv("config/sgdp_ancient.csv")

V = a %>% select(sample, chrom, map, map_end, map_len, pos, pos_end, pos_len, cov) %>% 
    left_join(meta) %>% 
    mutate(label = ifelse(age==0, pop, sample)) %>%
    mutate(label2 = ifelse(age==0, region, sample))
Q = V %>% filter(region!="Africa", region!="America", map_len > 100000) %>% 
    group_by(sample, region, label ,age ) %>% 
    summarize(m=mean(map_len)) %>% 
    group_by(label, age, region) %>% 
    summarize(m=mean(m))
Q2 = Q%>% 
    ggplot(aes(x=age, y=1/m, label=label, color=region)) + geom_text()

load_frags <- function(mode1='gtmode',
                       mode2='error',
                       frag_cutoff=0.1,
                       state='NEA',
                       states='AFR_NEA_DEN',
                       bin_size=5000,
                       panel1='simons',
                       panel2='paper5c'){
    fname = sprintf("stats/frags2panel/%s_%s/state_%s_%s/%s/%s/%s_%s_archaicadmixture.frags",
                    mode1, mode2, state, frag_cutoff, states, bin_size, panel1, panel2)
    a = read_csv(fname)
    a$bin_size = bin_size
    a$mode1=mode1
    a$mode2=mode2
    a$state=state
    a$frag_cutoff = frag_cutoff
    a$states = states
    a$panel1 = panel1
    a$panel2 = panel2
    return(a)
}

load_res <- function(sample="Muierii", 
                     bins = c(2000, 5000, 10000), 
                     mode='error',
                     str='admixfrog/%s/%s/AFR_NEA_DEN/%s_archaicadmixture.res.xz'){
    B = lapply(bins, function(b){
        infile = sprintf(str, mode, b, sample)
        print(infile)
        read_csv(infile) %>% mutate(len = len * b / 1e6)
                     }) %>%
        bind_rows(.id='bin_size') %>%
        mutate(bin_size = as.factor(bins[as.integer(bin_size)])) %>%
        rename(map_len = len) %>%
        group_by(map_len, bin_size, it, state) %>%
        tally
    return(B)
}

CFG = yaml::yaml.load_file("config/panels.yaml")
samples_simons = CFG$panels$simons
samples_paper = CFG$panels$paper5c


frags_called <- lapply(c(2000, 5000, 10000), function(b)load_frags(bin_size=b)) %>% 
    bind_rows  %>%
    left_join(meta)

frags_resampled_sgdp = sapply(samples_simons, load_res, mode='gtmode', simplify=F) %>%  
    bind_rows(.id="sample") %>% 
    filter(state==STATE) %>%
    select(-age) %>%
    left_join(meta) %>% 
    mutate(sample=fct_reorder(sample, age)) %>% 
    write_csv("stats/sgdp_frags_resampled.xz")

frags_resampled_ancient = sapply(samples_paper, load_res, mode='error', simplify=F) %>%  
    bind_rows(.id="sample") %>% 
    left_join(x) %>%  
    mutate(sample=fct_reorder(sample, age)) %>%
    write_csv("stats/sgdp_ancient_resampled.xz")

V2 = V %>%mutate(len=map_len, bin_size='5000') %>% filter(len>=.10, age>9000) %>% 
    mutate(sample=fct_reorder(sample, age))

P = ggplot(mapping=aes(x=len, color=bin_size)) + 
    stat_ecdf(data=L2, lty=1) + 
    stat_ecdf(data=V2, lty=2) + 
    scale_x_reverse() + 
    scale_y_log10() + 
    coord_cartesian(ylim=c(0.01, 1), xlim=c(2, .2), expand=c(0,0,0,0)) + 
    facet_wrap(~sample, ncol=2, strip='l', dir='v')

#figure with tracts from samples
