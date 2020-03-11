library(tidyverse)
source("scripts/plotting/lib.R")
#save.image("salkhit.rdebug")

bin = read_csv("admixfrog/error/5000/AFR_NEA_DEN/Salkhit_archaicadmixture.bin.xz")
runs = read_csv("admixfrog/error/5000/AFR_NEA_DEN/Salkhit_archaicadmixture.rle.xz")
snp = read_csv("admixfrog/error/5000/AFR_NEA_DEN/Salkhit_archaicadmixture.snp.xz")
ref = read_csv("ref/ref_archaicadmixture.csv.xz")

ubin = read_csv("admixfrog/error/5000/AFR_NEA_DEN/UstIshim_archaicadmixture.bin.xz")
uruns = read_csv("admixfrog/error/5000/AFR_NEA_DEN/UstIshim_archaicadmixture.rle.xz")
usnp = read_csv("admixfrog/error/5000/AFR_NEA_DEN/UstIshim_archaicadmixture.snp.xz")

#longest_run = runs %>% filter(target=="AFRDEN") %>% arrange(-len) %>% head(1)
longest_run = runs %>% filter(target=="AFRDEN") %>% arrange(-len) %>% head(1) %>% tail(1)
longest_run_snp = snp %>% 
    filter(chrom==longest_run$chrom, pos>=longest_run$pos - 2e5, pos < longest_run$pos_end + 2e5) %>%
    mutate(in_region = pos >= longest_run$pos & pos < longest_run$pos_end) %>%
    left_join(ref %>% select(-map)) 

u_longest_run_snp = usnp %>% 
    filter(chrom==longest_run$chrom, pos>=longest_run$pos - 2e5, pos < longest_run$pos_end + 2e5) %>%
    mutate(in_region = pos >= longest_run$pos & pos < longest_run$pos_end) %>%
    left_join(ref %>% select(-map)) 

u_lrs0 = u_longest_run_snp %>% mutate(pafr= AFR_alt / (AFR_alt + AFR_ref),
                                 pnea= NEA_alt / (NEA_alt + NEA_ref),
                                 #peur= EUR_ref / (EUR_alt + EUR_ref),
                                 pden= DEN_alt / (DEN_alt + DEN_ref)) 

lrs0 = longest_run_snp %>% mutate(pafr= AFR_alt / (AFR_alt + AFR_ref),
                                 pnea= NEA_alt / (NEA_alt + NEA_ref),
                                 #peur= EUR_ref / (EUR_alt + EUR_ref),
                                 pden= DEN_alt / (DEN_alt + DEN_ref)) 
lrs = lrs0 %>%
    gather(pafr:pden, key='pop', value='f') %>% 
    select(chrom, snp_id, pos, map, bin, pop,tref, talt, in_region, f) %>%
    mutate(delta_f = abs(talt / (tref+talt) - f))

fix = lrs0 %>% filter(abs(pafr-pden)==1 | abs(pafr-pnea)==1) %>%
    mutate(snp_id=1:nrow(.)) %>%
    gather(pafr:pden, key='pop', value='f') %>% 
    select(chrom, snp_id, pos, map, bin, pop,tref, talt, in_region, f) %>%
    mutate(delta_f = abs(tref / (tref+talt) - f))

lrsb = lrs %>% group_by(chrom, bin, map, pop, in_region) %>% 
    summarize(f=mean(f), delta_f=mean(delta_f))

P1 = lrs %>% ggplot(aes(x=snp_id, y=delta_f, fill=in_region)) + 
    facet_wrap(~pop, ncol=3) + geom_col(width=1) + coord_flip()
P2 = fix %>% ggplot(aes(x=as.factor(pos), y=delta_f, fill=in_region)) + 
    facet_wrap(~pop, ncol=3) + geom_col(width=1) + coord_flip()
P3 = lrs %>% ggplot(aes(x=snp_id, y=delta_f, fill=in_region)) + 
     geom_smooth(span=.1, method='loess', aes(color=pop, fill=pop))
RP = runs %>% filter(type=="state", target != "AFR") %>% 
    mutate(chrom=factor(chrom, levels=1:22)) %>%
    rle_plot_map(minlen=0.2)

ggsave("figures/paper/salkhit_longest_run.png", P2, width=4, height=8)
ggsave("figures/paper/salkhit_longest_run.pdf", P2, width=4, height=8)
ggsave("figures/paper/salkhit_runs.png", RP, width=8, height=4)
ggsave("figures/paper/salkhit_runs.pdf", RP, width=8, height=4)




#TRIALS
x = lrs0 %>% filter(abs(pafr-pden)==1) %>% 
    mutate(p_=talt/(tref+talt)) %>% 
    select(snp_id, pos, p, p_, pafr, pnea, pden) %>% 
    gather('k', 'v', p:pden) %>% 
    mutate(target=pos>longest_run$pos & pos < longest_run$pos_end,
           pos=scales::comma(pos))
x %>% ggplot() + 
    geom_col(aes(x=pos, y=-.05, fill=target)) + 
    geom_col(aes(x=pos, y=v))  + 
    facet_wrap(~k, ncol=5) + coord_flip()

