source("scripts/plotting/meancol.R")
library(tidyverse)
#files = list.files("admixfrog/2000/AFR_VIN_DEN/", "archaicadmixture.res.xz", full=T)
files = snakemake@input$rle
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
bin_size =  as.integer(snakemake@wildcards$bin_size) / 1e6
n_reps=  as.integer(snakemake@params$n_reps)
states = strsplit(snakemake@wildcards$states, "_")[[1]]


 b = lapply(files, read_csv, col_types='ddddd') %>% 
     bind_rows(.id="sample")  %>% 
     mutate(sample= names[as.integer(.$sample)])

P = b %>% 
    mutate(len=len * bin_size, n=n/n_reps / bin_size) %>% 
    group_by(sample, len, state) %>% 
    summarize(n=sum(n)) %>% 
    filter(state>0) %>% 
    mutate(state = states[state+1]) %>%
    ggplot(aes(x=len, y=n, color=as.factor(state))) + geom_step() + 
    facet_wrap(~sample, scale="free_y") + 
    coord_cartesian(xlim=c(0,2.5)) + 
    ylim(0, 5 / bin_size) + 
    col_scale()

ggsave(snakemake@output$rle, P, width=20, height=11)
P2 = P + scale_y_log10()


ggsave(snakemake@output$logrle, P2, width=20, height=11)
