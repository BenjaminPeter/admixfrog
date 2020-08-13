library(tidyverse)
source("scripts/plotting/lib.R")



infile = snakemake@input$rle
pmin = snakemake@params$pmin
pmax = snakemake@params$pmax
lmin = snakemake@params$lmin
lmax = snakemake@params$lmax
outfile = snakemake@output$mapplot

save.image("rdebug")


R = read_rle(infile) %>% filter(type=='state')

v = R %>% group_by(target, type) %>% 
    summarize(l=sum(map_len)) %>% 
    filter(type=='state') %>% 
    ungroup %>% 
    mutate(l=l/sum(l))

states = v %>% filter(l>pmin, l<pmax) %>% select(target) %>% unlist

P= R %>% filter(target %in% states) %>%
	rle_plot_map(minlen=lmin, maxlen=lmax)

ggsave(outfile, P, width=20, height=11.5)
