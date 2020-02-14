library(tidyverse)
library(scales)
source("scripts/plotting/lib.R")


bin_size = as.integer(snakemake@wildcards$bin_size)
infile = snakemake@input$bin
snpfile = snakemake@input$snp
names = snakemake@wildcards$sample
p_max = snakemake@params$pmax
p_min = snakemake@params$pmin
save.image("rdebug")


data = load_bin_data(infile, names)
TRACK = get_track(data, snakemake@wildcards$TRACK, p_min, p_max)

d2 = bin_to_long(data) %>% filter( variable %in% TRACK)  %>%
    filter(value > 1e-3)


P1 = bin_colplot_pos(d2) + 
        facet_wrap(~chrom, ncol=1, strip.position="left")

ggsave(snakemake@output$posplot, P1, width=20, height=11)
P2 = bin_colplot_map(d2) + 
        facet_wrap(~chrom, ncol=1, strip.position="left")
ggsave(snakemake@output$mapplot, P2, width=20, height=11)

