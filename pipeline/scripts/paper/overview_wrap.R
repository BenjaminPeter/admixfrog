library(yaml)
library(tidyverse)
source("scripts/plotting/lib.R")
#save.image("bla")

bin_size = as.integer(snakemake@wildcards$bin_size)
infile = snakemake@input$bin
names = snakemake@wildcards$sample
p_max = snakemake@params$pmax
p_min = snakemake@params$pmin
base_size = ifelse(snakemake@params$type == 'paper', 9, 12)
print(c(p_min, p_max))

data = load_bin_data(infile, names)
TRACK = get_track(data, snakemake@wildcards$TRACK, p_min, p_max)
print(TRACK)

d2 = bin_to_long(data) %>% 
   #filter(value > .2)  %>% #' debug track'
    filter( variable %in% TRACK) 

P2 = bin_colplot_wrap(d2, base_size=base_size) 
if(snakemake@params$type == 'paper'){
    P2 =  P2 + 
        theme(
              legend.title=element_blank(),
              legend.position='none',
              panel.spacing.y = unit(.1, "lines")
        )
    ggsave(snakemake@output$mapplot, P2, width=3.5, height=2.0, dpi=1000)

} else {
    P2 = P2 +
        theme(
              legend.title=element_blank(),
              legend.position='bottom'
        )  
    ggsave(snakemake@output$mapplot, P2, width=3.6*2, height=2.4*2)
}

