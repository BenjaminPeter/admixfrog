library(scales)
library(yaml)
source("scripts/plotting/lib.R")
#save.image("bla")

bin_size = as.integer(snakemake@wildcards$bin_size)
infile = snakemake@input$bin
names = snakemake@wildcards$sample
p_max = snakemake@params$pmax
p_min = snakemake@params$pmin

data = load_bin_data(infile, names)
TRACK = get_track(data, snakemake@wildcards$TRACK, p_min, p_max)

d2 = bin_to_long(data) %>% 
    filter( variable %in% TRACK)  %>%
    filter(value > 1e-1) 

P2 = bin_colplot_map(d2, add_chrom=T) + 
    facet_wrap(~chrom, ncol=1, strip='l') +
    scale_y_continuous("Posterior Probability", breaks=c(), expand=c(0,0, 0, 0)) + 
    theme(
          strip.text.y = element_text(angle=180, size=7),
          legend.title=element_blank(),
          legend.position='right',
          legend.text = element_text(size=7),
          legend.key.size = unit(7,"pt"),
          panel.spacing.y = unit(1, "pt"),
          panel.background=element_rect(color='white')
        )
ggsave(snakemake@output$mapplot, P2, width=6, height=3.0, dpi=600)

