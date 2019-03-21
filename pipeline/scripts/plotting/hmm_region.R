source("scripts/plotting/lib.R")


bin_size = as.integer(snakemake@wildcards$bin_size)
infile = snakemake@input$bin
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
cutoff = as.numeric(snakemake@wildcards$cutoff)
region = snakemake@wildcards$region
R = snakemake@config$region[[region]]
chrom_ = as.integer(R$chrom)
start = as.numeric(R$start)
end = as.numeric(R$end)
n_samples = length(names)
p_max = snakemake@params$pmax
p_min = snakemake@params$pmin

if(is.null(start)){
stop("region not found. check config/regions.yaml")
}

data = load_bin_data(infile, names) %>%
    filter(chrom==chrom_, pos < end, pos > start)

TRACK =	get_track(data, snakemake@wildcards$TRACK, p_min, p_max)

d2 = bin_to_long(data) %>% 
	filter( variable %in% TRACK, value>1e-2) 

P2 = bin_colplot_map(d2, add_chrom=F) + 
	facet_wrap(~sample, ncol=1, strip.position="left") +
        ggtitle(sprintf("[%s]%d-%d : %s", chrom_, start, end, region)) 
ggsave(snakemake@output$mapplot, P2, width=10, height=1 * (n_samples + 1), limitsize=F)
P1 = bin_colplot_pos(d2, add_chrom=F) + 
	facet_wrap(~sample, ncol=1, strip.position="left") +
        ggtitle(sprintf("[%s]%d-%d : %s", chrom_, start, end, region)) 

ggsave(snakemake@output$posplot, width=10, height=1 * (n_samples + 1), limitsize=F)
