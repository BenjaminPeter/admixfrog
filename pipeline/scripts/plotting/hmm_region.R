source("scripts/comparison_plot.R")

save.image(".rdebug")

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

if(is.null(start)){
stop("region not found. check config/regions.yaml")
}

p_max = snakemake@params$pmax
p_min = snakemake@params$pmin

if(!is.null(snakemake@wildcards$TRACK )){
    TRACK = strsplit(snakemake@wildcards$TRACK, "_")[[1]]
} else{
    TRACK = NULL
}

data = load_data(infile, names) %>%
    filter(chrom==chrom_, pos < end, pos > start)
if(is.null(TRACK)){
    v <- data %>% 
	dplyr::select(-1:-9)  %>% 
	summarize_all(.funs=mean)  %>% 
	unlist 
    TRACK <- names(v[v<=p_max& v>=p_min])                                                           
    TRACK = TRACK[TRACK != "GARBAGE"]
    print(v)
} 
print(TRACK)
print(c(p_max, p_min))


d2 = get_long_data(data) %>% filter( variable %in% TRACK) 

P =  ggplot() + 
	#geom_line(data=d2, aes(x=bin_pos/1e6, y=value, color=variable, fill=variable), lwd=.3) + 
	geom_col(data=d2, mapping=aes(x=pos/1e6, y=value, fill=variable), width=bin_size/1e6) + 
#	geom_point(data=snps, mapping=aes(x=pos/1e6, y=as.numeric(deam)), pch=".", size=.1) +
        facet_wrap(~sample, ncol=1, strip.position="left") +
        ggtitle(sprintf("[%s]%d-%d : %s", chrom_, start, end, region))

ggsave(snakemake@output$plot, width=10, height=1 * (n_samples + 1), limitsize=F)
