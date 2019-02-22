source("scripts/comparison_plot.R")
source("scripts/meancol.R")

col_scale = readRDS("scripts/plotting/cols.rds")


bin_size = as.integer(snakemake@wildcards$bin_size)
infile = snakemake@input$bin
snpfile = snakemake@input$snp
names = snakemake@wildcards$sample
p_max = snakemake@params$pmax
p_min = snakemake@params$pmin
save.image("rdebug")

if(!is.null(snakemake@wildcards$TRACK )){
    TRACK = strsplit(snakemake@wildcards$TRACK, "_")[[1]]
} else{
    TRACK = NULL
}
print(TRACK)


data = load_data(infile, names)
if(is.null(TRACK)){
    v <- data %>% 
	select(-1:-9) %>%
	summarize_all(.funs=mean)  %>% 
	unlist 
    TRACK <- names(v[v<p_max& v>p_min])                                                           
} 

d2 = get_long_data(data) %>% filter( variable %in% TRACK)  %>%
    filter(value > 1e-2)
if(F){
snps = read_csv(snpfile) %>% 
    mutate(chrom=factor(chrom, levels=unique(chrom))) %>%
    mutate(deam=!endsWith(lib, "nodeam")) %>% 
    select(chrom, pos, tref, talt, deam) %>% 
    filter(tref+talt>0) %>%
    select(chrom, pos, deam) %>% distinct()
}


P =  ggplot() + 
	#geom_line(data=d2, aes(x=bin_pos/1e6, y=value, color=variable, fill=variable), lwd=.3) + 
	geom_col(data=d2, mapping=aes(x=map, y=value, fill=variable), width=bin_size/1e6) + 
	#geom_point(data=snps, mapping=aes(x=map, y=as.numeric(deam)), pch=".", size=.1) +
        facet_wrap(~chrom, ncol=2, strip.position="left") + col_scale +
        theme_classic()

ggsave(snakemake@output$png, P, width=20, height=11)
