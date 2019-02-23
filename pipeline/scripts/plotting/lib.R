library(tidyverse)
source("scripts/plotting/meancol.R")
get_track <- function(bin_data, TRACK, p_min, p_max){
    if(!is.null(TRACK)){
        TRACK = strsplit(TRACK, "_")[[1]]
    } else{
        TRACK = NULL
    }
    v <- bin_data %>% 
    select(-1:-9) %>%
    summarize_all(.funs=mean)  %>% 
    unlist 
    TRACK <- names(v[v<p_max& v>p_min])                                                           
    return(TRACK)
}
FACET_THEME = theme(panel.background = element_rect(color='#eeeeee'),
                    strip.background = element_rect(size = .3), 
                    panel.spacing.y = unit(.2, "lines"))

THEME = theme_classic() + FACET_THEME

read_binout <- function(fname){
    COL = cols(
               chrom=col_factor(),
               pos=col_integer(),
               id=col_integer(),
               chrom_id=col_integer(),
               viterbi=col_factor(),
               n_snps=col_integer()
               )
    #COL$cols = strsplit("cdiiilci","")[[1]]
    a <- read_csv(fname, col_types=COL) %>%
        mutate(chrom = sort_chroms(chrom))
}

read_rle <- function(fname){
    read_csv(fname, col_types=cols(chrom=col_factor())) %>%
        mutate(pos = as.integer(pos), 
               chrom = sort_chroms(chrom),
               pos_end = as.integer(pos_end),
               n_snps = as.integer(n_snps),
               gap = as.integer(gap),
               len = as.integer(len),
               type=as.factor(type),
               target=as.factor(target),
               block = as.integer(block)) 
}

make_chrom_limits <- function(){
    ref = read_csv("~/proj/100a/basic_processing/bed/rec/twomillion.bed", col_types='ciiccddddddd') 
    mins = ref %>% select(-ref, -alt) %>% group_by(chrom) %>% summarize_all(min, na.rm=T)
    max = ref %>% select(-ref, -alt) %>% group_by(chrom) %>% summarize_all(max, na.rm=T)
    inner_join(mins, max, by='chrom', suffix=c("_min", '_max')) %>% write_csv("ref/chrom_limits.csv")


}

bg_chrom <- function(ref=NULL, map="AA_Map"){
    x = read_csv("ref/chrom_limits.csv", col_types=cols(chrom=col_factor())) %>% 
        select('chrom', min=sprintf('%s_min', map), max=sprintf("%s_max", map))
    if(!is.null(ref)) x = filter(x, chrom %in% unique(ref$chrom))

    return(geom_rect(data=x, mapping=aes(xmin=min, xmax=max, ymin=-Inf, ymax=Inf), fill='#efefef'))

}

load_rle_data <- function(rlefiles, name){
    a <- lapply(infiles, read_rle)
    names(a) <- name
    a <- bind_rows(a, .id="sample")
}

sort_chroms <- function(chrom)
    chrom <- factor(as.character(chrom), levels=c(1:22,"X", "Y", "mt"))

load_bin_data <- function(infiles, name, widths=TRUE){
    a <- lapply(infiles, read_binout)
    names(a) <- name
    a <- bind_rows(a, .id="sample")
    if(!widths) return(a)
	a %>% group_by(sample, chrom) %>%
	arrange(sample, chrom, pos, map)  %>%
	mutate( pwidth = diff(c(pos, max(pos)))) %>%
	mutate( bwidth = diff(c(map, max(map)))) %>%
	ungroup 
}

bin_to_long <- function(data){
    b <- data %>% 
        gather(variable, value, -sample:-n_snps, -pwidth, -bwidth)
}

bin_colplot_map <- function(d2){
    d2 %>% ggplot() +
        bg_chrom(d2) +
        geom_col(mapping=aes(x=map, 
                       y=value, 
                       fill=variable,
                       width=bwidth)) + 
        scale_x_continuous(expand=c(0,0), name='Position (cM)') + 
        scale_y_continuous(expand=c(0,0), name='Probability') +
        coord_cartesian(ylim=0:1) +
        col_scale() + THEME
}
bin_colplot_pos <- function(d2){
    d2 %>% ggplot() +
        bg_chrom(d2, map='pos' ) +
        geom_col(mapping=aes(x=pos, 
                       y=value, 
                       fill=variable,
                       width=pwidth)) + 
        scale_x_continuous(expand=c(0,0), name='Position (cM)') + 
        scale_y_continuous(expand=c(0,0), name='Probability') +
        coord_cartesian(ylim=0:1) +
        col_scale() + THEME
}

rle_plot_map <- function(data, minlen=0.1, maxlen=1, bin_size=1){
    data %>%
	mutate(len = len * bin_size) %>%
	ggplot() +
	bg_chrom(data) + 
	geom_tile(aes(x =(map +map_end)/ 2, width = map_end - map, 
		   y=target, height=1, fill=pmin(len,maxlen))) + 
	facet_wrap(~chrom, ncol=1, strip='left') + 
	THEME + 
	scale_fill_distiller(type='seq',limits=c(minlen, maxlen), na.value="#ffffff", palette='OrRd', direction=1, name='Length(cM)') + 
	#scale_fill_viridis_c(limits=c(minlen, maxlen), na.value="#ffffff", option='A', direction=-1, name='Length(cM)') + 
	scale_x_continuous(expand=c(0,0), name='Position (cM)') + 
	scale_y_discrete(expand=c(0,0), name='state')
}

