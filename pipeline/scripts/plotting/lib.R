library(tidyverse)
source("scripts/plotting/meancol.R")
get_track <- function(data, TRACK, p_min, p_max){
    if(!is.null(TRACK)){
        TRACK = strsplit(TRACK, "_")[[1]]
    } else{
        TRACK = NULL
    }
    v <- data %>% 
    select(-1:-9) %>%
    summarize_all(.funs=mean)  %>% 
    unlist 
    TRACK <- names(v[v<p_max& v>p_min])                                                           
    return(TRACK)
}
FACET_THEME = theme(strip.background = element_rect(size = .3), panel.spacing = unit(.4, "lines"))

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
    a <- read_csv(fname, col_types=COL)
}

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
    d2 %>% 
    ggplot(mapping=aes(x=map, 
                       y=value, 
                       fill=variable,
                       width=bwidth)) + 
        xlab("Position (cM)") + ylab("Probability") +
        ylim(0,1) +
        geom_col() + col_scale() + THEME
}
bin_colplot_pos <- function(d2){
    d2 %>% mutate(pos = pos / 1e6, 
                  pwidth = pwidth /1e6) %>%
    ggplot(mapping=aes(x=(pos+pwidth / 2),
                       y=value, 
                       fill =variable,
                       width=pwidth)) + 
        xlab("Position (Mb)") + ylab("Probability") +
        ylim(0,1) +
        geom_col() + col_scale() + THEME
}
