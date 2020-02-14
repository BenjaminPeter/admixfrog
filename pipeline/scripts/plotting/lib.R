suppressPackageStartupMessages({
library(scales)
library(tidyverse)
library(yaml)
source("scripts/plotting/fit.R")
source("scripts/plotting/meancol.R")
})
get_track <- function(bin_data, TRACK, p_min, p_max){
    if(!is.null(TRACK)){
        TRACK = strsplit(TRACK, "_")[[1]]
    } else{
        TRACK = NULL
    }
    v <- bin_data %>% 
    select(-1:-7) %>%
    summarize_all(.funs=mean)  %>% 
    unlist 
    TRACK <- names(v[ (v<p_max& v>p_min) |  str_length(names(v)) >4   ])                                                           
    return(TRACK)
}
FACET_THEME = theme(panel.background = element_rect(color='#eeeeee'),
                    strip.background = element_rect(size = .3), 
                    panel.spacing.y = unit(.2, "lines"))

THEME = theme_classic() + FACET_THEME

read_binout <- function(fname){
    COL = cols(
               chrom=readr::col_factor(),
               pos=col_integer(),
               id=col_integer(),
               viterbi=readr::col_factor(),
               n_snps=col_integer()
               )
    #COL$cols = strsplit("cdiiilci","")[[1]]
    a <- read_csv(fname, col_types=COL) %>%
        mutate(chrom = sort_chroms(chrom))
}
read_snpout <- function(fname){
    COL = cols(
               chrom=col_factor(),
               pos=col_integer(),
               map=col_double(),
               bin=col_integer(),
               tref=col_integer(),
               talt=col_integer(),
               G0=col_double(),
               G0=col_double(),
               G1=col_double()
               )
    #COL$cols = strsplit("cdiiilci","")[[1]]
    a <- read_csv(fname, col_types=COL) %>%
        mutate(chrom = sort_chroms(chrom))
}

read_rle <- function(fname){
    read_csv(fname, col_types=cols(chrom=readr::col_factor())) %>%
        mutate(pos = as.integer(pos), 
               chrom = sort_chroms(chrom),
               pos_end = as.integer(pos_end),
               len = as.integer(len),
               type=as.factor(type),
               target=as.factor(target))
}

make_chrom_limits <- function(){
    ref = read_csv("~/proj/100a/basic_processing/bed/rec/twomillion.bed", col_types='ciiccddddddd') 
    mins = ref %>% select(-ref, -alt) %>% group_by(chrom) %>% summarize_all(min, na.rm=T)
    max = ref %>% select(-ref, -alt) %>% group_by(chrom) %>% summarize_all(max, na.rm=T)
    inner_join(mins, max, by='chrom', suffix=c("_min", '_max')) %>% write_csv("ref/chrom_limits.csv")


}
bg_chrom_wrap <- function(ref=NULL, map="AA_Map"){
    x = read_csv("annotation/chrom_limits.csv", col_types=cols(chrom=readr::col_factor())) %>% 
        select('chrom', min=sprintf('%s_min', map), max=sprintf("%s_max", map))

    x_wrap = x %>% mutate(c_wrap=high_chrom(chrom),
                          min2 = ifelse(as.character(chrom)==as.character(c_wrap), min, 360 - max),
                          max = ifelse(as.character(chrom)==as.character(c_wrap), max, 360 - min)) %>%
            select(-min) %>%
            rename(min=min2)


    return(geom_rect(data=x_wrap, color=NA, mapping=aes(xmin=min, xmax=max, ymin=-Inf, ymax=Inf), fill='#efefef'))
}

bg_chrom <- function(ref=NULL, map="AA_Map", flip=F){
    x = read_csv("ref/chrom_limits.csv", col_types=cols(chrom=col_factor())) %>% 
        select('chrom', min=sprintf('%s_min', map), max=sprintf("%s_max", map))
    if(!is.null(ref)) x = filter(x, chrom %in% unique(ref$chrom))
    if(flip){
	return(geom_rect(data=x, color=NA, mapping=aes(ymin=min, ymax=max, xmin=-Inf, xmax=Inf), fill='#efefef'))
    }

    return(geom_rect(data=x, color=NA, mapping=aes(xmin=min, xmax=max, ymin=-Inf, ymax=Inf), fill='#efefef'))
}

read_run <- function(rfile){
    read_csv(rfile, col_types='ddddd') %>%
	    group_by(state, len) %>%
	    summarize(n=sum(n))
}
read_run2 <- function(rfile){ #new format
    read_csv(rfile, col_types='iiific') %>%
	    group_by(state, len) %>%
        tally
}

load_runs_data <- function(files, name){
    a <- lapply(files, read_run2)
    names(a) <- name
    a <- bind_rows(a, .id="sample")
}
load_cont_data <- function(contfiles, name){
    a <- lapply(contfiles, read_csv, col_types='cdcici')
    names(a) <- name
    a <- bind_rows(a, .id="sample")
}

load_rle_data <- function(rlefiles, name){
    a <- lapply(rlefiles, read_rle)
    names(a) <- name
    a <- bind_rows(a, .id="sample")
}

sort_chroms <- function(chrom)
    #chrom <- factor(as.character(as.integer(chrom)), levels=c(1:22,"X", "Y", "mt"))
    chrom <- factor(as.character(chrom), levels=c(1:22,"X", "Y", "mt"))

load_bin_data <- function(infiles, name, widths=TRUE){
    a <- lapply(infiles, read_binout)
    names(a) <- name
    a <- bind_rows(a, .id="sample")
    a$sample = factor(a$sample, level=name)
    if(!widths) return(a)
	a %>% group_by(sample, chrom) %>%
	arrange(sample, chrom, pos, map)  %>%
	mutate( pwidth = diff(c(pos, max(pos)))) %>%
	mutate( bwidth = diff(c(map, max(map)))) %>%
	ungroup 
}

load_snp_data <- function(infiles, name, widths=TRUE){
    a <- lapply(infiles, read_snpout)
    names(a) <- name
    a <- bind_rows(a, .id="sample")
    if(!widths) return(a)
	a %>% group_by(sample, chrom) %>%
	arrange(sample, chrom, pos, map)  %>%
	ungroup 
}

bin_to_long <- function(data){
    b <- data %>% 
        gather(variable, value, -sample:-n_snps, -pwidth, -bwidth)
}

bin_colplot_wrap <- function(d2, add_chrom=T, base_size=8){
    print(nrow(d2))
    d_wrap = d2 %>% mutate(c_wrap=high_chrom(chrom), 
                           m_wrap = ifelse(as.character(chrom)==as.character(c_wrap), map, 360 - map))
    d_label = data.frame(c_wrap=levels(d_wrap$c_wrap), 
                         lab=c(paste0(23-c(as.integer(levels(d_wrap$c_wrap)[1:11])), " "), NA))
    d_label$xintercept=c(rep(360, 11), 0)
    P = d_wrap %>% ggplot() 
    if(add_chrom) P = P + bg_chrom_wrap(d_wrap)
    P +
        geom_col(mapping=aes(x=m_wrap, 
                       y=value, 
                       fill=variable,
                       width=bwidth * 1.001),
                        lwd=0) + 
        facet_wrap(~c_wrap, ncol=1, strip='l') + 
        geom_text(aes(label=paste0(" ",c_wrap, " ")), x=-0, y=0.5, data=d_label, hjust=1, size=base_size / 3) + 
        geom_text(aes(label=lab), x=Inf, y=0.5, data=d_label, hjust=.5, size=base_size / 3) +
        geom_vline(aes(xintercept=xintercept), lwd=.33, color='black', data=d_label) + 
        geom_vline(xintercept=0, color='black', lwd=0.33) + 
        scale_y_continuous(name="Pr(Local Ancestry)\n", breaks=c(), expand=expand_scale(0,0)) +
        scale_x_continuous(name="Position (cM)", 
			    breaks=seq(0, 360, 40), expand=expand_scale(c(0,0.04),0.00)) +
        coord_cartesian(ylim=0:1, clip='off') +
        col_scale() + theme_classic(base_size) +
        theme(
          strip.background = element_blank(),

          strip.text = element_blank()
        )
}

bin_colplot_map <- function(d2, add_chrom=T){
    P = d2 %>% ggplot() 
    if(add_chrom) P = P + bg_chrom(d2)
    P +
        geom_col(mapping=aes(x=map, 
                       y=value, 
                       fill=variable,
                       #color=variable,
                       width=bwidth)) + 
        scale_x_continuous(expand=c(0,0), name='Position (cM)', labels=comma) + 
        scale_y_continuous(expand=c(0,0), name='Probability') +
        coord_cartesian(ylim=0:1) +
        col_scale() + THEME
}
bin_colplot_pos <- function(d2, add_chrom=T){
    P = d2 %>% group_by(sample, pos, chrom) %>%
        mutate(value=cumsum(value)) %>%
        mutate(v0 = lag(value, default=0)) %>% 
        ggplot()
    if(add_chrom) P = P + bg_chrom(d2, map='pos')
    P + 
        geom_rect(mapping=aes(xmin=pos, xmax=pos+pwidth,
                       ymin=v0, ymax=value,
                       fill=variable,
                       #width=pwidth
                       )) + 
        scale_x_continuous(expand=c(0,0), name='Position (Mbp)', labels=function(x)x/1e6) + 
        scale_y_continuous(expand=c(0,0), name='Probability') +
        coord_cartesian(ylim=0:1, expand=c(0,0,0,0)) +
        col_scale() + THEME
}

rle_plot_map <- function(data, minlen=0.1, maxlen=1){
    data %>%
	filter(map_len > minlen) %>%
	ggplot() +
	bg_chrom(data) + 
	geom_tile(aes(x =(map +map_end)/ 2, width = map_end - map, 
		   y=0.5, height=1, fill=target )) + 
	facet_wrap(~chrom, ncol=1, strip='left') + 
	THEME + 
	col_scale() +
	scale_alpha_continuous(range=c(0.3,1), limits=c(minlen, maxlen), name='Length(cM)') + 
	scale_x_continuous(expand=c(0,0), name='Position (cM)') + 
	scale_y_discrete(expand=c(0,0), name='State') + 
    theme(strip.text.y = element_text(angle=180))
}


plot_m_gamma <- function(R, generation_time){
    SCALING = generation_time * 100 #* 2
    tmax = pmin(100000 / SCALING, max(R$semean, R$slmean, na.rm=T))
    t <- seq(0, tmax, 0.01)
    n_samples = length(unique(R$sample))


    gpred = R %>% 
        rowwise %>% 
        do(l = data.frame(time=t+.$scaled_age, 
			  mig=dgamma(t, shape=.$shape, scale=1/.$scale))) %>%
        bind_cols(R) %>% 
        filter(!is.na(lmean + emean), slmean < tmax * 1.1) %>%
        unnest(l) 

    is_signif = R$sample[R$delta_ll  > qchisq(.95, 1)]

    M = gpred %>%
        ggplot( aes(x=time*SCALING, y=mig, color=sample %in% is_signif)) + 
        geom_line(lty=2) + 
        scale_x_continuous(name="time (y)", expand=expand_scale(0,0)) + 
        coord_cartesian(xlim=c(0, tmax * SCALING)) + 
        geom_vline(aes(xintercept=emean *SCALING + age), data=gpred) + 
        geom_vline(aes(xintercept=lmean * SCALING + age), lty=2, data=gpred) + 
        facet_wrap(~sample, scale="free_y", ncol=1, strip.position="left") +
        geom_rect(aes(xmin=0, ymin=0, ymax=Inf, xmax=age), fill='white')
}



high_chrom = function(chrom){
    chrom = as.character(chrom)
    chrom[chrom == 'X'] = 0
    chrom = as.integer(chrom)
    chrom[is.na(chrom)] = 0
    chrom[chrom>11] = 23 - chrom[chrom>11]
    chrom = as.factor(chrom)
    levels(chrom)[1] = 'X'
    chrom = factor(chrom, levels=c(1:11, 'X'))
}

