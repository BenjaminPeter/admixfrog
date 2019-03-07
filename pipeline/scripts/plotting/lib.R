suppressPackageStartupMessages({
library(tidyverse)
library(VGAM)
source("scripts/plotting/meancol.R")
})
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

bg_chrom <- function(ref=NULL, map="AA_Map"){
    x = read_csv("ref/chrom_limits.csv", col_types=cols(chrom=col_factor())) %>% 
        select('chrom', min=sprintf('%s_min', map), max=sprintf("%s_max", map))
    if(!is.null(ref)) x = filter(x, chrom %in% unique(ref$chrom))

    return(geom_rect(data=x, mapping=aes(xmin=min, xmax=max, ymin=-Inf, ymax=Inf), fill='#efefef'))

}

read_run <- function(rfile){
    read_csv(rfile, col_types='ddddd') %>%
	    group_by(state, len) %>%
	    summarize(n=sum(n))
}

load_runs_data <- function(files, name){
    a <- lapply(files, read_run)
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

rle_plot_map <- function(data, minlen=0.1, maxlen=1){
    data %>%
	filter(map_len > minlen) %>%
	ggplot() +
	bg_chrom(data) + 
	geom_tile(aes(x =(map +map_end)/ 2, width = map_end - map, 
		   y=0.5, height=1, fill=target, alpha=pmin(map_len,maxlen))) + 
	facet_wrap(~chrom, ncol=1, strip='left') + 
	THEME + 
	col_scale() +
	scale_alpha_continuous(range=c(0.3,1), limits=c(minlen, maxlen), name='Length(cM)') + 
	scale_x_continuous(expand=c(0,0), name='Position (cM)') + 
	scale_y_discrete(expand=c(0,0), name='state')
}


fit_lomax=function(l, trunc=0.01){
    try({
    o = tlomax_opt2(l[l>trunc], trunc)
    res = optim(c(4, 4), o, method="L-BFGS-B", lower=c(1e-5, 1))
    pars = res$par
    return(data.frame(scale=pars[1], shape=pars[2], ll_lomax=res$value))
    })
    return(data.frame(scale=NA, shape=NA, ll_lomax=NA))
}

fit_exp=function(l, trunc=4){
    try({
    o = texp_opt(l, trunc)
    res = optim(c(-1), o, method="L-BFGS-B")
    par =exp(res$par)
    return(data.frame(rate=par, ll_exp = res$value))
    })
    return(data.frame(rate=NA, ll_exp=NA))
}
#' lomax copy stufff 
dtlomax <- function(x, scale, shape, trunc, log=T){
    if(log){
        k <- VGAM::dlomax(x, scale=scale, shape=shape, log=T) - 
        VGAM::plomax(trunc, scale=scale, shape=shape, lower=F, log=T)
        k[x<trunc] <- -Inf
        return(k)
    } else {
        print("lomax normal")
        k <- VGAM::dlomax(x, scale=scale, shape=shape, log=F) / 
        VGAM::plomax(trunc, scale=scale, shape=shape, lower=F, log=F)
        k[x<trunc] <- NA#0
        return(k)
    }
}

texp_opt <- function(x, trunc=3){
    f <- function(par){
        -sum(dtexp(x, exp(par[1]), trunc=trunc))
    }
    return(f)
}
dtexp <- function(x, rate=1, trunc=0, log=T){
    x = x[x>=trunc]
    if(log){
	dexp(x, rate, log=T) - pexp(trunc, rate, lower=F, log=T)
    } else {
	dexp(x, rate, log=F) / pexp(trunc, rate, lower=F, log=F)
    }
}
tlomax_opt2 <- function(x, trunc=0){
    f <- function(par){
        x = x[x>=trunc]
        a=-sum(dtlomax(x, scale=par[1], shape=par[2], trunc=trunc))
	return(a)
    }
}

rle_fit_pars <- function(data, trunc=0.05){
    df <- data %>% group_by(sample) 

    x1 = df %>% do( p2=fit_exp(.$map_len, trunc=trunc)) %>% 
	    unnest(p2) %>%
	    mutate(emean=rate)
    x2 = df %>% do( p2=fit_lomax(.$map_len, trunc=trunc)) %>% 
	    unnest(p2) %>%
	    mutate(lmean=(shape-1)/scale)

    return(inner_join(x1, x2) %>% mutate(delta_ll = ll_lomax - ll_exp))
}

rle_fit_plot <- function(data, R, trunc=0.049, xmax=6){
    data2 = data %>% filter(map_len > trunc)
    P = data2 %>%
        ggplot(aes(x=map_len, group=sample)) +
        geom_point(aes(y=1-..y..), stat='ecdf', pad=F) +
        scale_y_log10(expand=expand_scale(0,0), name='Quantile') +
        scale_x_continuous(expand=expand_scale(0,0), name = "Length (cM)") +
        coord_cartesian(xlim=c(trunc, xmax), ylim=c(1e-3, 1)) +
        facet_wrap(~sample, scale='free') 

    s = seq(trunc, max(data$map_len), .01)
    lpred = R %>% rowwise %>% 
        do(l = data.frame(lengths=s, 
			  lomax=plomax(s, shape=.$shape, scale=.$scale, lower=F)/
                  plomax(trunc, shape=.$shape, scale=.$scale, lower=F)))
    epred = R %>% rowwise %>% 
        do(e = data.frame(lengths=s, 
			  exp=exp(-(s - trunc) * .$rate)))

    pred = inner_join(unnest(bind_cols(R, lpred)), unnest(bind_cols(R, epred))) %>%
        select(lengths, sample, lomax, exp) %>%
        gather(k, v, lomax, exp)

     
    P + geom_path(data=pred, mapping=aes(color=k, x=lengths, y=v, group=sample))
}
plot_m_gamma <- function(R, generation_time){
    SCALING = generation_time * 100 #* 2
    tmax = pmin(150000 / SCALING, max(R$semean, R$slmean, na.rm=T))
    t <- seq(0, tmax, 0.01)
    n_samples = length(unique(R$sample))


    gpred = R %>% 
        rowwise %>% 
        do(l = data.frame(time=t+.$sage, 
			  mig=dgamma(t, shape=.$shape, scale=1/.$scale))) %>%
        bind_cols(R) %>% 
        filter(!is.na(lmean + emean), slmean < tmax * 1.1) %>%
        unnest(l) 

    M = gpred %>%
        ggplot( aes(x=time*SCALING, y=mig)) + 
        geom_line(lty=2) + 
        #coord_cartesian(ylim=c(0,20), xlim=c(0,1000)) + xlab("time (gen)") + 
        scale_x_continuous(name="time (y)", expand=expand_scale(0,0)) + 
        coord_cartesian(xlim=c(0, tmax * SCALING)) + 
        geom_vline(aes(xintercept=semean * SCALING), data=gpred) + 
        facet_wrap(~sample, scale="free_y", ncol=1, strip.position="left") +
        geom_rect(aes(xmin=0, ymin=0, ymax=Inf, xmax=age), fill='white')
}
