require(VGAM)
LOWER=1e-9
require(tidyverse)
require(reshape2)
#source("lomax.R")


read_binout <- function(fname){
    a <- read_csv(fname) %>% mutate(chrom=factor(chrom, levels=unique(chrom)))
}
read_snpout <- function(fname){
    a <- read_csv(fname) %>% select(-X1) %>% mutate(chrom=as.integer(chrom_id+1))
}


to_wide <- function(fname){
    b <- a %>% select(-chrom_id) %>% 
        melt(id.vars=c("chrom", "bin_id", "bin_pos", "viterbi")) %>% 
        as_tibble %>% 
        mutate(viterbi=as.factor(names(a)[viterbi+5]))
}


get_rle <- function(cond, smooth=30000/bin_size){
    r <- rle(cond); 
    r$values[r$lengths <=smooth & !r$values] <- T; 
    r <- rle(rep(r$values, r$lengths))

    df = tibble(starts=cumsum(c(0,r$lengths[-length(r$lengths)]))+1,
                ends = cumsum(r$lengths),
                lengths = r$lengths,
                values = r$values)

    return(df)
}

load_data <- function(infiles, name){
    a <- lapply(infiles, read_binout)
    names(a) <- name
    a <- bind_rows(a, .id="sample")
}

get_long_data <- function(data){
    b <- data %>% 
        select(-id, -chrom_id, -hap, -n_snps) %>%
        melt(id.vars=c("sample", "chrom", "map", "pos", "pwidth",  "viterbi")) %>% 
        as_tibble 
}

basic_plot <- function(a, b, lvl=NULL, p_max = .5, p_min = 5e-3){
    if(is.null(lvl)){
    v <- a %>% 
        select(-chrom, -chrom_id, -bin_pos, -bin_id, -sample, -viterbi)  %>% 
        summarize_all(.funs=mean)  %>% 
        unlist 
    lvl <- names(v[v<p_max& v>p_min])                                                           
		} 

    b = b  %>% filter( variable %in% lvl) 
    ggplot() + 
	geom_line(data=b, aes(x=bin_pos/1e6, y=value, color=variable, fill=variable), lwd=.3) + 
        #geom_col(width=bin_size/1e6) + 
        facet_wrap(~chrom*sample, ncol=2, strip.position="left")
}

rle_plot <- function(data){
    df <- data %>% group_by(sample) %>%
        arrange(sample, chrom) %>%
        do(rle=get_rle(.$TRACK)) 
    df %>% unnest(rle) %>%
        filter(values) %>%
        ggplot(aes(x=lengths * bin_size / 1e6)) + 
        stat_bin(binwidth=2*bin_size/1e6, geom="point") +
        facet_wrap(~sample) +
        scale_y_log10(name="# fragments")# + 
}

olo=function(l, trunc=4){
    try({
    o = tlomax_opt2(l, trunc)
    pars =(optim(c(4, 4), o, method="L-BFGS-B", lower=c(0, 1e-2))$par)
    return(data.frame(scale=pars[1], shape=pars[2]))
    })
    return(data.frame(scale=NA, shape=NA))
}
oe=function(l, trunc=4){
    try({
    o = texp_opt(l, trunc)
    par =exp(optim(c(-1), o, method="L-BFGS-B")$par)
    return(data.frame(rate=par))
    })
    return(data.frame(rate=NA))
}

rle_fit_pars <- function(data, trunc=4){
    df <- data %>% group_by(sample) %>%
        arrange(sample, chrom) %>%
        do(rle=get_rle(.$TRACK)) %>%         
        unnest(rle) %>% 
        filter(values, lengths > trunc) %>%
        group_by(sample)

    x1 = df %>% do( p2=oe(.$lengths, trunc=trunc)) %>% unnest(p2) %>% 
        mutate(emean=rate * 1e4)
    x2 = df %>% do( p2=olo(.$lengths, trunc=trunc)) %>% unnest(p2) %>%
        mutate(lmean=(shape-1)/scale * 1e4)

    return(inner_join(x1, x2))
}
rle_fit_plot <- function(data, R){
    df <- data %>% group_by(sample) %>%
        arrange(sample, chrom) %>%
        do(rle=get_rle(.$TRACK)) %>%         
        unnest(rle) %>%
        filter(values) %>%
        filter(lengths >= 4) %>%
        group_by(sample, lengths) %>% 
        tally %>%
        mutate(d=n / sum(n))

    s = seq(4, 50e6/bin_size, .1)
    lpred = R %>% rowwise %>% 
        do(l = data.frame(lengths=s, 
			  lomax=dtlomax(s, shape=.$shape, scale=.$scale, trunc=4, log=F)))
    epred = R %>% rowwise %>% 
        do(e = data.frame(lengths=s, 
			  exp=dtexp(s, rate=.$rate, trunc=4, log=F)))

    pred = inner_join(unnest(bind_cols(R, lpred)), unnest(bind_cols(R, epred))) 
    df1 = df %>% select(sample, lengths, d) %>% 
	melt(id.vars=c("sample", "lengths")) 

    df2 = pred %>% select(sample, lengths, exp, lomax) %>% 
	melt(id.vars=c("sample", "lengths")) %>%
	filter(value > min(df1$value))
                                                                                                         
    #df3 = bind_rows(df1, df2) %>% as_tibble                                                              

     
    ggplot() + 
        geom_point(data=df1, aes(x=lengths * bin_size / 1e6, y=value)) + 
        geom_line(data=df2, aes(x=lengths * bin_size / 1e6, y=value, color=variable), lwd=1) + 
        facet_wrap(~sample, scale="free") + 
        xlab("Lengths (Mb)") + 
        ylab("# fragments") 
}

plot_m_gamma <- function(R, bin_size, generation_time){
    t <- seq(1, 100000) / 4000
    gpred = R %>% rowwise %>% 
        do(l = data.frame(time=t, 
			  mig=dgamma(t, shape=.$shape, scale=1/.$scale))) %>%
	bind_cols(R) %>% 
        unnest(l) %>% 
	filter(mig>1e-4) %>%
        ggplot( aes(x=time*bin_size*generation_time, y=mig, color=sample)) + 
        geom_line(lty=2) + 
        #coord_cartesian(ylim=c(0,20), xlim=c(0,1000)) + xlab("time (gen)") + 
        xlab("time (y)") + 
        geom_vline(aes(xintercept=rate*bin_size*generation_time, color=sample), data=R) + 
	facet_wrap(~sample, scale="free_y", ncol=1, strip.position="left") +
	coord_cartesian(xlim=c(0, 1e5))
}

get_rundf <- function(data, min_length=0){
    df <- data %>% group_by(sample) %>%         
        arrange(sample, chrom) %>%              
        do(rle=get_rle(.$TRACK)) %>%         
        unnest(rle) %>%
        mutate(values=ifelse(values * lengths > min_length, TRUE, FALSE))

    df_run = df %>% group_by(sample) %>% 
        do(runs=data.frame(id=1:sum(.$lengths)-1, 
                           runs=rep(.$values, .$lengths))) %>% 
        unnest(runs) %>% 
        left_join(data, iby=c("sample", "id")) %>%
        select(sample, id, runs) %>%
        dcast(id ~ sample) %>%
        as_tibble
}

#r200 = get_rundf(data, 200)                                                           
#CC = cov(r200[,-c(1,5)]); diag(CC) <- NA                                              
#heatmap.2(CC, symm=T, trace='n', breaks=100, col="viridis", scale="none", symbreaks=F)


#' lomax copy stuff
dtlomax <- function(x, scale, shape, trunc, log=T){
    if(log){
	k <- VGAM::dlomax(x, scale=scale, shape=shape, log=T) - 
	VGAM::plomax(trunc, scale=scale, shape=shape, lower=F, log=T)
	k[x<trunc] <- -Inf
	return(k)
    } else {
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
        a=-sum(dtlomax(x, scale=par[1], shape=par[2], trunc=trunc))
	return(a)
    }
}
