require(tidyverse)
require(reshape2)
source("scripts/lomax.R")

bin_size = 10000
panel = "AFR_NEA_CHA_ALT"
samples <- c("UstIshim", "Kostenki14", "Tianyuan", "Salkhit", "Yana1", "Yana2", 
             "Kolyma", "Malta", "LBK", "Sunghir3", "Sunghir4")
samples <- c("UstIshim", "Tianyuan", "Salkhit", "Yana1", "Yana2",  "Oase",
             "Kolyma",  "LBK", "Sunghir3", "Sunghir4", "Sunghir1", "Kostenki14sg")
ascertainment <- "archaicadmixture"

#bin_size = 10000
#panel = "ALT_NEA_DEN"
#samples <- sprintf("Denisova%s", c("2", "3", "4", "8", "11", "13"))
#ascertainment <- "allsites"


#bin_size=20000
#panel="ALT_NEA_CHA"
#samples <- sprintf("Chagyrskaya%s", c("02", "07", "09", "12", "13", "17", "19", "41", "60"))
#samples <- c(samples, "vindija3316","vindija3325", "vindija3326", 
#             "vindijag1")
#ascertainment = "allsites"


#infiles = sprintf("data/%s/%d/%s/%s_nofilter.bin.out.csv.gz", samples, bin_size, panel, ascertainment)
infiles = sprintf("data/%d/%s/%s_%s_nofilter.bin.out.csv.gz",  bin_size, panel, samples,ascertainment)
snpfiles = sprintf("data/%d/%s/%s_%s_nofilter.snp.out.csv.gz",  bin_size, panel, samples,ascertainment)

read_binout <- function(fname){
    a <- read_csv(fname) %>% select(-X1) %>% mutate(chrom=as.integer(chrom_id+1))
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

load_data <- function(){
    a <- lapply(infiles, read_binout)
    names(a) <- samples
    a <- bind_rows(a, .id="sample")
    a %>% mutate(viterbi=as.factor(names(a)[viterbi+5]))
}

get_long_data <- function(data){
    b <- a %>% select(-chrom_id) %>% 
        melt(id.vars=c("sample", "chrom", "bin_id", "bin_pos",  "viterbi")) %>% 
        as_tibble 
}

basic_plot <- function(a, b){
    v <- a %>% 
        select(-chrom, -chrom_id, -bin_pos, -bin_id, -sample, -viterbi)  %>% 
        summarize_all(.funs=mean)  %>% 
        unlist
    lvl <- names(v[v<.5& v>1e-5])                                                           

    b  %>% filter(chrom==1, bin_pos > 223e6, bin_pos < 231e6, variable %in% lvl) %>% 
        ggplot(aes(x=bin_pos/1e6, y=value, color=variable, fill=variable)) + 
        geom_col(width=bin_size/1e6) + 
        facet_wrap(~chrom*sample, ncol=1, strip.position="left")
}

rle_plot <- function(data){
    df <- data %>% group_by(sample) %>%
        arrange(sample, chrom) %>%
        #do(rle=get_rle(.$viterbi!="AFR"))
        do(rle=get_rle(.$AFR < .3))
    df %>% unnest(rle) %>%
        filter(values) %>%
        ggplot(aes(x=lengths * bin_size / 1e6)) + 
        stat_bin(binwidth=2*bin_size/1e6, geom="point") +
        stat_density(dth=2*bin_size/1e6, geom="point") +
        facet_wrap(~sample) +
        scale_y_log10(name="# fragments")# + 
        #scale_x_log10(name="length (Mb)")
}

olo=function(l, trunc=4){
    o = tlomax_opt2(l, trunc)
    pars =(optim(c(4, 4), o, method="L-BFGS-B", lower=c(0, 1e-2))$par)
    return(data.frame(scale=pars[1], shape=pars[2]))
}
oe=function(l, trunc=4){
    o = texp_opt(l, trunc)
    par =exp(optim(c(-1), o, method="L-BFGS-B")$par)
    return(data.frame(rate=par))
}

rle_fit_pars <- function(data, trunc=4){
    df <- data %>% group_by(sample) %>%
        arrange(sample, chrom) %>%
        #do(rle=get_rle(.$viterbi!="AFR"))
        do(rle=get_rle(.$AFR < .3)) %>%
        unnest(rle) %>% 
        filter(values, lengths > trunc) %>%
        group_by(sample)

    x1 = df %>% do( p2=oe(.$lengths, trunc=trunc)) %>% unnest(p2) %>% 
        mutate(emean=rate * 1e4)
    x2 = df %>% do( p2=olo(.$lengths, trunc=trunc)) %>% unnest(p2) %>%
        mutate(lmean=(shape-1)/scale * 1e4)

    return(inner_join(x1, x2))
}
rle_fit_plot <- function(data, R, lmax=100){
    df <- data %>% group_by(sample) %>%
        arrange(sample, chrom) %>%
        #do(rle=get_rle(.$viterbi!="AFR"))
        do(rle=get_rle(.$AFR < .3)) %>% 
        unnest(rle) %>%
        filter(values) %>%
        filter(lengths >= 4, lengths < lmax) %>%
        group_by(sample, lengths) %>% 
        tally %>%
        mutate(d=n / sum(n))

    s = seq(4, lmax, .1)
    lpred = R %>% rowwise %>% 
        do(l = data.frame(lengths=s, lomax=dtlomax(s, shape=.$shape, scale=.$scale, trunc=4, log=F)))
    epred = R %>% rowwise %>% 
        do(e = data.frame(lengths=s, exp=dtexp(s, rate=.$rate, trunc=4, log=F)))

    pred = inner_join(unnest(bind_cols(R, lpred)), unnest(bind_cols(R, epred)))
    df1 = df %>% select(sample, lengths, d) %>% melt(id.vars=c("sample", "lengths"))                     
    df2 = pred %>% select(sample, lengths, exp, lomax) %>% melt(id.vars=c("sample", "lengths"))          
                                                                                                         
    df3 = bind_rows(df1, df2) %>% as_tibble                                                              

     
    ggplot() + 
        geom_point(data=df1, aes(x=lengths, y=value)) + 
        geom_line(data=df2, aes(x=lengths, y=value, color=variable)) + 
        facet_wrap(~sample) + 
        xlab("Lengths (x10kb)") + 
        ylab("# fragments") 
}

plot_m_gamma <- function(R){
    t <- seq(1, 1000) / 4000
    gpred = R %>% rowwise %>% 
        do(l = data.frame(time=t, mig=dgamma(t, shape=.$shape, scale=1/.$scale)))

    bind_cols(gpred, R) %>% 
        unnest(l) %>% 
        ggplot( aes(x=time*1e4*29, y=mig, color=sample)) + 
        geom_line(lwd=2) + 
        #coord_cartesian(ylim=c(0,20), xlim=c(0,1000)) + xlab("time (gen)") + 
        coord_cartesian(ylim=c(0,20), xlim=c(0,30000)) + xlab("time (y)") + 
        geom_vline(aes(xintercept=rate*1e4*29, color=sample), data=R) + facet_wrap(~sample)
}

get_rundf <- function(data, min_length=0){
    df <- data %>% group_by(sample) %>%         
        arrange(sample, chrom) %>%              
        do(rle=get_rle(.$AFR < .3)) %>%         
        unnest(rle) %>%
        mutate(values=ifelse(values * lengths > min_length, TRUE, FALSE))

    df_run = df %>% group_by(sample) %>% 
        do(runs=data.frame(bin_id=1:sum(.$lengths)-1, 
                           runs=rep(.$values, .$lengths))) %>% 
        unnest(runs) %>% 
        left_join(data) %>%
        select(sample, bin_id, runs) %>%
        dcast(bin_id ~ sample) %>%
        as_tibble
}

#r200 = get_rundf(data, 200)                                                           
#CC = cov(r200[,-c(1,5)]); diag(CC) <- NA                                              
#heatmap.2(CC, symm=T, trace='n', breaks=100, col="viridis", scale="none", symbreaks=F)


