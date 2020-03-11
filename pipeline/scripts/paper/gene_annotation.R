library(GenomicRanges)
library(tidyverse)
library(valr)
library(cowplot)


get_desert_overlap_true <- function(frags, bval){
    bval$is_introgressed = bval %>% 
        makeGRangesFromDataFrame(start.field='id', end.field='id') %>% 
        overlapsAny(frags %>% makeGRangesFromDataFrame)
    return(mean(bval$in_desert & bval$is_introgressed))
}

get_desert_overlap <- function(frags, bval){
    bval$is_introgressed = bval %>% 
        makeGRangesFromDataFrame %>% 
        overlapsAny(frags %>% makeGRangesFromDataFrame)
    return(mean(bval$in_desert & bval$is_introgressed))
}

get_bbin_overlap_true <- function(frags, bval){
    bval$is_introgressed = bval %>% 
        makeGRangesFromDataFrame(start.field='id', end.field='id') %>% 
        overlapsAny(frags %>% makeGRangesFromDataFrame)
    x = bval %>% filter(!is.na(bbin)) %>% group_by(bbin) %>% summarize(m=mean(is_introgressed))
    return(mutate(x, bbin=as.character(bbin)))
}

get_X_overlap_true <- function(frags, bval){
    bval$is_introgressed = bval %>% 
        makeGRangesFromDataFrame(start.field='id', end.field='id') %>% 
        overlapsAny(frags %>% makeGRangesFromDataFrame)
    x = bval %>% group_by(X=chrom=='X') %>% 
        summarize(m=mean(is_introgressed))
    return(x$m[x$X])
}

get_bbin_overlap <- function(frags, bval){
    bval$is_introgressed = bval %>% 
        makeGRangesFromDataFrame(keep=T) %>% 
        overlapsAny(frags %>% makeGRangesFromDataFrame)
    x = bval %>% filter(!is.na(bbin)) %>% group_by(bbin) %>% summarize(m=mean(is_introgressed))
    return(mutate(x, bbin=as.character(bbin)))
}

get_X_overlap <- function(frags, bval){
    bval$is_introgressed = bval %>% 
        makeGRangesFromDataFrame(keep=T) %>% 
        overlapsAny(frags %>% makeGRangesFromDataFrame)
    x = bval %>% group_by(X=chrom=='X') %>% 
        summarize(m=mean(is_introgressed))
    return(x$m[x$X])
}

test_bbin_overlap <- function(frags, genome, bval, n=100){
    res = data.frame(bbin=levels(bval$bbin))
    for(i in 1:n){
        #print(sprintf("%s/%s", i, n))
        ff = bed_shuffle(frags, genome)
        bbin_overlap <- get_bbin_overlap(ff, bval)
        names(bbin_overlap)[2] = sprintf("B%s", i)
        res <- left_join(res, bbin_overlap, by='bbin')
    }
    return(res)
}

test_desert_overlap <- function(frags, genome, bval, n=100){
    res = c()
    for(i in 1:n){
        #print(sprintf("%s/%s", i, n))
        ff = bed_shuffle(frags, genome)
        res <- c(res, get_desert_overlap(ff, bval))
    }
    return(res)
}

test_X_overlap <- function(frags, genome, bval, n=100){
    res = c()
    for(i in 1:n){
        ff = bed_shuffle(frags, genome)
        res <- c(res, get_X_overlap(ff, bval))
    }
    return(res)
}


#' load all identified introgression tracts
frags = read_csv("tables/paper/clean/tableS2_frags.csv.gz") %>% 
    rename(start=id, end=id_end) %>% tbl_interval

#' load all bins with b-values
bval_file = "tables/paper/bvals.csv.xz"
bval = read_csv(bval_file, col_types=cols(chrom=col_character()))  %>%
    group_by(chrom) %>%
    mutate(start=row_number(), end=row_number()) %>%
    ungroup %>% 
    mutate(bbin=cut(bval, quantile(bval, 0:5/5, na.rm=T), include.lowest=T))
genome = bval %>% group_by(chrom) %>% 
    summarize(min=min(id), max=max(id), size=max-min) %>%
    tbl_genome



set.seed(10)
n_resamples = 1000

obs_deserts <- get_desert_overlap_true(frags, bval)
obs_X <- get_X_overlap_true(frags, bval)
obs_bbin <- get_bbin_overlap_true(frags, bval) %>% 
    mutate(bbin=factor(bbin, levels=levels(bval$bbin))) %>%
    rename(obs=m)

background_deserts <- test_desert_overlap(frags, genome, bval, n=n_resamples)
background_bbin <- test_bbin_overlap(frags, genome, bval, n=n_resamples) %>%
    mutate(bbin=factor(bbin, levels=levels(bval$bbin))) %>%
    gather(resample, bval, -bbin) %>% select(-resample) %>% as_tibble  %>%
    left_join(obs_bbin)
background_X <- test_X_overlap(frags, genome, bval, n=n_resamples)


p_deserts <- mean(obs_deserts < background_deserts)
p_X <- mean(obs_X < background_X)


p_bbin <- background_bbin %>% group_by(bbin) %>% summarize(p = mean(obs < bval))
p_bbin %>% write_csv("tables/paper/p_bbin.csv")
p_deserts %>% as.data.frame %>%  write_csv("tables/paper/p_deserts.csv")
p_X %>% as.data.frame %>%  write_csv("tables/paper/p_X.csv")


df_deserts = data.frame(resample=background_deserts, obs=obs_deserts) %>% as_tibble
P_deserts = df_deserts %>% ggplot(aes(y=resample)) + 
    geom_boxplot() + 
    geom_boxplot(aes(y=obs), color='red', lty=2) +
    theme_classic(10) +
    ylab("desert introgressed (cM)") +
    theme(axis.text.x=element_blank())

P_bbin = background_bbin %>% ggplot(aes(x=bbin, y=bval)) + 
    geom_boxplot() + 
    geom_boxplot(aes(x=bbin, y=obs), color='red') +
    theme_classic(9) +
    scale_x_discrete(NULL) + 
    scale_y_continuous("Pr(introgressed)") + 
    theme(axis.text.x=element_text(angle=90))

P = plot_grid(P_deserts, P_bbin, rel_widths=c(2, 5), align='h', label_size=9) 
ggsave("figures/paper/s_resample.png", P, width=3.5, height=2.5)


