library(tidyverse)
source("scripts/plotting/lib.R")
library(yaml)

#' run is a data frame with chrom, pos and pos_end
#' ext: how many bp to add left/right
#'



plot_snp_run <- function(snp, ref, run, ext=c(4e5, 5e4),
                     filter_ambiguous=T,
                     filter_multi=F,
                     filter_fixed_strict=F,
                     plot_coverage=T,
                     large_top=1,
                     ignore_bins=F,
                     plot_est=T,
                     min_cov = 1,
                     min_freq = 0,
                     one_snp_per_bin = F,
                     pops = c("AFR", "VIN", "DEN"),
                     dont_plot_pop=c(),
                     base_pop= 'DEN', #for fixed,
                     p_read_name = 'p_read'
                     ){ 
    run_snp = snp %>% 
            filter(chrom==min(run$chrom), pos>=min(run$pos) - ext[1], 
                   pos < max(run$pos_end) + ext[2]) %>%
            left_join(ref %>% select(-map)) 
    
    a1 = c(sprintf("%s_alt", pops), "talt")
    a2 = c(sprintf("%s_ref", pops), "tref")
    freqs = lapply(1: (length(pops) + 1), 
                   function(i) run_snp[a1[i]] / (run_snp[a1[i]] + run_snp[a2[i]])) %>% 
        bind_cols %>% 
        as.data.frame
    names(freqs) = c(pops, 'p_read')
    freqs = freqs %>% 
        bind_cols(run_snp %>% select(p)) %>%
        mutate_all(function(x)ifelse(.[,base_pop] > 0.5, 1-x, x)) %>%
        mutate_all(function(x)ifelse(.[,base_pop] == .5 & rowSums(.[,1:(ncol(.)-2)], na.rm=T) > 0.5, 1-x, x))
    lsr0 = bind_cols(run_snp %>% select(-p), freqs)

    lsr0 = lsr0 %>% filter(p_read >= min_freq)
    

    #TRIALS
    is_fixed = sapply(pops, function(p)p=abs(lsr0[base_pop] - lsr0[p]) == 1,
                      simplify=F, USE.NAMES=T)
    lsr0$fixed = 'none'
    for(p in pops) lsr0$fixed[is_fixed[[p]]] = sprintf('fixed in %s' ,p)
    lsr0$fixed[is_fixed %>% bind_rows %>% rowSums(na.rm=T) > 1] = 'multiple'

    lsr0$strict_fixed = rowSums(lsr0[,pops] == 0) == length(pops) - 1 & rowSums(lsr0[,pops] == 1) == 1

    
    x= lsr0 %>% 
        mutate(depth=tref+talt) %>% 
        filter(depth >= min_cov) %>%
        select(snp_id, bin, fixed, strict_fixed, depth, pos, p_est=p, all_of(pops), p_read) %>% 
        gather('k', 'v', c('p_read', pops), factor_key=T) 

    x$target=apply(run[,c('pos', 'pos_end')], 1, 
                   function(row)row[1] < x$pos & x$pos < row[2])  %>% apply(1, any)
    #    mutate(target=pos>run$pos & pos < run$pos_end,
    #           pos=scales::comma(pos)) %>%
    x = x %>% mutate(pos=scales::comma(pos)) %>%
        filter(!is.na(p_est))

    if(filter_ambiguous){
        x = x %>% filter(!fixed %in% c('none'))
    }

    if(filter_multi){
        x = x %>% filter(!fixed %in% c('multiple'))
    }

    if(filter_fixed_strict){
        x = x %>% filter(strict_fixed)
    }

    if(one_snp_per_bin){
        x = x %>% group_by(bin) %>%
            filter(snp_id == min(snp_id))
    }

    if(large_top){
        x = x %>% mutate(
                         p_est = p_est * large_top,
                         v = ifelse(k=='p_read', v*large_top, v)
                         )
    }


    cols = yaml::read_yaml("colors.yaml")$colors %>% unlist
    names(cols) = sprintf("fixed in %s", names(cols))
    cols = c(cols, 'none'='darkgrey', 'multiple'='darkgrey')

    levels(x$k)[levels(x$k)=='p_read'] = p_read_name

    P = x %>% 
        filter(! k %in% dont_plot_pop) %>%
        ggplot() + 
        geom_hline(yintercept=0, color="lightgrey") + 
        geom_hline(yintercept=1, color='white') + 
        geom_col(aes(x=pos, y=v, alpha=target, fill=fixed))  + 
        theme_classic() + 
        theme(axis.text.x=element_text(angle=90, vjust=.5, size=6), 
              legend.position='none', panel.spacing.x=unit(0.05, 'lines'),
              strip.text.x=element_blank()) +
        scale_alpha_discrete(range=c(0.3, 1)) +
        scale_x_discrete("SNP position") +
        scale_y_continuous("frequency", labels = scales::rescale_max, 
                           breaks=function(x){c(mean(x), x[2])}, 
                                                     expand=c(0,0)) +
        scale_fill_manual(values=cols)


    if(ignore_bins){
        P = P +facet_grid(k~., space='free', scales='free', switch='x')
    } else {
        P = P + facet_grid(k~bin, space='free', scales='free', switch='x')
    } 
    if(plot_coverage){
        P = P + geom_text(aes(x=pos, y=large_top, label=depth), 
                          vjust=1, size=3, color='black', data=x %>% filter(k==p_read_name))
    } else {
        P = P + geom_hline(aes(yintercept=large_top), alpha=0, data = x %>% filter(k==p_read_name))
    }

    if(plot_est){
        P = P + geom_area(aes(x=pos, y=p_est, alpha=target), lwd=.3, lty=1, color='black', data=x %>% filter(k==p_read_name))
    }
    return(P)
}





run_string <- function(run){
    require(scales)
    require(glue)
    run %>% select(chrom, pos, pos_end) %>% 
        summarize(pos=min(pos), pos_end=max(pos_end), chrom=first(chrom)) %>%
        mutate_all(.f=function(x)comma(as.integer(x))) %>% glue_data('chr{chrom}:{pos}-{pos_end}')
}
