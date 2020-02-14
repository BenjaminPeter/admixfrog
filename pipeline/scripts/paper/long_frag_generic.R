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
        select(snp_id, bin, fixed, strict_fixed, depth, pos, p_est=p, one_of(pops), p_read) %>% 
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
        ggplot() + 
        #geom_col(aes(x=pos, y=-.05, fill=target)) + 
        geom_hline(yintercept=0, color="lightgrey") + 
        geom_hline(yintercept=1, color='white') + 

        geom_col(aes(x=pos, y=v, alpha=target, fill=fixed))  + 
        theme_classic() + 
        theme(axis.text.x=element_text(angle=90, vjust=.5, size=6), #strip.text.x=element_text(size=0),
              legend.position='none', panel.spacing.x=unit(0.05, 'lines'),
              strip.text.x=element_blank()) +
        scale_alpha_discrete(range=c(0.2, 1)) +
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

test1 <- function(){
    ref = read_csv("ref/ref_archaicadmixture.csv.xz")
    base="admixfrog/error/5000/AFR_NEA_DEN/Salkhit_archaicadmixture"
    base_ty="admixfrog/error/5000/AFR_NEA_DEN/Tianyuan_archaicadmixture"

    TARGET = 'AFRDEN'


    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base))
    
    run = runs %>% filter(target==TARGET) %>% arrange(-len) %>% head(1) %>% tail(1)
    P = plot_snp_run(snp, ref, run, ext=c(1e4, 1e4), filter_ambiguous=T, pops=c("AFR", "VIN", "DEN"))
}

test2 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz")
    base="admixfrog/error2CAFR/5000/NEA_DEN/denisova8_hcneaden"

    TARGET = 'NEADEN'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base))
    
    run = runs %>% filter(target==TARGET) %>% arrange(-len) %>% head(1) %>% tail(1)
    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     min_cov = 2,
                     ext=c(4e4, 4e4), filter_ambiguous=T, pops=c("AFR", "VIN", "DEN"), base_pop='DEN')
}

test3 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/altai_hcneaden"

    TARGET = 'NEADEN'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    
    run = runs %>% filter(target==TARGET) %>% arrange(-len) %>% head(1) %>% tail(1)
    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     min_cov = 50,
                     ext=c(4e4, 4e4), filter_ambiguous=T, pops=c("AFR", "VIN", "DEN"), base_pop='VIN')
}

test4 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CEUR/5000/A=ALT+D12_D=DEN+D11/denisova3_hcneaden"

    TARGET = 'AD'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, chrom != '6') %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 20,
                     ext=c(1e5, 1e5), filter_ambiguous=F, 
                     #pops=c("CHA", "VIN", "ALT", "D11", "DEN"), base_pop='ALT')
                     pops=c("ALT", "DEN"), base_pop='DEN')
}

paper_d11 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/denisova11_hcneaden"

    TARGET = 'NEA'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='homo' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     plot_coverage=T,
                     plot_est = T,
                     min_cov = 2,
                     large_top = 2,
                     one_snp_per_bin=F,
                     ext=c(4e5, 4e5), filter_ambiguous=T, 
                     pops=c("NEA", "DEN"), base_pop='NEA')
    ggsave("figures/paper/longest/d11_run.png", P, width=7.2, height=2)
}

paper_d3 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CEUR/5000/A=ALT+D12_D=DEN+D11/denisova3_hcneaden"

    TARGET = 'A'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom == '1' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    run = data.frame(chrom=6, pos=32798795, pos_end =36534853)
    run = data.frame(chrom=6, pos=32798795, pos_end =51823616)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 3,
                     ignore_bins=T,
                     large_top = 2,
                     one_snp_per_bin=F,
                     ext=c(1e6, 1e6), filter_ambiguous=F, 
                     pops=c("ALT", "DEN"), base_pop='DEN',
                     p_read_name='Denisova3') +
        theme(axis.text.x = element_blank())
    ggsave("figures/paper/longest/d3_run6.png", P, width=7.2, height=2)
}

paper_d3_c6 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CEUR/5000/A=ALT+D12_D=DEN+D11/denisova3_hcneaden"

    TARGET = 'A'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom == '6' , map_len > 0.1) 
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run,
                     filter_multi=F,
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 3,
                     ignore_bins=T,
                     large_top = 2.6,
                     one_snp_per_bin=F,
                     ext=c(5e6, 5e6), filter_ambiguous=F, 
                     pops=c("NEA", 'DEN'), base_pop='DEN',
                     p_read_name='Denisova3') +
        theme(axis.text.x = element_blank())
    ggsave("figures/paper/longest/d3_run_c6.png", P, width=7.2, height=1.75)
}

paper_d2 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="error2CAFR/20000/NEA_DEN/denisova2_hcneaden"


    TARGET = 'NEA'

    runs = read_csv(sprintf("rle/%s.rle0.4.xz", base))
    snp = read_csv(sprintf("admixfrog/%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     filter_fixed_strict=T,
                     plot_coverage=T,
                     plot_est = T,
                     ignore_bins=T,
                     min_cov = 1,
                     large_top = 2.5,
                     min_freq = 0.1,
                     one_snp_per_bin=F,
                     p_read_name='Denisova 2',
                     ext=c(2e6, 2e6), 
                     filter_ambiguous=T, 
                     pops=c("NEA", "DEN", "AFR"), base_pop='AFR')
    ggsave("figures/paper/longest/d2_run.png", P, width=7.2, height=2.5)
}

paper_d2_chr11 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="error2CAFR/5000/NEA_DEN/denisova2_hcneaden"


    TARGET = 'NEA'

    snp = read_csv(sprintf("admixfrog/%s.snp.xz", base), col_types=cols(chrom=col_character()))
    runs = read_csv(sprintf("rle/%s.rle0.5.xz", base))
    run = runs %>% filter(target==TARGET, type=='state' ) %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     filter_fixed_strict=T,
                     plot_coverage=T,
                     plot_est = T,
                     ignore_bins=T,
                     min_cov = 1,
                     large_top = 2.5,
                     min_freq = 0.1,
                     one_snp_per_bin=F,
                     p_read_name='Denisova 2',
                     ext=c(2e6, 2e6), 
                     filter_ambiguous=T, 
                     pops=c("NEA", "DEN", "AFR"), base_pop='AFR')
    ggsave("figures/paper/longest/d2_run11.png", P, width=7.2, height=2.5)
}
paper_d8 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/denisova8_hcneaden"

    TARGET = 'NEA'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom == 'X') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_fixed_strict=F,
                     filter_ambiguous=F, 
                     plot_coverage=T,
                     plot_est = F,
                     ignore_bins=T,
                     min_cov = 1,
                     large_top = 2.5,
                     min_freq = 0.1,
                     one_snp_per_bin=F,
                     p_read_name='Denisova 8',
                     ext=c(2e6, 2e6), 
                     pops=c("NEA", "DEN", "AFR"), base_pop='AFR') 
    ggsave("figures/paper/longest/d8_run.png", P, width=7.2, height=2.5)
}

paper_vindija <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/vindija3319_hcneaden"

    TARGET = 'NEA'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom=='X' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     filter_ambiguous=T, 
                     plot_coverage=T,
                     plot_est = T,
                     min_cov = 1,
                     large_top = 3,
                     min_freq = 0.1,
                     one_snp_per_bin=T,
                     ext=c(4e6, 3e6), 
                     pops=c("NEA", "DEN", "AFR"), base_pop='AFR')
    ggsave("figures/paper/longest/vindija_run.png", P, width=7.2, height=2.0)
}

paper_altai <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/altai_hcneaden"

    TARGET = 'DEN'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom!='X' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=T,
                     filter_ambiguous=T, 
                     plot_coverage=F,
                     plot_est = F,
                     min_cov = 1,
                     ignore_bins=T,
                     large_top = 2,
                     min_freq = 0.0,
                     one_snp_per_bin=F,
                     ext=c(5e5, 3e5), 
                     p_read_name = 'Denisova 5',
                     pops=c("VIN", "DEN"), base_pop='VIN') +
        theme(axis.text.x = element_blank())
    ggsave("figures/paper/longest/d5_run.png", P, width=7, height=1.75)
}

paper_spy <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/spy1_hcneaden"

    TARGET = 'DEN'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom!='X' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_ambiguous=T, 
                     plot_coverage=T,
                     plot_est = T,
                     min_cov = 1,
                     ignore_bins=F,
                     large_top = 2,
                     min_freq = 0.0,
                     one_snp_per_bin=F,
                     ext=c(1e5, 1e5), 
                     pops=c("NEA", "DEN"), base_pop='NEA')
    ggsave("figures/paper/longest/spy_run.png", P, width=7.2, height=2)
}

paper_mez1 <- function(){
    ref = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
    base="admixfrog/error2CAFR/5000/NEA_DEN/spy1_hcneaden"

    TARGET = 'DEN'

    runs = read_csv(sprintf("%s.rle.xz", base))
    snp = read_csv(sprintf("%s.snp.xz", base), col_types=cols(chrom=col_character()))
    run = runs %>% filter(target==TARGET, type=='state', chrom!='6' ) %>% arrange(-len) %>% head(1) %>% tail(1)
    #run = runs %>% filter(target==TARGET, chrom == '6') %>% arrange(-len) %>% head(1) %>% tail(1)

    P = plot_snp_run(snp, ref, run, 
                     filter_multi=F,
                     filter_ambiguous=T, 
                     plot_coverage=T,
                     plot_est = T,
                     min_cov = 1,
                     ignore_bins=F,
                     large_top = 2,
                     min_freq = 0.0,
                     one_snp_per_bin=F,
                     ext=c(1e5, 1e5), 
                     pops=c("VIN", "ALT", "CHA", "DEN"), base_pop='VIN')
    ggsave("figures/paper/longest/mez_run.png", P, width=7.2, height=2)
}

run_string <- function(run){
    require(scales)
    require(glue)
    run %>% select(chrom, pos, pos_end) %>% 
        summarize(pos=min(pos), pos_end=max(pos_end), chrom=first(chrom)) %>%
        mutate_all(.f=function(x)comma(as.integer(x))) %>% glue_data('chr{chrom}:{pos}-{pos_end}')
}
