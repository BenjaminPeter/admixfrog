#I want full-page figure for ~20 samples:
# - 1 row is ~1cM
# panels are
#     1. prop archaic (2cM)
#     2. prop cont (2cM)
#     3. run-length-distribution (4 cm)
#     4. est. age (3 cm)
#     5. sample fragment (4 cm)



library(tidyverse)
library(cowplot)
library(reshape2)
source("scripts/plotting/lib.R")
save.image("salkhit.rdebug")
infiles = snakemake@input$bins
bin_size = as.numeric(snakemake@wildcards$bin_size) / 1e6
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
ref_file=snakemake@input$ref

TARGET = snakemake@wildcards$target
STATE = 2


n_snps = 1349147
col_lomax = 'blue'
col_exp = 'red'
TRUNC = .2
G_TIME = 29
SCALING = G_TIME * 100



#load data

snp = load_snp_data(snakemake@input$snp, names) 

snp2 = snp %>% select(sample, chrom, pos, p) %>% spread(sample, p)

frags = read_csv(snakemake@input$frags)


# restrict fragments to fragments found in Salkhit
frags2 = frags %>% filter(sample=="Salkhit") %>% 
    select(chrom, frag) %>% filter(!is.na(frag)) %>% 
    unique %>% left_join(frags)  %>% 
    mutate(n_snps=round(score/nscore)) %>% 
    select(sample, chrom, frag, start, end, map_len, len, n_snps) %>% 
    arrange(frag, -map_len)  
frags2ty = frags %>% filter(sample=="Tianyuan") %>% 
    select(chrom, frag) %>% filter(!is.na(frag)) %>% 
    unique %>% left_join(frags)  %>% 
    mutate(n_snps=round(score/nscore)) %>% 
    select(sample, chrom, frag, start, end, map_len, len, n_snps) %>% 
    arrange(frag, -map_len)  

#get bins overlapping each fragment
frags2= frags %>% filter(!is.na(frag)) %>%
    group_by(frag) %>%
    mutate(start=min(start), end=max(end)) %>%
    select(frag, sample, start, end)
bin_tol = 100
frags3 = frags2 %>% rowwise %>% do(bin=seq(.$start-bin_tol,.$end+bin_tol)) %>%  
    bind_cols(frags2, .) %>%
    left_join(frags %>% select(sample, frag, chrom, map_len)) %>%
    mutate(chrom=as.factor(chrom)) %>% unnest(bin)

frags4 =frags3 %>% 
    inner_join(select(snp, sample, bin, chrom, pos, tref, talt, p), by=c("bin", 'chrom', 'sample'))


ref= read_csv(ref_file, col_types=cols(chrom=col_factor()))
ref = ref %>% mutate(AFR = AFR_alt / (AFR_ref + AFR_alt),
                     NEA = NEA_alt / (NEA_ref + NEA_alt),
                     DEN = DEN_alt / (DEN_ref + DEN_alt)) %>%
    select(chrom, pos, ref, alt, AFR, NEA, DEN)

frags5 = frags4 %>% filter(!is.na(pos)) %>% left_join(ref, by=c('chrom', 'pos'))

write_csv(frags5, "tables/denisova_frags.csv.gz")

frags6 = frags5 %>% 
    filter(AFR != NEA | DEN != NEA) %>%
    mutate(den1= ifelse(DEN > AFR & DEN > NEA, talt, 0),
           den2= ifelse(DEN < AFR & DEN < NEA, tref, 0),
           nea1= ifelse(NEA > AFR & NEA > DEN, talt, 0),
           nea2= ifelse(NEA < AFR & NEA < DEN, tref, 0),
           afr1= ifelse(AFR > DEN & AFR > NEA, talt, 0),
           afr2= ifelse(AFR < DEN & AFR < NEA, tref, 0)) %>%
    mutate(ad1= ifelse(DEN == AFR & DEN > NEA, talt, 0),
           ad2= ifelse(DEN == AFR & DEN < NEA, tref, 0),
           an1= ifelse(NEA == AFR & AFR > DEN, talt, 0),
           an2= ifelse(NEA == AFR & AFR < DEN, tref, 0),
           nd1= ifelse(NEA == DEN & NEA > AFR, talt, 0),
           nd2= ifelse(NEA == DEN & NEA < AFR, tref, 0)) %>%
    mutate(afr= afr1+afr2, nea=nea1+nea2, den=den1+den2,
           ad = ad1 + ad2, nd=nd1+nd2, an=an1 +an2,
           A=afr+an, D=den+nd, N=nea, O=ad) %>%
    select(sample, frag, chrom, pos, A, D, N, O)
    #select(sample, frag, chrom, pos, afr, nea, den, ad, nd, an, AFR, NEA ,DEN, tref, talt)



frag_plot <- function(frag_id){
    frags6 %>% filter(frag==frag_id) %>% 
        gather(pop, v, A:O) %>% filter(v>0) %>% 
        ggplot(aes(x=factor(pos), y=v, fill=pop)) + 
        geom_col() + facet_wrap(~sample, ncol=1) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = .5, vjust=.5)) +
        scale_fill_manual(values=c( "#89bff7", "#ffb000","#444444","#efefef"))
}
frag_plot2 <- function(frags, frag_id, w=1000){
    frags6 = frags %>% 
        filter(AFR != NEA | DEN != NEA) %>%
        mutate(den1= ifelse(DEN > AFR & DEN > NEA, talt, 0),
               den2= ifelse(DEN < AFR & DEN < NEA, tref, 0),
               nea1= ifelse(NEA > AFR & NEA > DEN, talt, 0),
               nea2= ifelse(NEA < AFR & NEA < DEN, tref, 0),
               afr1= ifelse(AFR > DEN & AFR > NEA, talt, 0),
               afr2= ifelse(AFR < DEN & AFR < NEA, tref, 0)) %>%
        mutate(ad1= ifelse(DEN == AFR & DEN > NEA, talt, 0),
               ad2= ifelse(DEN == AFR & DEN < NEA, tref, 0),
               an1= ifelse(NEA == AFR & AFR > DEN, talt, 0),
               an2= ifelse(NEA == AFR & AFR < DEN, tref, 0),
               nd1= ifelse(NEA == DEN & NEA > AFR, talt, 0),
               nd2= ifelse(NEA == DEN & NEA < AFR, tref, 0)) %>%
        mutate(afr= afr1+afr2, nea=nea1+nea2, den=den1+den2,
               ad = ad1 + ad2, nd=nd1+nd2, an=an1 +an2,
               A=afr+an, D=den+nd, N=nea, O=ad) %>%
        select(sample, frag, chrom, pos, A, D, N, O)
    frags6 %>% filter(frag==frag_id) %>% 
        gather(pop, v, A:O) %>% filter(v>0) %>% 
        ggplot(aes(x=pos, y=v, fill=pop)) + 
        geom_col(width=w) + facet_wrap(~sample, ncol=1) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = .5, vjust=.5)) +
        scale_fill_manual(values=c( "#89bff7", "#ffb000","#444444","#efefef"))
}




save.image("salkhit2.rdebug")

