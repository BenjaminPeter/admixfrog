CHROMS = c(1:22, "X", "Y")
GEN_TIME = 29
MIN_LENGTH = 0.2
library(GenomicRanges)
library(Homo.sapiens)
library(tidyverse)
library(scales)
source("scripts/plotting/fit.R")
annotate_genes <- function(frags){
    gfrags = frags %>% 
    mutate(chrom2=sprintf("chr%s", chrom)) %>%
    makeGRangesFromDataFrame(keep.extra.columns=T, 
                             seqnames='chrom2', start.f='pos', end.f='pos_end')
    all_genes  = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

    overlaps = findOverlaps(all_genes, gfrags) %>% as_tibble
    df_genes = all_genes %>% as_tibble %>% mutate(queryHits=1:nrow(.)) %>% 
        left_join(as.data.frame(org.Hs.egSYMBOL)) %>%
        select(queryHits, gene=symbol)
    df_frags = gfrags %>% as_tibble %>% mutate(subjectHits=1:nrow(.)) %>%
        select(-(1:5))
    df = overlaps %>% left_join(df_genes) %>% left_join(df_frags) %>% 
        group_by(chrom, target, type, map, map_end) %>% 
        summarize(genes=paste(gene, collapse=":")) %>% 
        right_join(frags) %>% 
        mutate(genes=replace_na(genes, ""))
} 

a_nea = read_csv(snakemake@input$a_nea)
a_den = read_csv(snakemake@input$a_den)
h_nea = read_csv(snakemake@input$h_nea)
h_den = read_csv(snakemake@input$h_den)

x_nea = read_csv(snakemake@input$x_nea) %>% mutate(target=ifelse(target=='ALT', 'A', target))
x_den = read_csv(snakemake@input$x_den) %>% mutate(target=ifelse(target=='DEN', 'D', target))

nice_names = rbind(
                   c('altai', 'Denisova 5'),
                   c('chagyrskaya08', 'Chagyrskaya 8'),
                   c('denisova11', 'Denisova 11'),
                   c('denisova2', 'Denisova 2'),
                   c('denisova3', 'Denisova 3'),
                   c('denisova4', 'Denisova 4'),
                   c('denisova8', 'Denisova 8'),
                   c('goyet', 'Goyet Q56-1'),
                   c('hst', 'HST'),
                   c('lescottes', 'Les Cottés Z4-1514'),
                   c('mez1', 'Mezmaiskaya 1'),
                   c('mez2', 'Mezmaiskaya 2'),
                   c('scladina', 'Scladina'),
                   c('spy1', 'Spy 1'),
                   c('vindija3319', 'Vindija 33.19')
                   ) %>% as_tibble
names(nice_names) = c('sample', 'Sample')

targets = rbind(
                c('altai', 'DEN', 'state'),
                c('chagyrskaya08', 'DEN', 'state'),
                c('mez1', 'DEN', 'state'),
                c('mez2', 'DEN', 'state'),
                c('spy1', 'DEN', 'state'),
                c('goyet', 'DEN', 'state'),
                c('lescottes', 'DEN', 'state'),
                c('vindija3319', 'DEN', 'state'),
                c('scladina', 'DEN', 'state'),
                c('hst', 'DEN', 'state'),
                c('denisova2', 'NEA', 'state'),
                c('denisova4', 'NEA', 'state'),
                c('denisova8', 'NEA', 'state'),
                c('denisova11', 'NEA', 'homo'),
                c('denisova11', 'DEN', 'homo'),
                c('denisova3', 'A', 'state'),
                c('altai', 'D', 'state')
                ) %>%
    apply(1, paste, collapse=':')

df0 = bind_rows(a_nea, a_den, h_nea, h_den, x_den, x_nea) 
df = df0 %>%
    mutate(chrom=factor(chrom, levels=CHROMS)) %>%
    arrange(chrom, pos, map) %>%
    mutate(sample = as.factor(sample)) %>%
    mutate(target = as.factor(target)) %>%
    mutate(type = as.factor(type)) %>%
    filter(too_big == 'NO', map_len >= MIN_LENGTH) %>% 
    mutate(sex=chrom=='X') %>%
    select(sample, target, type, chrom, pos, pos_end, map, map_end, id=start, id_end=end, map_len, pos_len, sex, len)

tbl = df %>% group_by(sample, target, type, .drop=F) %>% 
    summarize(n_auto = sum(!sex), 
              n_sex=sum(sex),
              tot_cm = sum(map_len),
              tot_pos = sum(pos_len) / 1e6,
              longest_cm = max(map_len),
              longest_pos = max(pos_len) / 1e6,
              )

frags = df %>% rowwise %>% 
    mutate(x=paste(c(as.character(sample), as.character(target), as.character(type)),collapse=':')) %>%
    filter(x %in% targets)  %>%
    select(-x, -sex) %>%
    arrange(sample, target, type, chrom, pos) %>% annotate_genes %>%
    select(sample, target, type, chrom, pos, pos_end, map, map_end, id, id_end, pos_len, map_len, len, genes) 


frags %>% write_csv(snakemake@output$tableS2)


tbl2 = tbl %>% rowwise %>% 
    mutate(x=paste(c(as.character(sample), as.character(target), as.character(type)),collapse=':')) %>%
    filter(x %in% targets) %>%
    select(-x) %>% arrange(sample, -tot_cm) %>%
    mutate(tot_pos = round(tot_pos, 2))

tbl3 = frags %>% group_by(sample, target, type) %>% 
    rle_fit_grp(trunc=MIN_LENGTH) %>% mutate(gen=emean * 100, years=gen * 29) %>% 
    arrange(years) %>% 
    select(sample, years, gen, type, target) %>% left_join(tbl2, .) %>%
    left_join(read_tsv(snakemake@input$ages)) %>%
    mutate(se_years = years / sqrt(n_auto + n_sex), se_gen = gen / sqrt(n_auto + n_sex)) %>%
    right_join(nice_names, .) %>% select(-sample)  %>%
    mutate(Years = sprintf("%s±%s", 100 * round(years / 100), 100 * round(2*se_years/ 100))) %>% 
    mutate(Years=ifelse(n_auto + n_sex >=10, Years, '-')) %>%
    mutate(Gen = sprintf("%s±%s", round(gen), round(2*se_gen))) %>% 
    mutate(Gen=ifelse(n_auto + n_sex >=10, Gen, '-')) %>%
    mutate(l1=ifelse(longest_cm>0, round(longest_cm, 2), '-'), 
           l2=ifelse(longest_pos>0, round(longest_pos, 2), '-')) %>%
    mutate(t1=round(tot_cm, 1),
           t2=round(tot_pos, 1),
           n=sprintf("%s/%s", n_auto, n_sex)) %>%
    rename('Max(cM)' = l1, 'Max(Mb)' = l2,
           'Tot(cM)' = t1, 'Tot(Mb)' = t2,
           'n(A/X)' = n)
levels(tbl3$type) = c('(h)', '', '*')
tbl3 = tbl3 %>% 
    mutate(type = ifelse(str_length(target)==1, '*', as.character(type))) %>%
    mutate(target=as.character(target)) %>%
    mutate(target=ifelse(target == 'A', 'NEA', target)) %>%
    mutate(target=ifelse(target == 'D', 'DEN', target)) %>%
    mutate(type = sprintf("%s%s", target, type)) %>%
    rename(Type=type, 'Age Range'='age_range')

tbl3_order = c(6,9,1,2,7,4,5,3,13,14,10,17,11,12,16,15,8)
tbl3[tbl3_order,] %>%  select(1, 3, 12, 22, 20, 21, 18, 19, 17, 16, 13) %>% 
    write_csv(snakemake@output$tableS1)
 


