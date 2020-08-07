require(tidyverse)
base2 = 'error2CAFR/5000/NEA_DEN/denisova2_hcneaden'
base8 = 'error2CAFR/5000/NEA_DEN/denisova8_hcneaden'
base11 = 'error2CAFR/5000/NEA_DEN/denisova11_hcneaden'
MHC=c(28477797, 33448354)

filter_mhc <- function(x){
    x %>% filter(! (chrom=='6' & pos > MHC[1] & pos < MHC[2]))
}
filter_x <- function(x)
    x %>% filter(chrom !='X')

d2 = read_csv(sprintf("admixfrog/%s.bin.xz", base2), col_types=cols(chrom=col_character())) %>% mutate(intro2 = DEN < 0.5) %>% filter_mhc 
d8 = read_csv(sprintf("admixfrog/%s.bin.xz", base8), col_types=cols(chrom=col_character())) %>% mutate(intro8 = DEN < 0.5) %>% filter_mhc
d11 = read_csv(sprintf("admixfrog/%s.bin.xz", base11), col_types=cols(chrom=col_character())) %>% mutate(intro11 = NEA > 0.5) %>% filter_mhc

comp28 = d2 %>% select(id, intro2) %>% inner_join(select(d8, id, intro8), by='id')
tbl28 = table(comp28[,-1])
ff28 = fisher.test(tbl28)

comp211 = d2 %>% select(id, intro2) %>% inner_join(select(d11, id, intro11), by='id')
tbl211 = table(comp211[,-1])
ff211 = fisher.test(tbl211)

comp811 = d8 %>% select(id, intro8) %>% inner_join(select(d11, id, intro11), by='id')
tbl811 = table(comp811[,-1])
ff811 = fisher.test(tbl811)
