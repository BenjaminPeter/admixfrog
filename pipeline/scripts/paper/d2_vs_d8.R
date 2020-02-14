#require(tidyverse)
#base2 = 'error2CAFR/5000/NEA_DEN/denisova2_hcneaden'
#base8 = 'error2CAFR/5000/NEA_DEN/denisova8_hcneaden'

#d2 = read_csv(sprintf("admixfrog/%s.bin.xz", base2), col_types=cols(chrom=col_character())) %>% mutate(intro2 = DEN < 0.5)
#d8 = read_csv(sprintf("admixfrog/%s.bin.xz", base8), col_types=cols(chrom=col_character())) %>% mutate(intro8 = DEN < 0.5)

#comp = d2 %>% select(id, intro2) %>% inner_join(select(d8, id, intro8), by='id')
#tbl = table(comp[,-1])
#ff = fisher.test(tbl)

require(tidyverse)
#base3 = 'error2CAFR/5000/A=ALT+D12_D=DEN+D11/denisova3_hcneaden'
#base5 = 'error2CAFR/5000/A=ALT+D12_D=DEN+D11/altai_hcneaden'

#d3 = read_csv(sprintf("admixfrog/%s.bin.xz", base3), col_types=cols(chrom=col_character())) %>% mutate(intro3 = AD > 0.5)
#d5 = read_csv(sprintf("admixfrog/%s.bin.xz", base5), col_types=cols(chrom=col_character())) %>% mutate(intro5 = AD > 0.5)

#comp = d3 %>% select(id, intro3) %>% inner_join(select(d5, id, intro5), by='id')
#tbl35 = table(comp[,-1])
#ff35 = fisher.test(tbl35)

base11 = 'error2CAFR/5000/NEA_DEN/denisova11_hcneaden'
#d11 = read_csv(sprintf("admixfrog/%s.bin.xz", base11), col_types=cols(chrom=col_character())) %>% mutate(intro11 = NEA > 0.5)

ft112 = d11 %>% select(id, chrom, intro11) %>% inner_join(select(d2, id, intro2), by='id') %>% 
    select(-id) %>% filter(chrom != 6) %>% select(-chrom) %>% table %>% fisher.test
ft118 = d11 %>% select(id, intro11) %>% inner_join(select(d8, id, intro8), by='id') %>% select(-id) %>% table %>% fisher.test
