L2 = L %>% filter(state=="NEA", len>=.20)  %>% mutate(sample=fct_reorder(sample, age))
D = L2 %>% 
    rename(map_len=len) %>% 
    group_by(sample, map_len,bin_size, it) %>% 
    tally  %>% 
    group_by(sample, bin_size, it) 
Z = D %>% rle_fit_grp(0.2) %>% left_join(meta) %>%
    mutate(sample=fct_reorder(sample, age))


V2 = V %>% filter(state=="NEA", map_len>=.20)  %>% mutate(sample=fct_reorder(sample, age))
DV2 = V2 %>% 
    group_by(sample, map_len,bin_size) %>% 
    tally  %>% 
    group_by(sample, bin_size) 
ZV2 = DV2 %>% rle_fit_grp(0.2) %>% left_join(meta) %>%
    mutate(sample=fct_reorder(sample, age))
