require(tidyverse)
a = read_csv("stats/frags2panel/gtmode_error/state_NEA_0.1/AFR_NEA_DEN/10000/simons_paper5b_archaicadmixture.frags")

ff = a %>% filter(!is.na(frag), age>0) %>% group_by(sample) %>% tally %>% filter(n>100)

X = ff %>% left_join(a) %>% group_by(sample, frag) %>% tally %>% mutate(n=n>0) %>% spread(sample, n, fill=0) %>%
    select(-frag)
 D = dist(t(X), 'binary') %>% as.matrix
 diag(D) <- NA
