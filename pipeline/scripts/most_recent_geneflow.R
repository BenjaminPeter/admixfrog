library(dplyr)

n=1:100
G=37.04
p = 0.03

get_max_frag_len <- function(n, G, p, ...){
    x = 100 * rpois(1, n*G) %>%  #number of recombination events
        runif(0, G) %>% sort %>% #rec_pos
        diff %>% data.frame %>% sample_frac(p) #observed frags

    #    max, if any
    if(nrow(x) == 0) { return(NA)}
    return(max(x))
}

get_null_dist <- function(n=1:50, G=30, p=0.03, reps=10){
    v = sapply(n, function(n0)sapply(1:reps, get_max_frag_len, n=n0, G=30, p=.03))
    return(v)
}

v = get_null_dist(n=n, G=30, p=0.03, reps=10000) 
df = v %>% t %>% as_tibble %>% mutate(G=n) %>%
    pivot_longer(cols=-G) %>%
    select(-name) %>% 
    mutate(vbin=round(value, 1))

df %>% write_csv("most_recentdf.csv.xz")
