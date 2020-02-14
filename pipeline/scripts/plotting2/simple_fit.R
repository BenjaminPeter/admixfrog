source("~benjamin_peter/pipeline/scripts/plotting/fit.R")

MIN_LEN = 0.1

d1 = read_csv("admixfrog/error/5000/AFR_NEA_DEN/Salkhit_archaicadmixture.rle.xz") %>%
    mutate(sample="Salkhit") %>%
    filter(map_len >= MIN_LEN, type=="state", target!="AFR") 
d2 = read_csv("admixfrog/error/5000/AFR_NEA_DEN/Tianyuan_archaicadmixture.rle.xz") %>%
    mutate(sample="Tianyuan") %>%
    filter(map_len >= MIN_LEN, type=="state", target!="AFR") 
d3 = read_csv("admixfrog/error/5000/AFR_NEA_DEN/Yana1_archaicadmixture.rle.xz") %>%
    mutate(sample="Yana1") %>%
    filter(map_len >= MIN_LEN, type=="state", target!="AFR") 
d4 = read_csv("admixfrog/error/5000/AFR_NEA_DEN/Yana2_archaicadmixture.rle.xz") %>%
    mutate(sample="Yana2") %>%
    filter(map_len >= MIN_LEN, type=="state", target!="AFR") 
d5 = read_csv("admixfrog/error/5000/AFR_NEA_DEN/Malta_archaicadmixture.rle.xz") %>%
    mutate(sample="Malta") %>%
    filter(map_len >= MIN_LEN, type=="state", target!="AFR") 

data = bind_rows(d1, d2, d3, d4, d5)

pars = data %>%
    group_by(sample, target, map_len)  %>%
    tally %>%
    group_by(sample, target) %>%
    rle_fit_grp(MIN_LEN)



P = data %>% ggplot(aes(color=target, x=map_len)) + 
    stat_ecdf() + 
    scale_x_reverse() + 
    scale_y_log10(lim=c(0.10, 1)) + 
    facet_grid(sample~., scales='free', space='free')
    

