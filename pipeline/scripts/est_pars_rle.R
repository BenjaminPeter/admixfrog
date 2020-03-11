args = commandArgs(T)
print(args[1])
TARGET = "NEA"
TYPE = "state"
bin_size = 2000
N = 2000
TRUNC = 0.1

if(!is.na(args[1]) ) fname = as.character(args[1])
if(!is.na(args[2]) ) bin_size = as.numeric(args[2])
if(!is.na(args[3]) ) N = as.numeric(args[3])
if(!is.na(args[4]) ) TRUNC = as.numeric(args[4])

infile = sprintf("%s.rle.xz", fname)
outfile = sprintf("%s.est.csv", fname)

print("estimating parameters")
print(sprintf("infile:  %s", infile))
print(sprintf("outfile: %s", outfile))



suppressPackageStartupMessages({
require(tidyverse)
})


#a = read_csv("admixfrog/2000/AFR_VIN_DEN/Oase_archaicadmixture.res.xz")
a = read_csv(infile, progress=F) %>% 
    filter(target==TARGET, TYPE==type)

that <-function(r, N=10000, m=0.025)
    2 * N - 2 * sqrt (N* (r + N * (m-1)) / (m-1))

m = 0.04

#migration rate from posterior draws
m_post = a %>% 
    mutate(len = len * bin_size / 1e6) %>%
    filter(len >= TRUNC) %>%
    group_by(state, it) %>% 
    summarize(n=sum(len*n), 
              l=weighted.mean(len, weight=n) - TRUNC ) %>% 
    group_by(state) %>% 
    summarize(sdl = sd(l),
              sdm = sd(n),
              m=mean(n), 
              l=mean(l)
              ) %>%
    mutate(sdm=sdm/sum(m), m=m/sum(m) ) %>% 
    mutate(t=that(100/l, N, m), 
           tmin=that(100 / (l + 2* sdl ), N, (m - 2 * sdm)),
           tmax=that(100 / (l - 2* sdl ), N, (m + 2 * sdm)),
           T=100/l / (1-m),
           Tmin=100/ (l + 2 * sdl) / ( 1 - m - 2 * sdm),
           Tmax=100/ (l - 2 * sdl) / ( 1 - m + 2 * sdm))

print(m_post)

write_csv(m_post, outfile)


