suppressPackageStartupMessages({
require(tidyverse)
require(VGAM)
})

bin_size = 2000
N = 5000
m = 0.03

#a = read_csv("admixfrog/2000/AFR_VIN_DEN/Oase_archaicadmixture.res.xz")
#a = read_csv("admixfrog/2000/AFR_VIN_DEN/Kostenki14sg_archaicadmixture.rle.xz")
#a = read_csv("admixfrog/2000/AFR_VIN_DEN/LBK_archaicadmixture.rle.xz") %>%
a = read_csv("admixfrog/delta0.2/2000/AFR_NEA_DEN/French_archaicadmixture.rle.xz") %>%
    filter(target=='NEA', type=='state') %>%
    select(map_len) %>% arrange(-map_len)
           
b = read_csv("admixfrog/delta0.2/2000/AFR_NEA_DEN/French_archaicadmixture.bin.xz") 
b1 = read_csv("admixfrog/delta0.1/2000/AFR_NEA_DEN/French_archaicadmixture.bin.xz") 
b2 = read_csv("admixfrog/basic/2000/AFR_NEA_DEN/French_archaicadmixture.bin.xz") 
B = read_csv("admixfrog/delta0.2/2000/AFR_NEA_DEN/UstIshim_archaicadmixture.bin.xz") 
O = read_csv("admixfrog/delta0.2/2000/AFR_NEA_DEN/Oase_archaicadmixture.bin.xz") 
K = read_csv("admixfrog/delta0.2/2000/AFR_NEA_DEN/Kostenki14sg_archaicadmixture.bin.xz") 
k = read_csv("admixfrog/basic/10000/AFR_NEA_DEN/French_archaicadmixture.res.xz") %>% filter(state==1)
s = read_csv("admixfrog/delta0.2/2000/AFR_NEA_DEN/Sardinian_archaicadmixture.res.xz")


fit <- function(...){
    fl = fit_lomax(...)
    fe = fit_exp(...)
    x= (cbind(fl, fe) %>% mutate(delta_ll = ll_lomax-ll_exp))
    print(x)
    return(x)
}

fit_lomax = function(len, n=1, lmax=200, trunc=1e-1, init=c(-4, -4)){
    if(length(n)==1)n = rep(n, length(len))
    n = n[len>=trunc & len < lmax]
    len = len[len>=trunc & len < lmax]

    if(length(len) == 0) return(data.frame(ll_lomax=NA, scale=NA, shape=NA))
    f <- function(pars){
        ll = sum(dlomax(len / 100, exp(pars[1]), exp(pars[2]), log=T) * n)
        den = (plomax(lmax / 100, exp(pars[1]), exp(pars[2])) - 
               plomax(trunc/ 100, exp(pars[1]), exp(pars[2])))
        ll - log(den) * sum(n)
    }
    o = optim(init, f, method='BFGS', control=list(fnscale=-1))
    ll = o$value
    pars = exp(o$par)
    return(data.frame(ll_lomax=ll, lrate=pars[1], shape=pars[2], lmean = pars[2] / pars[1]))
}
fit_exp = function(len, n=1, lmax=1, trunc=.05, ...){
    if(length(n)==1)n = rep(n, length(len))
    n = n[len>=trunc & len < lmax]
    len = len[len>=trunc & len < lmax]
    if(length(len) == 0) return(data.frame(ll_lomax=NA, scale=NA, shape=NA))
    f <- function(pars){
        ll = sum(dexp(len /100, pars[1], log=T) * n)
        den = (pexp(lmax /100 , pars[1]) - pexp(trunc /100, pars[1]))
        return(ll - log(den) * sum(n))
    }
    o = optim(c(1), f, method='L-BFGS-B', lower=1e-6, upper=3e5, control=list(fnscale=-1))
    ll = o$value
    pars = o$par
    return(data.frame(ll_exp=ll, erate=pars[1]))
}
lx = a %>% 
    filter(it<20) %>%
    group_by(it, state) %>% 
    mutate(len = len * bin_size / 1e6) %>% 
    do( fit(.$len, .$n, lmax=100, trunc=2e-2)) 
#    do( fit_exp(.$len, .$n, lmax=200))
