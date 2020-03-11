require(rstan)
require(tidyverse)

#options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
later:::ensureInitialized()
K = snakemake@params$K
trunc = snakemake@params$trunc

save.image("runstanfrag.rdebug")

#frags = readRDS("frags.rds") %>%
frags = read_csv(snakemake@input$frags) %>%
    filter(map_len >=trunc) %>% 
    mutate(sample_id = as.integer(as.factor(sample)),
        fid= as.integer(interaction(chrom, frag, sep='|', drop=T))) 

stan_data0 = frags %>% 
    select(sample_id, map_len, fid, age) %>%
    filter(!is.na(sample_id+map_len + fid + age))

stan_data = as.list(stan_data0)
stan_data$trunc_ = trunc
stan_data$N = length(stan_data$age)
stan_data$t = exp(seq(log(35000), log(120000), length.out=K))
stan_data$K = K
stan_data$F = length(unique(stan_data0$fid))


#fit model
model = stan_model(file='frag.stan')
o = optimizing(model, stan_data, hessian=F, verbose=T, tol_rel_obj=1,
               init=list(theta=as.list(rep(1/K), K),
                         gtime = 4500,
                         t = as.list(seq(35000, 100000, length.out=K))))
theta = o$par[1:K]
#times = o$par[(K+2):(2*K+1)]
gtime = o$par[K+1]
pis = o$par[-1:-(1*K+1)]


pi_mat = exp(matrix(pis, ncol=K)) ; 
pi_mat[pi_mat==Inf] = 1e10; 
pi_mat = pi_mat / rowSums(pi_mat); 
colnames(pi_mat) = sprintf("pi%s", 1:K); 

M = colSums(t(pi_mat) * stan_data$t)
res = stan_data0 %>% group_by(fid) %>% tally %>% arrange(fid) %>% bind_cols(tibble(M=M))  %>% left_join(frags)
res %>% write_csv(snakemake@output$res)

