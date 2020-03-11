require(rstan)
require(tidyverse)

stan_script = snakemake@input$stan
n_chains = snakemake@params$n_chains
n_iter = snakemake@params$n_iter
options(mc.cores = n_chains)
rstan_options(auto_write = TRUE)
later:::ensureInitialized()
K = snakemake@params$K
#trunc = snakemake@params$trunc
trunc = as.numeric(snakemake@wildcards$trunc)

#save.image("runstanfrag.rdebug")

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

stan_data$tmin = stan_data0 %>% group_by(fid) %>% 
    arrange(fid) %>%
    summarize(max_age=max(age)) %>%
    mutate(tmin = sapply(max_age, function(a)min(which(stan_data$t > a))))  %>%
    select(tmin) %>% unlist



#fit model
#frag_fit = stan("frag.stan", verbose=T, data=stan_data, 
frag_fit = stan(stan_script, verbose=T, data=stan_data, 
                chains=n_chains, iter=n_iter)

e = rstan::extract(frag_fit)

pmean = apply(e$pi, 2:3, mean)
M = colSums(t(pmean) * stan_data$t)
sdM = apply(t(pmean) * stan_data$t, 2, sd)
G = colMeans(e$gtime)
sdG = apply(e$gtime, 2, sd)

times = as.data.frame(t(stan_data$t)) %>% t %>% as.data.frame %>% 
    rownames_to_column(var='time_bin')

res1 = stan_data0 %>% group_by(fid) %>% tally %>% arrange(fid) %>% 
    bind_cols(tibble(
                     M=M,
                     sdM=sdM,
                     G=G,
                     sdG=sdG,
                     ))  %>% left_join(frags)
res1 %>% write_csv(snakemake@output$res)

res2 = stan_data0 %>% group_by(fid) %>% tally %>% arrange(fid) %>% 
   bind_cols(as.data.frame(pmean)) %>% 
   left_join(frags) %>%
   gather(time_bin, p, sprintf("V%s", 1:K)) %>% left_join(times) %>%
   rename(time=V1)
res2 %>% write_csv(snakemake@output$res2)
