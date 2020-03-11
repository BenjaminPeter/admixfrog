require(rstan)
require(tidyverse)

#options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
later:::ensureInitialized()


load("fig1.rdebug")

trunc = 0.05
K = 10

bin_size = as.numeric(snakemake@wildcards$bin_size) * 1e-6
rle_files = snakemake@input$runs
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]



#pre-processing
runs_data = load_runs_data(rle_files, names) %>% 
    ungroup %>%
    mutate(map_len= len * bin_size) %>%
    filter(state==1)  %>%
    select(sample, map_len, n)
ages = read_table2("config/ages.yaml",
           col_names=c('sample', 'age'))
stan_data0 = runs_data %>% left_join(ages) %>%
    mutate(age=replace_na(age, 0)) %>%
    filter(map_len >=trunc) 

stan_data = as.list(stan_data0)
stan_data$trunc_ = trunc
stan_data$N = length(stan_data$n)
stan_data$t = exp(seq(log(35000), log(120000), length.out=K))
stan_data$K = K


#fit model
model = stan_model(file='test.stan')
o = optimizing(model, stan_data, hessian=F, verbose=T, tol_rel_obj=1,
               init=list(theta=as.list(rep(1/K), K),
                         gtime = 2600,
                         t = as.list(seq(35000, 120000, length.out=K))))
theta = o$par[1:K]
times = o$par[(K+2):(2*K+1)]
gtime = o$par[K+1]
pis = o$par[-1:-(2*K+1)]

#  post-processing
pars = data.frame(theta, times) %>% arrange(times)

pi_mat = exp(matrix(pis, ncol=K)) ; 
pi_mat[pi_mat==Inf] = 1e10; 
pi_mat = pi_mat / rowSums(pi_mat); 
colnames(pi_mat) = sprintf("pi%s", 1:K); 
G = cbind(stan_data0, pi_mat) 
H = G %>% group_by(sample, age) %>% summarize_at(vars(starts_with("pi")), weighted.mean, weight= .$n )
Z = H %>% gather(k,v, -sample, -age) %>% left_join(data.frame(k=sprintf("pi%s",1:K), t=times))
Z  %>% ggplot(aes(x=t, y=v, group=sample, color=sample)) + geom_point()  + geom_line(direction='vh')



