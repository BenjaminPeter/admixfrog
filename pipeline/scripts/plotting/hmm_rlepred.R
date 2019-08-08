#' main prediction script

suppressPackageStartupMessages({
library(tidyverse)
source("scripts/plotting/lib.R")
})

save.image("rlepred.rdebug")



rle_files = snakemake@input$rle
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
n_samples = length(names)
target_ = snakemake@wildcards$target
trunc = as.numeric(snakemake@wildcards$trunc)
xmax = snakemake@params$xmax
gtime = snakemake@params$generation_time
ages = read_table2(snakemake@input$ages, 
		   col_names=c('sample', 'age'))
states=strsplit(snakemake@wildcards$states, "_")[[1]]
bin_size = as.numeric(snakemake@wildcards$bin_size) * 1e-6

gtime = 51.4
SCALING = 100 * gtime



runs_data = load_runs_data(rle_files, names) %>%
	ungroup() %>%
	mutate(map_len=len*bin_size) %>%
	filter(state==target_) %>% 
	select(sample, map_len, n)


R = runs_data %>%
	rle_fit_pars(trunc) %>%
	left_join(ages) %>%
	mutate(age=replace_na(age, 0)) %>%
#	mutate(scaled_age=age / SCALING) %>%
#	mutate(semean = emean + scaled_age,
#	       slmean = lmean + scaled_age) %>%
#        arrange(semean) %>% 
        mutate(sample=factor(sample, levels=unique(sample)),
               delta_ll = delta_ll /200)

write_csv(R, snakemake@output$fit)

data = runs_data %>% 
    group_by(sample) %>% 
	do(map_len=rep(.$map_len, .$n)) %>% 
	unnest


P2 =  plot_m_gamma(R, gtime) + ggtitle(snakemake@output$rleplot)
ggsave(snakemake@output$gammaplot, width=10, height=1 + (n_samples)  * .8, limitsize=F)

data = data %>% mutate(sample=factor(sample, levels=levels(R$sample)))
P = rle_fit_plot(data, R, trunc, xmax) + ggtitle(snakemake@output$rleplot)
ggsave(snakemake@output$rlelogplot, width=20, height=10)

P3= P + scale_y_continuous(name="") + ggtitle(snakemake@output$rleplot)
ggsave(snakemake@output$rleplot, width=20, height=10)


