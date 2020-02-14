suppressPackageStartupMessages({
#source("~/programs/admixfrog/plotting/comparison_plot.R")
library(tidyverse)
source("scripts/plotting/lib.R")
})


rle_files = snakemake@input$rle
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
n_samples = length(names)
target_ = snakemake@wildcards$target

trunc = as.numeric(snakemake@wildcards$trunc)
#trunc = snakemake@params$trunc
xmax = snakemake@params$xmax
gtime = snakemake@params$generation_time
type_ = snakemake@params$type
ages = read_table2(snakemake@input$ages, 
		   col_names=c('sample', 'age'))


data = load_rle_data(rle_files, names) %>%
	filter(target==target_, type==type_, map_len >= trunc)


SCALING = 100 * gtime #* 2
R = rle_fit_pars(data, trunc) %>%
	left_join(ages) %>%
	mutate(age=replace_na(age, 0)) %>%
	mutate(sage=age / SCALING) %>%
	mutate(semean = emean + age / SCALING, 
	       slmean = lmean + age / SCALING) %>%
        arrange(semean) %>% 
        mutate(sample=factor(sample, levels=unique(sample)))
write_csv(R, snakemake@output$fit)
save.image("rdebug8")

P2 =  plot_m_gamma(R, gtime) + ggtitle(snakemake@output$rleplot)
ggsave(snakemake@output$gammaplot, width=10, height=1 + (n_samples)  * .8, limitsize=F)

data = data %>% mutate(sample=factor(sample, levels=levels(R$sample)))
P = rle_fit_plot(data, R, trunc, xmax) + ggtitle(snakemake@output$rleplot)
ggsave(snakemake@output$rlelogplot, width=20, height=10)

P3= P + scale_y_continuous(name="") + ggtitle(snakemake@output$rleplot)
ggsave(snakemake@output$rleplot, width=20, height=10)


