library(tidyverse)
library(cowplot)
library(reshape2)
source("scripts/plotting/lib.R")
#save.image("rdebug")
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
bin_size = as.numeric(snakemake@wildcards$bin_size) / 1e6


n_snps = 1349147
col_lomax = 'blue'
col_exp = 'red'
STATE = 1
TRUNC = .2
G_TIME = 29
SCALING = G_TIME * 100



ages = read_table2("config/ages.yaml", col_names=c("sample", "age"))

runs = load_runs_data(snakemake@input$runs, names) %>%
    left_join(ages) %>% 
    mutate(age=replace_na(age, 0)) %>%
    mutate(sample=fct_reorder(sample, -age))

# some plotting elements
NONE = rep(0, 4)
BLANK = element_blank()
STUFF = list(
    THEME , 
    facet_wrap(~sample, scale='free_y', ncol=1, strip='left'),
    theme(legend.position="none", 
          panel.spacing=unit(0, 'lines'),,
          strip.background=BLANK,
          strip.text=BLANK,
          axis.title.y=BLANK))



D_runs = runs %>%
    mutate(map_len = len * bin_size) 

#RUNSPLOT
runs_plot <- function(runs){
    P_runs = runs %>% 
	mutate(map_len = len * bin_size) %>%
	filter(state==STATE, map_len>=TRUNC, map_len<=7) %>% 
	group_by(sample, state) %>% 
	mutate(n=n/sum(n) / bin_size) %>% 
	ungroup %>%
	ggplot(aes(x=map_len, y=n)) + 
	geom_step() + 
	THEME +
	STUFF + 
	facet_wrap(~sample, ncol=1, strip='left', scale='free_y') + 
	scale_y_log10(expand=NONE, breaks=c(1e-3, 1))  +
	scale_x_continuous(expand=NONE, name='length (cM)')  +
	geom_text(aes(label=sprintf("%s: %s ky", sample, round(age/1000, 1))), y=Inf, x=Inf, hjust=1, vjust=1, data=runs %>% select(sample, age) %>% distinct, size=7)

}
P_runs = runs_plot(runs)

P1 = runs_plot(runs %>% filter(sample=="UstIshim")) + xlab("length (cM)")
ggsave("figures/paper/runs_UI.png", P1, width=8, height=4)
P2 = P1 + scale_y_continuous(expand=NONE)
ggsave("figures/paper/runs_UI2.png", P2, width=8, height=4)


pars = D_runs %>% filter(state==STATE) %>% 
    rle_fit_pars(trunc=TRUNC) %>% 
    mutate(delta_ll=delta_ll/200, p=pchisq(-2 * delta_ll, 1, lower=F)) %>% 
	left_join(ages) %>%
    mutate(age=replace_na(age, 0)) %>%
	mutate(scaled_age=age / SCALING) %>%
	mutate(semean = emean + scaled_age,
	       slmean = lmean + scaled_age,
           sample = fct_reorder(sample, age)
           ) 

fit_exp = pred_dexp(D_runs, pars, TRUNC) %>% bind_cols(pars %>% select(sample)) 
fit_lomax = pred_dlomax(D_runs, pars) %>% bind_cols(pars %>% select(sample)) 
fit_runs = inner_join(unnest(fit_exp), unnest(fit_lomax))

col_exp="red"
col_lomax="blue"
P_fit = P_runs + 
    geom_line(data=fit_runs %>% filter(exp>1e-3), aes(y=exp), color=col_exp, lty=2, lwd=1) + 
    #geom_line(data=fit_runs, aes(y=lomax), color=col_lomax, lty=2) + 
    coord_cartesian(xlim=c(TRUNC, 7), ylim=c(1e-3,10))  +
    scale_x_continuous(expand=NONE)




#P2 = P_prop + theme(axis.text.y=BLANK, axis.title.y=BLANK)


#P = plot_grid(P_prop, P_runs, P2, P_region, P_cov,
#          align = 'h', nrow=1,
#          axis='lb',
#          rel_widths = c(4, 4, 3, 4, 2)
#          )
P = plot_grid(P_runs, P_fit, 
          align = 'h', nrow=1,
          axis='lb',
          rel_widths = c(6, 6)
          )
ggsave(snakemake@output$plot, width=14, height=10)

