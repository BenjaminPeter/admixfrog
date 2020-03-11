#I want full-page figure for ~20 samples:
# - 1 row is ~1cM
# panels are
#     1. prop archaic (2cM)
#     2. prop cont (2cM)
#     3. run-length-distribution (4 cm)
#     4. est. age (3 cm)
#     5. sample fragment (4 cm)



library(tidyverse)
library(cowplot)
library(reshape2)
source("scripts/plotting/lib.R")
#save.image("fig1.rdebug")
infiles = snakemake@input$bins
bin_size = as.numeric(snakemake@wildcards$bin_size) / 1e6
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]

TARGET = snakemake@wildcards$target
STATE = 1


n_snps = 1349147
col_lomax = 'blue'
col_exp = 'red'
TRUNC = .2

#G_TIME = 29
#SCALING = G_TIME * 100



#load data
cont = load_cont_data(snakemake@input$cont, names)
ages = read_table2("config/ages.yaml", col_names=c("sample", "age"))
rle = load_rle_data(snakemake@input$rle, names) %>%
    filter(target==TARGET, type=='state') %>%
    left_join(ages) %>% 
    mutate(age=replace_na(age, 0)) %>%
    mutate(sample=fct_reorder(sample, -age))
runs = load_runs_data(snakemake@input$runs, names) %>%
    left_join(ages) %>% 
    mutate(age=replace_na(age, 0)) %>%
    mutate(sample=fct_reorder(sample, -age))
data = load_bin_data(infiles, names, widths=F) %>%
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




# PROP plot
D_prop = data %>% group_by(sample) %>%
    select(-sample:-n_snps, -age) %>%
    summarize_all(mean)  %>% melt %>%
    rename(state=variable, prob=value)  %>%
    mutate(state=fct_reorder(state, -prob))


P_prop <- D_prop %>%
	ggplot(aes(x=sample, y=prob, fill=state)) +
	geom_col() + 
	coord_flip(ylim=c(0, .05)) +
    col_scale() +
    STUFF + 
    scale_x_discrete(expand=NONE) +
    scale_y_continuous(expand=NONE, name="proportion") +
    theme(legend.position="none", 
          panel.spacing=unit(0, 'lines'),,
          strip.background=BLANK,
          strip.text=BLANK,
          axis.title.y=BLANK)



#CONTAMINATION
D_cont = cont %>% 
    left_join(ages) %>% 
    mutate(sample=fct_reorder(sample, -age)) %>% 
    mutate(deam=ifelse(endsWith(lib, 'nodeam'), 'nodeam', 'deam')) %>%
    mutate(deam=replace_na(deam, 'nodeam')) %>%
    group_by(sample, deam) %>% 
    summarize(p_cont=weighted.mean(cont, weight=n_snps),
              cov=sum(n_snps))  %>%
    mutate(cont=p_cont*cov, endogenous=cov - cont) %>%
    ungroup %>%
    complete(sample, deam, fill=list(cont=0))

P_cov = D_cont %>% ggplot(aes(x=deam, y= 100 *cov/n_snps, fill=deam)) +
    geom_col() + facet_wrap(~sample, ncol=1, strip='left') +
    coord_flip(ylim=c(1,5000)) + 
    scale_y_log10(name="Cov", expand=NONE, breaks=c(1, 10, 100, 1000), labels=c(0.01, 0.1, 1, 10)) +
    THEME + STUFF +
    geom_text(aes(y=3, x = deam, label=deam), hjust=0, size=2) +
    geom_hline(yintercept=100, col='grey', lty=2, lwd=.5, alpha=0.5) +
    theme(axis.text.y=BLANK, axis.title.y=BLANK)

P_cov2 = D_cont %>% 
    select(-cont) %>% rename(cont=cov) %>% gather(k, v, cont:endogenous) %>%
    ggplot(aes(x=deam, y= 100 * v/n_snps, fill=k)) +
    geom_col(position='identity') + 
    facet_wrap(~sample, ncol=1, strip='left') +
    coord_flip(ylim=c(1,5000)) + 
    scale_y_log10(name="Cov", expand=NONE, breaks=c(1, 10, 100, 1000), labels=c(0.01, 0.1, 1, 10)) +
    THEME + STUFF +
    geom_text(aes(y=3, x = deam, label=deam), hjust=0, size=2) +
    geom_hline(yintercept=100, col='grey', lty=2, lwd=.5, alpha=0.5) +
    theme(axis.text.y=BLANK, axis.title.y=BLANK)


REGION_CHROM = 9
REGION_START = 94
REGION_END = 99

D_region = data %>%
    filter(chrom==REGION_CHROM, REGION_START < map , map < REGION_END)  %>%
    group_by(sample, chrom) %>%
    mutate(pos_end = lag(pos)) %>%
    ungroup %>%
    gather(variable, value, -sample:-n_snps, -age, -pos_end) %>%
    filter(variable != "AFR", value > 0.01)

P_region = D_region %>% 
    ggplot(aes(x=map, y=value, fill=variable)) +
    geom_col(width=bin_size) + 
    facet_wrap(~sample, ncol=1, strip='l') +
    STUFF + 
    col_scale() +
    scale_x_continuous(expand=NONE, name='chr9 (cM)' ) +
    scale_y_continuous(breaks=NULL, expand=NONE) +
    theme(legend.position="none", 
          panel.spacing=unit(0, 'lines'),,
          strip.background=BLANK,
          strip.text=BLANK,
          axis.title.y=BLANK)

P_region2 = D_region %>% 
    group_by(chrom, map, pos, sample) %>% mutate(value=cumsum(value)) %>%
    ggplot(aes(xmin=pos / 1e6, xmax=pos_end/1e6, ymin=0, ymax=value, fill=variable)) +
    geom_rect() +
    facet_wrap(~sample, ncol=1, strip='l') +
    STUFF + 
    col_scale() +
    scale_x_continuous(expand=NONE, name='chr9 (Mb)' ) +
    scale_y_continuous(breaks=NULL, expand=NONE) +
    theme(legend.position="none", 
          panel.spacing=unit(0, 'lines'),,
          strip.background=BLANK,
          strip.text=BLANK,
          axis.title.y=BLANK) +
    coord_cartesian(xlim=c(94, 99.2))


D_runs = runs %>%
    mutate(map_len = len * bin_size) 

#RUNSPLOT
D_runs2 = runs %>% 
    mutate(map_len = len * bin_size) %>%
    filter(state==STATE, map_len>=TRUNC, map_len<=4) %>% 
    group_by(sample, state) %>% 
    mutate(n=n/sum(n) / bin_size) 
rmin = round(min(D_runs2$n), 3)
P_runs = D_runs2 %>%
    ggplot(aes(x=map_len, y=n)) + 
    geom_step() + 
    THEME +
    STUFF + 
    theme(legend.position="none", 
          panel.spacing=unit(0, 'lines'),,
          strip.background=BLANK,
          strip.text=BLANK,
          axis.title.y=BLANK) +
    facet_wrap(~sample, ncol=1, strip='left', scale='free_y') + 
    scale_y_log10(expand=NONE, breaks=c(round(rmin, 3), 1)) 

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
    geom_line(data=fit_runs %>% filter(exp>rmin), aes(y=exp), color=col_exp, lty=2, lwd=.5) + 
    #geom_line(data=fit_runs, aes(y=lomax), color=col_lomax, lty=2) + 
    coord_cartesian(xlim=c(TRUNC, 4), ylim=c(rmin,10))  +
    scale_x_continuous(expand=NONE, name="fragment length (cM)")


D_frags = read_csv(snakemake@input$frags) %>%
    filter(map_len >= TRUNC) %>%
    mutate(sample=fct_reorder(sample, -age)) %>%
    mutate(discovery_prob = exp(- (M-age) / G * TRUNC))

P_frags = D_frags %>%
    ggplot(aes(x=M/1000, weight=map_len / discovery_prob)) +
    geom_rect(aes(xmin=0, xmax=age/1000, ymin=0, ymax=Inf), fill='lightgrey') + 
    geom_histogram(binwidth=1, aes(y=stat(ncount))) + 
    facet_wrap(~sample, ncol=1, strip='left') +
    THEME + STUFF +
    scale_x_continuous(expand=rep(0,4), name="time (kya)", breaks=seq(0, 80, 10))+ 
    scale_y_continuous(expand=rep(0,4), breaks=c(0.5, 1)) +
    coord_cartesian(xlim=c(0, 80))  +
    theme(legend.position="none", 
          panel.spacing=unit(0, 'lines'),
          strip.background=BLANK,
          strip.text=BLANK,
          axis.title.y=BLANK)






P2 = P_prop + theme(axis.text.y=BLANK, axis.title.y=BLANK)


P = plot_grid(P_prop, P_fit, P_frags, P_region2, P_cov2,
          align = 'h', nrow=1,
          axis='lb',
          labels=LETTERS[1:5],
          rel_widths = c(2, 3.5, 3.5,  3, 1.5)
          )
#P = plot_grid(P_prop, P_fit, P_region, 
#          align = 'h', nrow=1,
#          axis='lb',
#          rel_widths = c(6, 6, 4)
#         )
ggsave(snakemake@output$plot, width=8, height=8)

