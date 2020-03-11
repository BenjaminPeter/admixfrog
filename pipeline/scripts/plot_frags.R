source("scripts/plotting/lib.R")


frag = read_csv(snakemake@input$times)
#min_length = snakemake@params$min_length
min_length = as.numeric(snakemake@wildcards$trunc)

d1  = frag %>% filter(map_len > min_length) %>%
    mutate(sample=fct_reorder(sample, age)) %>%
    mutate(discovery_prob = 1 - exp(- (M-age) / G * min_length))



do_plot <- function(P){
    P + 
    geom_rect(aes(xmin=0, xmax=age, ymin=0, ymax=Inf), fill='lightgrey') + 
    geom_histogram(binwidth=1000, aes(y=stat(ncount))) + 
    facet_wrap(~sample, ncol=1, strip='left') +
    THEME + 
    scale_x_continuous(expand=rep(0,4))+ 
    scale_y_continuous(expand=rep(0,4)) +
    coord_cartesian(xlim=c(0, 80000), ylim=c(0, 40))
}

save.image("plot_frags.rdebug")

P_unw = ggplot(d1, aes(x=M)) %>% do_plot
P_mapw = ggplot(d1, aes(x=M, weight=map_len)) %>% do_plot
P_mdw = ggplot(d1, aes(x=M, weight=map_len / discovery_prob)) %>% do_plot

ggsave(snakemake@output$fragplot_unw, P_unw, width=5, height=18)
ggsave(snakemake@output$fragplot_mapw, P_mapw, width=5, height=18)
ggsave(snakemake@output$fragplot_mdw, P_mdw, width=5, height=18)
