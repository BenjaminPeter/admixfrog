source("~/programs/admixfrog/plotting/comparison_plot.R")

#save.image(".rdebug")
infiles = snakemake@input$bins
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
data = load_data(infiles, names)

posterior = data %>% group_by(sample) %>% 
    select(-1:-9) %>%
    summarize_all(mean)  %>% melt %>%
    rename(state=variable, prob=value)

max = data %>% group_by(sample, viterbi) %>% 
    tally() %>%
    mutate(n=n/sum(n)) %>% 
    rename(state=viterbi, max=n)


df = full_join(posterior, max) %>% 
    melt %>%
    filter(state != "GARBAGE") %>%
    group_by(state) %>% 
    arrange(-value) %>%
    ungroup %>% 
    mutate(state=factor(state, levels=unique(state)))

P = ggplot(df, aes(x=variable, y=value, fill=state)) + 
    geom_col() + 
    facet_grid(~sample)  + 
    coord_cartesian(ylim=c(0, 1))  + 
    theme(legend.position="bottom")



n_samples = length(names)
ggsave(snakemake@output$plot, width=.8 * n_samples, height=4, limitsize=F)
P = ggplot(df, aes(x=variable, y=value, fill=state)) + 
    geom_col() + 
    facet_grid(~sample)  + 
    coord_cartesian(ylim=c(0, .15))  + 
    theme(legend.position="bottom")
ggsave(snakemake@output$plot2, width=.8 * n_samples, height=4, limitsize=F)


