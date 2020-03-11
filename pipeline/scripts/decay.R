require(tidyverse)
#infile="stats/frags2panel/gtmode_error/state_NEA_0.2/AFR_NEA/5000/simons_paper5c_archaicadmixture.frags"
#outplot="stats/frags2panel/gtmode_error/state_NEA_0.2/AFR_NEA/5000/simons_paper5c_archaicadmixture.png"
#outdata="stats/frags2panel/gtmode_error/state_NEA_0.2/AFR_NEA/5000/simons_paper5c_archaicadmixture.csv.gz"
infile = snakemake@input$frags
outplot = snakemake@output$plot
outdata = snakemake@output$plot_data
sgdp = read_csv(snakemake@input$meta)

a = read_csv(infile) %>% 
    select(-type) %>% 
    left_join(sgdp) %>%
    mutate(id=ifelse(type=='modern', pop, sample))


data = a %>% filter(map_len > 0.2) %>% 
    group_by(id, region, age) %>% 
    summarize(n_samples=n_distinct(sample),
              len=mean(map_len), 
              n_tracts=n()/n_samples,
              tot=sum(map_len)/n_samples)


P = data %>% filter(region != "Africa") %>% 
    ggplot(aes(x=age, y=1/len, label=id, color=region)) + 
    geom_text()

P %>% ggsave(outplot, ., width=8, height=6)


data %>% write_csv(outdata)

