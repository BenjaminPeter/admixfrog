require(tidyverse)
save.image("cont.rdebug")
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
data = lapply(snakemake@input$samples, read_csv) %>% 
	bind_rows(.id="sample") %>%
	mutate(sample=names[as.integer(sample)])

data %>% write_csv(snakemake@output[[1]])

