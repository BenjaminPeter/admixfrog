library(tidyverse)

ldata = snakemake@input$res2 %>% lapply(read_csv)
names(ldata) = snakemake@config$panels[[snakemake@wildcards$panel]]
ldata %>% bind_rows(.id='sample') %>% write_csv(snakemake@output$res2_res)
save.image("merge.rdebug")
