library(tidyverse)
cont = snakemake@input$cont
outfile = snakemake@output$cont
out2 = snakemake@output$cov
names = snakemake@config$panels[[snakemake@wildcards$panel]]

n_sites = snakemake@input$n_sites[1] %>% read_csv %>% nrow

cont_files = lapply(cont, read_csv)
names(cont_files) = names

save.image("test.rdebug")

cont_files %>% bind_rows(.id="sample") %>% write_csv(outfile) %>%
    group_by(sample) %>% 
    summarize(error=weighted.mean(error, n_snps), 
              cont=weighted.mean(cont, n_snps), 
              n_reads=sum(n_snps), 
              cov=n_reads/n_sites) %>% 
    write_csv(out2)



