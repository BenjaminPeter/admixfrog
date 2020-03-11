require(tidyverse)
lapply(snakemake@input, read_csv) %>%
    bind_rows %>% 
    select(sample,cov, cont) %>%
    write_csv(snakemake@output[[1]])
