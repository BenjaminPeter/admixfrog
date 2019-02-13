suppressPackageStartupMessages({
library(admixr)
library(glue)
library(reshape2)
library(tidyverse)
source("scripts/load.R")
})

map_ = snakemake@params$recmap

print(snakemake@output$csv)



bed <- load_bed(snakemake@input$bed) %>%
    select(chrom, pos, ref, alt, map=snakemake@params$recmap)

target <- read_csv(snakemake@input$sample_stats,
		   col_names=c("chrom", "pos", "rg", "deam", "tref", "talt", "tdeam", "tother"),
		   col_types="cicciiiii", skip=1) %>%
		   mutate(pos = pos + 1) %>%
	    mutate(deam = as.factor(deam))
levels(target$deam) <- c("nodeam", "deam")
target <- target %>%  mutate(lib = paste(rg, deam, sep="_")) 


bed %>% 
    inner_join(target) %>%
    select(chrom, pos, map, lib, tref, talt) %>%
    write_csv(snakemake@output$csv)
