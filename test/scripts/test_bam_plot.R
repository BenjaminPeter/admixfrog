args = commandArgs(T)
require(tidyverse)
a = read_csv(sprintf("%s.bin.xz", args[1]))
a %>% select(-chrom, -pos:-n_snps) %>% 
    gather(pop, p, -map) %>% 
    ggplot(aes(x=map, y=p, fill=pop)) + geom_col(width=as.numeric(args[2]))
ggsave(sprintf("%s.png", args[1]))
