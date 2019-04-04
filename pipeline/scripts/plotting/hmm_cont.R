library(tidyverse)
library(reshape2)

infile = snakemake@input$cont


a = read_csv(infile) 
ref = read_csv(snakemake@input$ref)
tot_n_snps = nrow(ref)
print(tot_n_snps)

MIN = 100
a = a %>% filter(n_snps > MIN) %>% select(-lib)
a = a %>% arrange(-n_snps) %>% mutate(rg=factor(rg, unique(rg)))
b = a %>% mutate(contam.=n_snps*cont, endogenous =n_snps-contam.)  %>% 
    select(-n_snps, -cont)                      
b %>% melt(id.vars=c("rg", "len_bin", "deam")) %>% 
    ggplot(aes(x=len_bin, y=value, fill=variable)) + 
    geom_col() + 
    #facet_wrap(~rg * len_bin, ncol=1, strip.position="left") +
    facet_grid(rg ~deam) +
    coord_flip() + 
    theme(panel.spacing = unit(0, "lines")) +
    theme_classic(7) + 
    xlab("length bin") + ylab("# reads") +
    ggtitle(snakemake@wildcards$sample)
ggsave(snakemake@output$png, width=6, height=nrow(a) * .4 +1, limitsize=F)

P2 = a %>% group_by(deam) %>% 
    summarize(cont=weighted.mean(cont, n_snps), 
              contaminant=sum(n_snps*(cont) / tot_n_snps), 
              endogenous=sum(n_snps*(1-cont)) / tot_n_snps)%>%
    select(deam, contaminant, endogenous) %>%
    gather(k, v, -deam) %>%
    ggplot(aes(x=deam, y=v, fill=k)) +
    geom_col() +
    theme_classic(14) +
    ylab("coverage") +
    theme(legend.position='bottom',
          legend.title=element_blank(),
          axis.title.x=element_blank())


ggsave(snakemake@output$png2,P2, width=4, height=4, limitsize=F)

