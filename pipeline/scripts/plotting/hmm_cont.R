library(tidyverse)
library(reshape2)

infile = snakemake@input$cont


a = read_csv(infile) 
ref = read_csv(snakemake@input$ref)
tot_n_snps = nrow(ref) - 261587 - 415546

print("WARNING, nsnp for hcneaden panel")
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

tbl =  a %>% group_by(deam) %>%                              
    summarize(cont=weighted.mean(cont, n_snps),            
          contaminant=sum(n_snps*(cont) / tot_n_snps), 
          endogenous=sum(n_snps*(1-cont)) / tot_n_snps) %>%
    mutate(deam=as.factor(deam))
levels(tbl$deam) = c('deam', 'other')

 
P2 = tbl %>%
    select(deam, contaminant, endogenous) %>%
    gather(k, v, -deam) %>%
    ggplot(aes(x=deam, y=v, fill=k)) +
    geom_col() +
    theme_classic(8) +
    scale_y_continuous(NULL) +
    theme(legend.position='none',
          axis.text.y=element_text(angle=90, hjust=0.5),
          legend.title=element_blank(),
          axis.title.x=element_blank())


ggsave(snakemake@output$png2,P2, width=.9, height=2)


a %>% 
#    group_by(deam) %>%
    summarize(cont=weighted.mean(cont, n_snps), 
                cov=sum(n_snps) / tot_n_snps) %>%
    mutate(sample=snakemake@wildcards$sample) %>%
    write_csv(snakemake@output$table)

