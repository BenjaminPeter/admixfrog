source("scripts/comparison_plot.R")
library(corrplot)
library(viridis)

bin_size = as.integer(snakemake@wildcards$bin_size)
panel = snakemake@wildcards$panel
infile = snakemake@input$bin
snpfile = snakemake@input$snp
names = snakemake@config$panels[[panel]]
cutoff = as.numeric(snakemake@wildcards$cutoff)
l_cutoffs = snakemake@params$lengths / bin_size * 1000

save.image(".rdebug")


TRACK = strsplit(snakemake@wildcards$TRACK, "_")[[1]]
data = load_data(infile, names) 
if(cutoff > 0){
    data$TRACK = rowSums(data[,TRACK]) > cutoff
}else{
    data$TRACK = rowSums(data[,TRACK]) <  (-cutoff)
}

coords <- data %>% select(chrom, bin_pos, bin_id)
mycov = function(...)cov(...) %>% cov2cor %>% replace_na(0)

df = lapply(l_cutoffs, get_rundf, data=data)
names(df) = l_cutoffs
df = df %>% bind_rows(.id="Length") %>% mutate(Length=as.integer(Length) * bin_size / 1000)
x = df %>% select(-bin_id) %>% group_by(Length) %>% do(c=mycov(.[,-1])) 
o = hclust(as.dist(1-x$c[[1]]))$order

png(filename=snakemake@output$pwplot, width=16, height=10, units="in", res=300)
par(mfrow=c(2,3))
for(i in 1:6)
    corrplot(x$c[[i]][o,o], diag=F, is.corr=F, main = sprintf("> %s kb",x$Length[i]), mar=c(0,0,2,0))

dev.off()


X = df %>% gather(sample, run, -1:-2)
Y = X %>% filter(run) %>% group_by(Length, bin_id)  %>% summarize(n=n()) %>% arrange(-n) 
Z = Y %>% left_join(coords)                                                              
Z %>% filter( n>=1, Length > 0)  %>% 
    ungroup %>%
    arrange(-n) %>% 
    ggplot(aes(x=bin_pos, y=n, color=Length)) + 
    geom_col(position="identity") + 
    facet_wrap(~chrom, ncol=2, strip.position="left") +
    xlab("Position") +
    ylab("# individuals") + 
    scale_color_viridis_c(trans="log")

ggsave(snakemake@output$trackplot, width=20, height=11)

#save.image("pw.rdebug")

