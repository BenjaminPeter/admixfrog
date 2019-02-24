source("scripts/plotting/lib.R")
library(corrplot)
library(viridis)
save.image(".rdebug4")

bin_size = as.integer(snakemake@wildcards$bin_size)
panel = snakemake@wildcards$panel
infiles = snakemake@input$rle
snpfile = snakemake@input$snp
names = snakemake@config$panels[[panel]]
cutoff = as.numeric(snakemake@wildcards$cutoff)
state = snakemake@wildcards$target
type_ = snakemake@wildcards$type


data = load_rle_data(infile, names) %>%
    filter(target==state, type==type_)

coords <- data %>% select(chrom, pos, map)
mycov = function(...)cov(...) %>% cov2cor %>% replace_na(0)

#save.image(".rdebug")
#stop()
df = lapply(l_cutoffs, get_rundf, data=data)
names(df) = l_cutoffs
df = df %>% bind_rows(.id="Length") %>% mutate(Length=as.integer(Length) * bin_size / 1e6)
x = df %>% select(-id) %>% group_by(Length) %>% do(c=mycov(.[,-1])) 
o = hclust(as.dist(1-x$c[[1]]))$order




# the pairwise correlatin plot
png(filename=snakemake@output$pwplot, width=16, height=10, units="in", res=300)
par(mfrow=c(2,3))
for(i in 1:length(snakemake@params$lengths))
    corrplot(x$c[[i]][o,o], diag=F, is.corr=F, main = sprintf("> %s cM",x$Length[i]), mar=c(0,0,2,0))

dev.off()

#save.image("pw.rdebug")

# Run lengths
X = df %>% gather(sample, run, -1:-2)
Y = X %>% filter(run) %>% group_by(Length, id)  %>% summarize(n=n()) %>% arrange(-n) 
Z = Y %>% left_join(coords)                                                              
P = Z %>% filter( n>=1, Length > 0)  %>% 
    ungroup %>%
    arrange(-n) %>% 
    ggplot(aes(x=map, y=n, color=Length)) + 
    geom_col(position="identity") + 
    facet_wrap(~chrom, ncol=2, strip.position="left") +
    xlab("Position") +
    ylab("# individuals") + 
    scale_color_viridis_c(trans="log")

ggsave(snakemake@output$trackplot, width=20, height=11)


