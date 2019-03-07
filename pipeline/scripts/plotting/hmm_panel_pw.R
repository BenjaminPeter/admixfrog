source("scripts/plotting/lib.R")
library(corrplot)
library(viridis)
save.image(".rdebug4")

bin_size = as.integer(snakemake@wildcards$bin_size)
panel = snakemake@wildcards$panel
rlefiles = snakemake@input$rle
binfiles = snakemake@input$bins

names = snakemake@config$panels[[panel]]
penalty = as.numeric(snakemake@wildcards$penalty)
state = snakemake@wildcards$target
type_ = snakemake@wildcards$type

l_cutoffs = snakemake@params$lengths


data = load_rle_data(rlefiles, names) %>%
    filter(target==state, type==type_)
bins = load_bin_data(binfiles[1], names[1])
coords =  bins %>% select(chrom, id, pos, map) %>% distinct()

P = data %>% 
    ggplot(aes(x=0, ymin=map, ymax=map_end, color=sample)) +
    geom_linerange(lwd=2, position=position_dodge(.1)) +
    coord_flip() +
    facet_wrap(~chrom, ncol=1, strip='l') + THEME
ggsave(snakemake@output$tracksimple, width=20, height=11)


df =  data %>% 
    rowwise %>% 
    do(id=seq(.$id, .$id_end)) %>% 
    bind_cols(data %>% select(sample, score, map_len), .) %>%
    unnest %>% right_join(coords, .)

MAP_MIN = .2
# second overview plot: genomic positions and lengths, no indiv info
P2= df %>% arrange(-map_len) %>% 
    filter( map_len>MAP_MIN) %>% 
    ggplot(aes(x=pos/1e6, y=1, fill=map_len, color=map_len)) + 
    geom_col() + 
    facet_wrap(~chrom, ncol=1, strip='l') + 
    THEME + 
    xlab("Position (Mb)") +
    ylab("# individuals") + 
    scale_color_viridis_c(trans='log')
ggsave(snakemake@output$trackplot, P2,  width=20, height=11)

P3= df %>% arrange(-map_len) %>% 
    filter( map_len>MAP_MIN) %>% 
    ggplot(aes(x=map, y=1, group=map_len, color=map_len, fill=map_len)) + 
    geom_col() + 
    facet_wrap(~chrom, ncol=1, strip='l') + 
    THEME + 
    xlab("Position (Mb)") +
    ylab("# individuals") 
ggsave(snakemake@output$trackplot, P2,  width=20, height=11)


mycor = function(...)cov(...) %>% cov2cor %>% replace_na(0)
mycov = function(...)cov(...) %>% replace_na(0)
dbin = function(...)( 1 - (dist(t(...), 'binary') %>% as.matrix))
dbin2 = function(x){
    n = ncol(x)
    f = Vectorize(function(i, j) sum(x[,i] * x[,j]))
    m = outer(1:n, 1:n, FUN=f) / 1e6 * bin_size 
    colnames(m) = names(x)
    rownames(m) = names(x)
    m
}

pw_plot <- function(df, l_cutoffs, fun=mycov, ...){
    par(mfrow=c(2,3))
    for(l in l_cutoffs){
        g = df %>% select(id,sample, map_len) %>%
            mutate(run = map_len >l) %>%
            select(-map_len) %>%
            spread(key=sample, value=run, fill=0)  %>%
            do(c=fun(.[,-1], ...))
        cv = g$c[[1]]

        if(l == l_cutoffs[1]) o = hclust(as.dist(1-cv))$order
        corrplot(cv[o,o], diag=F, is.corr=F, 
                 main = sprintf("> %s cM",l),
                 mar=c(0,0,2,0))
    }
}

# the pairwise correlatin plot
png(filename=snakemake@output$pwplot, width=16, height=10, units="in", res=300)
pw_plot(df, l_cutoffs, mycor)
dev.off()
png(filename=snakemake@output$pwplot2, width=16, height=10, units="in", res=300)
pw_plot(df, l_cutoffs, dbin2)
dev.off()

#save.image("pw.rdebug")

