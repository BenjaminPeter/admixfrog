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
type_ = snakemake@wildcards$tt

l_cutoffs = snakemake@params$lengths

ages = read_table2("config/ages.yaml", col_names=c("sample", "age")) %>%
    left_join(data.frame(sample=names), .) %>% 
    mutate(age=replace_na(age, 0))


data = load_rle_data(rlefiles, names) %>%
    filter(target==state, type==type_)
bins = load_bin_data(binfiles[1], names[1])
coords =  bins %>% select(chrom, id, pos, map) %>% distinct()

P = ggplot(mapping=aes(x=0, ymin=map, ymax=map_end, color=sample)) +
    bg_chrom(flip=T) + 
    geom_linerange(data=data, lwd=2, position=position_dodge(.1)) +
    coord_flip() +
    facet_wrap(~chrom, ncol=1, strip='l') + THEME
ggsave(snakemake@output$tracksimple, width=20, height=11)

MAP_MIN = 0
MAP_MAX = 4
df =  data %>% 
    rowwise %>% 
    do(id=seq(.$id, .$id_end)) %>% 
    bind_cols(data %>% select(sample, score, map_len), .) %>%
    unnest %>% right_join(coords, .) %>%
    arrange(-map_len) %>% 
    filter(map_len >= MAP_MIN) %>%
    mutate(length = pmin(map_len, 2)) 

frags = df %>% filter(map_len >= MAP_MIN) %>% 
    mutate(too_big = map_len > MAP_MAX) %>% 
    mutate(too_big = ifelse(too_big, sample, too_big)) %>%
    arrange(too_big, chrom, pos) %>% 
    group_by(too_big, chrom) %>% 
    mutate(gap=id-lag(id) > 1, gap=replace_na(gap, 1)) %>%
    ungroup %>%
    mutate(
           frag=cumsum(gap), 
           frag=interaction(chrom, frag, too_big, drop=T),
           frag=as.integer(frag))

frag_df =  data %>% left_join(frags %>% select(chrom, id, sample, frag, too_big)) %>%
    left_join(ages)

frag_df %>% saveRDS("frags.rds")

#est age of frag
ll <- function(t, k, frag1, trunc=0.2){
    r=sum(dexp(frag1$map_len, pmax(0, (t-frag1$age) / k), log=T) - 
        pexp(trunc, pmax(0, (t-frag1$age) / k), log=T, lower = F)
    )
    return(pmax(r, -1e6))
}
opt_frag <- function(frag, trunc=0.2){
    o <- optim(fn=function(x)ll(x[1], x[2], frag1=frag, trunc=trunc),
               par=c(100000, 5000),
               control=list(fnscale=-1),
               method="L-BFGS-B",
               lower=c(00000, 2000), upper=c(200000, 8000))
    return(tibble(t=o$par[1], k=o$par[2], ll=o$value))
}
opt_frag2 <- function(frag, trunc=0.2, k=5000){
    o <- optim(fn=function(x)ll(x[1], k, frag1=frag, trunc=trunc),
               par=c(100000),
               control=list(fnscale=-1),
               method="L-BFGS-B",
               lower=c(30000), upper=c(200000))
    return(tibble(t=o$par[1], k=o$par[2], ll=o$value))
}





# second overview plot: genomic positions and lengths, no indiv info
P2= ggplot() + 
    bg_chrom() +
    geom_col(data=df, mapping=aes(x=map, y=1, fill=length, color=length)) + 
    facet_wrap(~chrom, ncol=1, strip='l') + 
    THEME + 
    xlab("Position (Mb)") +
    ylab("# individuals") + 
    scale_color_viridis_c(trans='log', name="length (cM)", aesthetics=c("color", "fill"))
ggsave(snakemake@output$trackplot, P2,  width=20, height=11)

P3= df %>% arrange(-map_len) %>% 
    filter( map_len>MAP_MIN) %>% 
    ggplot(aes(x=map, y=1, group=map_len, color=sample, fill=sample)) + 
    geom_col() + 
    facet_wrap(~chrom, ncol=1, strip='l') + 
    THEME + 
    xlab("Position (Mb)") +
    ylab("# individuals") 
ggsave(snakemake@output$trackfull, P3,  width=20, height=11)


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

get_cmat <- function(df, l, fun=mycov){
    g = df %>% select(id,sample, map_len) %>%
        mutate(run = map_len >l) %>%
        select(-map_len) %>%
        spread(key=sample, value=run, fill=0)  %>%
        do(c=fun(.[,-1]))
    cv = g$c[[1]]
}


pw_plot <- function(df, l_cutoffs, fun=mycov, ...){
    par(mfrow=c(2,3))
    for(l in l_cutoffs){
        cv = get_cmat(df, l, fun)

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

