source("scripts/plotting/lib.R")

panel = snakemake@wildcards$panel
rlefiles = snakemake@input$rle
binfiles = snakemake@input$bins

names = snakemake@config$panels[[panel]]
penalty = as.numeric(snakemake@wildcards$penalty)
state = snakemake@wildcards$target
type_ = snakemake@wildcards$tt

ages = read_table2("config/ages.yaml", col_names=c("sample", "age")) %>%
    left_join(data.frame(sample=names), .) %>% 
    mutate(age=replace_na(age, 0))

save.image("call_frags.rdebug")


data = load_rle_data(rlefiles, names) %>%
    filter(target==state, type==type_)
bins = load_bin_data(binfiles[1], names[1])
coords =  bins %>% select(chrom, id, pos, map) %>% distinct()

MAP_MIN = snakemake@params$MAP_MIN
MAP_MAX = snakemake@params$MAP_MAX
df =  data %>% 
    rowwise %>% 
    do(id=seq(.$id, .$id_end)) %>% 
    bind_cols(data %>% select(sample, score, map_len), .) %>%
    unnest %>% right_join(coords, .) %>%
    arrange(-map_len) %>% 
    filter(map_len >= MAP_MIN) %>%
    mutate(length = pmin(map_len, MAP_MAX)) 

frags = df %>% filter(map_len >= MAP_MIN) %>% 
    mutate(too_big = map_len > MAP_MAX) %>% 
    mutate(too_big = ifelse(too_big, sample, "NO")) %>%
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

frag_df %>% write_csv(snakemake@output$frags)
