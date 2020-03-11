library(GOfuncR)
library(tidyverse)
do_gofuncr <- function(frags){
    all_genes  = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

    gfrags = frags %>% 
        dplyr::select(-genes) %>%
        mutate(chrom2=sprintf("chr%s", chrom)) %>% 
        dplyr::rename(start=pos, end=pos_end) %>%
        makeGRangesFromDataFrame(keep.extra.columns=T, seqnames='chrom2')
    overlaps = gfrags %>% subsetByOverlaps(all_genes, .) %>%
        as_tibble %>% left_join(as.data.frame(org.Hs.egSYMBOL))

    gofunc_df = data.frame(gene_ids=overlaps$symbol, is_candidate=1) %>%
        filter(!is.na(gene_ids))
    
    gofunc_res = gofunc_df %>% go_enrich(n_randset=1000, gene_len=T)
    return(gofunc_res)
}

#' load all identified introgression tracts
frags = read_csv("tables/paper/clean/tableS2_frags.csv.gz") %>% 
    rename(start=id, end=id_end) %>% tbl_interval

gofunc_res = do_gofuncr(frags)

#' TODO save output permantently
