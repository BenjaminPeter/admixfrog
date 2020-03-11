library(Homo.sapiens)
library(GenomicRanges)
#library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg18)
library(tidyverse)


#' annotate bval. Got  this from Martin Petr
annotate_bval <- function(snps){
    bval_folder = '/mnt/expressions/mp/nea-over-time/data/bkgd'
    chain_path = '/mnt/expressions/mp/nea-over-time/data/hg18ToHg19.over.chain'
    bval_files <- list.files(bval_folder, full.names=TRUE, ".*.bkgd")
    
    bval_df_list <- lapply(bval_files, function(filename) {
        read.table(filename, col.names=c("bval", "length")) %>%
            mutate(chr=str_replace(basename(filename), ".bkgd", ""),
                   end=cumsum(length),
                   start=c(1, (end + 1)[-n()])) %>%
            dplyr::select(chr, start, end, bval)
    })

    # lift over to hg 19
    bval_regions_hg18 <- bind_rows(bval_df_list) %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    seqinfo(bval_regions_hg18) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)
    chain <- import.chain(chain_path)
    bval_regions <- liftOver(bval_regions_hg18, chain) %>% unlist


    gsnps = snps %>%
        mutate(chrom=sprintf("chr%s", chrom), end=pos, start=pos) %>% 
        makeGRangesFromDataFrame(T)


    bval_sites <- subsetByOverlaps(gsnps, bval_regions)
    hits <- findOverlaps(gsnps, bval_regions)
    bvals <- IntegerList(split(bval_regions$bval[subjectHits(hits)], queryHits(hits)))
    mcols(bval_sites)[['bval']] = bvals
    mcols(bval_sites[which(elementNROWS(bval_sites$bval) > 1)]) <- NA
    bval_sites$bval <- unlist(bval_sites$bval)

    b_control <- bval_sites %>% as_tibble %>% 
        dplyr::rename(chrom=seqnames) %>% 
        dplyr::select(-width, -strand, -start, -end)
    return(b_control)
    
}

get_bval_regions <- function(){
    bval_folder = '/mnt/expressions/mp/nea-over-time/data/bkgd'
    chain_path = '/mnt/expressions/mp/nea-over-time/data/hg18ToHg19.over.chain'
    bval_files <- list.files(bval_folder, full.names=TRUE, ".*.bkgd")
    
    bval_df_list <- lapply(bval_files, function(filename) {
        read.table(filename, col.names=c("bval", "length")) %>%
            mutate(chr=str_replace(basename(filename), ".bkgd", ""),
                   end=cumsum(length),
                   start=c(1, (end + 1)[-n()])) %>%
            dplyr::select(chr, start, end, bval)
    })

    # lift over to hg 19
    bval_regions_hg18 <- bind_rows(bval_df_list) %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    seqinfo(bval_regions_hg18) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)
    chain <- import.chain(chain_path)
    bval_regions <- liftOver(bval_regions_hg18, chain) %>% unlist
    bval_regions %>% saveRDS("tables/paper/bval_raw.rds")
}

#' annotate bval. Got  this from Martin Petr
annotate_bval_bins <- function(bins){
    bval_regions = readRDS("tables/paper/bval_raw.rds")

    gbins = bins %>%
        group_by(chrom) %>%
        mutate(end = lead(pos) - 1, start=pos) %>%
        ungroup %>%
        filter(!is.na(end)) %>%
        mutate(chrom=sprintf("chr%s", chrom)) %>% 
        makeGRangesFromDataFrame(T)

    hits <- findOverlaps(gbins, bval_regions)

    bval_regions$subjectHits=1:length(bval_regions)
    gbins$queryHits=1:length(gbins)

    h1 = hits %>% as_tibble %>% left_join(bval_regions %>% as_tibble %>% select(start, end, bval, subjectHits))
    h2 = h1 %>% left_join(gbins %>% as_tibble %>% select(queryHits, bin_start=start, bin_end=end, id))
    data = h2 %>% select(-queryHits, -subjectHits) %>% 
        mutate(start=pmax(start, bin_start), end=pmin(end, bin_end)) %>% 
        group_by(id) %>% summarize(bval=weighted.mean(bval, end-start))
    return(bins %>% select(chrom, map, pos, id, n_snps) %>% left_join(data))
}


#' load all SNP positions, and the bin they are in
panel = read_csv(snakemake@input$panel, col_types=cols(chrom='c'))
panel %>% annotate_bval %>% write_csv("tables/paper/bvals.snp.csv.xz")
#snp_file = "admixfrog/error2CAFR/5000/NEA_DEN/altai_hcneaden.snp.xz"
#snps = read_csv(snp_file, col_types=cols(chrom=col_character())) 
#snps2 = annotate_bval(snps)

#'  save output permantently
bin_file = snakemake@input$bin
bins = read_csv(bin_file, col_types=cols(chrom=col_character())) 


deserts <- rbind(
c('1',102200000,114900000),
c('2',201100000,211500000),
c('3',76500000,90500000),
c('7',106300000,124700000),
c('8',53900000,66000000),
c('18',25000000,41800000)
) %>% as.data.frame %>% dplyr::rename(chrom=V1, start=V2, end=V3) %>%
mutate(start=as.integer(as.character(start)), end=as.integer(as.character(end))) 

g_deserts = deserts %>% makeGRangesFromDataFrame

bins %>% annotate_bval_bins %>% 
    mutate(in_desert=bins2 %>% mutate(start=pos, end=pos) %>% 
           makeGRangesFromDataFrame %>% overlapsAny(g_deserts)) %>%
    write_csv(snakemake@output$bvals)
