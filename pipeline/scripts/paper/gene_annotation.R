library(Homo.sapiens)
library(GOfuncR)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg18)
library(stringr)
library(tidyverse)

frags = read_csv("tables/paper/all_frags.csv.gz")

all_genes  = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#gfrags = frags %>% filter(sample=='denisova2') %>% 
gfrags = frags %>% #filter(sample=='denisova3') %>% 
    mutate(chrom2=sprintf("chr%s", chrom)) %>% 
    dplyr::rename(start=pos, end=pos_end) %>%
    makeGRangesFromDataFrame(keep.extra.columns=T, seqnames='chrom2')
overlaps = gfrags %>% subsetByOverlaps(all_genes, .) %>%
    as_tibble %>% left_join(as.data.frame(org.Hs.egSYMBOL))

gofunc_df = data.frame(gene_ids=overlaps$symbol, is_candidate=1) %>%
    filter(!is.na(gene_ids))

#gofunc_res = gofunc_df %>% go_enrich(n_randset=1000, gene_len=T)

#gofunc_res %>% arrange(FWER_underrep) %>% head(10) %>% dplyr::select(3, 6)

# this is an unnecessarily huge dependency - but Seqinfo(genome="hg18")
# relies on the internet connection and did not seem to work in
# Jupyter's IRkernel

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

bval_regions_hg18 <- bind_rows(bval_df_list) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
seqinfo(bval_regions_hg18) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)

chain <- import.chain(chain_path)
bval_regions <- liftOver(bval_regions_hg18, chain) %>% unlist

panel = read_csv("ref/ref_hcneaden.csv.xz", col_types=cols(chrom=col_character()))
snps = panel %>% 
    mutate(chrom=sprintf("chr%s", chrom), end=pos) %>% 
    dplyr::select(chrom, start=pos, end) %>% makeGRangesFromDataFrame

snps2 = read_csv("admixfrog/error2CAFR/5000/NEA_DEN/altai_hcneaden.snp.xz", 
                 col_types=cols(chrom=col_character())) %>%
#    group_by(chrom, bin) %>% summarize(pos=dplyr::first(pos)) %>%
#    ungroup %>%
    mutate(chrom=sprintf("chr%s", chrom), end=pos, start=pos) %>% 
    makeGRangesFromDataFrame(T)

bval_sites <- subsetByOverlaps(snps2, bval_regions)
hits <- findOverlaps(snps2, bval_regions)
bvals <- IntegerList(split(bval_regions$bval[subjectHits(hits)], queryHits(hits)))
mcols(bval_sites)[['bval']] = bvals
mcols(bval_sites[which(elementNROWS(bval_sites$bval) > 1)]) <- NA
bval_sites$bval <- unlist(bval_sites$bval)

subsetByOverlaps(bval_sites, gfrags) %>% as_tibble  %>%
    group_by(chrom=seqnames, bin) %>% 
    summarize(bval=mean(bval), pos=dplyr::first(pos))%>% 
    write_csv("b_introgression.txt")

bval_sites %>% as_tibble %>% 
    group_by(chrom=seqnames, bin) %>% 
    summarize(bval=mean(bval), pos=dplyr::first(pos))%>% 
    write_csv("b_control.txt")



deserts <- rbind(
c('chr1',102200000,114900000),
c('chr2',201100000,211500000),
c('chr3',76500000,90500000),
c('chr7',106300000,124700000),
c('chr8',53900000,66000000),
c('chr18',25000000,41800000)
) %>% as.data.frame %>% dplyr::rename(chrom=V1, start=V2, end=V3) %>%
mutate(start=as.integer(as.character(start)), end=as.integer(as.character(end))) %>%
    makeGRangesFromDataFrame

B = b_introgression %>% mutate(is_introgressed=1) %>% right_join(b_control) %>% mutate(is_introgressed=replace_na(is_introgressed, 0)) %>% filter(!is.na(bval)) %>% mutate(bbin=cut(bval, quantile(bval, 0:5/5), include.lowest=T))
B %>% group_by(bbin) %>% summarize(m=mean(is_introgressed)) %>% group_by(bbin) %>% summarize(m=mean(m))

#check for overlap with deserts
