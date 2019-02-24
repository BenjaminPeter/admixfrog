library(tidyverse)
source("scripts/plotting/lib.R")

debug = F

infile = snakemake@input$rle
pmin = snakemake@params$pmin
pmax = snakemake@params$pmax
lmin = snakemake@params$lmin
lmax = snakemake@params$lmax
bin_size = as.numeric(snakemake@wildcards$bin_size) * 1e-6
outfile = snakemake@output$mapplot
maxgap = snakemake@params$maxgap



if(debug){
	infile = "admixfrog/10000/AFR_ZERO/SS6004480-dedup_twomillion.rle.xz"
	infile2 = "admixfrog/10000/AFR_ZERO/SS6004472-dedup_twomillion.rle.xz"
	infile3 = "admixfrog/10000/AFR_VIN_DEN//SS6004472-dedup_twomillion.rle.xz"
	infile4 = "admixfrog/10000/AFR_VIN_DEN//SS6004472-dedup_archaicadmixture.rle.xz"
	infile5 = "admixfrog/10000/AFR_VIN_DEN/Salkhit_archaicadmixture.rle.xz"
	#infile5 = "admixfrog/10000/AFR_VIN_DEN/Tianyuan_archaicadmixture.rle.xz"

	#R = read_rle(infile)
	R2 = read_rle(infile2)
	R3 = read_rle(infile3)
	#R4 = read_rle(infile4)
	R5 = read_rle(infile5)

	Z= R5 %>% mutate(chrom=sort_chroms(chrom), len=len * 0.01) %>% 
	filter(type!='state', target!="AFR", len>=0*0.01, 
	       gap/len < 11.1) 
} else {
    R = read_rle(infile) %>% filter(type=='state')
    print(levels(R$chrom))
    v = R %>% group_by(target, type) %>% 
    summarize(l=sum(len)) %>% 
    filter(type=='state') %>% 
    ungroup %>% 
    mutate(l=l/sum(l))
    states = v %>% filter(l>pmin, l<pmax) %>% select(target) %>% unlist
    print(states)
    P= R %>% filter(target %in% states, gap/len < maxgap) %>%
	rle_plot_map(minlen=lmin, maxlen=lmax, bin_size=bin_size)
    ggsave(outfile, P, width=20, height=11.5)
}
