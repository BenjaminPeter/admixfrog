require(tidyverse)
bed = read_tsv(snakemake@input$bed,
		 col_names=c("chrom", "pos0", "pos", "ref", "alt"),
		 col_types="ciicc")

chroms = c(1:22, "X")
read_map <- function(chrom){
	map_file = sprintf("recs/maps_b37/maps_chr.%s", chrom)
	v <- read_delim(map_file, delim=" ")
	v$chrom <- chrom
	v
}


rec = lapply(chroms, read_map) %>% bind_rows %>%
	rename(pos="Physical_Pos")

maps <- names(rec)[2:8]

approx_chrom <- function(bed, rec, map="deCODE"){
	s <- c()
	for(c in chroms){
		s <- c(s, approx(x=rec$pos[rec$chrom==c], 
			         y=rec[[map]][rec$chrom==c], 
			         xout=bed$pos[bed$chrom==c])$y)
	}
	bed[[map]] <- s
	bed
}
print(snakemake@output$bed)

for(m in maps)
	bed <- approx_chrom(bed, rec, m)


write_csv(bed, snakemake@output$bed)
