source("scripts/plotting/lib.R")
library(corrplot)
library(viridis)

joint_overview <- function(rlefiles, names,
                           outname
                           state="VIN",
                           type_="state"){

    data = load_rle_data(rlefiles, names) %>%
        filter(target==state, type==type_)
    
    P = data %>% 
        ggplot(aes(x=0, ymin=map, ymax=map_end, color=sample)) +
        geom_linerange(lwd=2, position=position_dodge(.1)) +
        coord_flip() +
        facet_wrap(~chrom, ncol=1, strip='l') + THEME
    ggsave(outname, width=20, height=11)
}

