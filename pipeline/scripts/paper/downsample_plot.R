library(glue)
library(scales)
source("scripts/plotting/lib.R")

options(scipen=999)
YSCALE = scale_y_continuous(name="Prob.", breaks=c(.5, 1), expand=expand_scale(0,0))
XSCALE = scale_x_continuous(NULL, breaks=seq(0, 283, 40), expand=expand_scale(0,0))
THEME2 = theme(strip.text.y = element_text(size = 7, angle=180),
               panel.spacing.y = unit(0.05, 'lines'))

fname = snakemake@input$bins
ds_names = snakemake@params$names
CHROM = snakemake@params$chrom
upper = snakemake@params$upper
lower = snakemake@params$lower

plot_fn = bin_colplot_map
if(!is.null(snakemake@params$plottype)){
    if(snakemake@params$plottype == 'pos')
        print("posplot")
        plot_fn = bin_colplot_pos
        XSCALE = list()
}

a = load_bin_data(fname, ds_names) %>% filter(chrom==CHROM) %>%
	bin_to_long %>% 
    mutate(sample=fct_rev(sample),
           variable=str_replace(variable, 'VIN', 'NEA'),
           variable=str_replace(variable, 'CHA', 'NEA'),
           variable=str_replace(variable, 'EUR', 'AFR'),
           variable=str_replace(variable, 'EAS', 'AFR'),
           variable=str_replace(variable, 'AFK', 'AFR'),
           ) %>%
	filter(value>0.001, variable != "AFR")

P = plot_fn(a) + 
    facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + THEME2 + XSCALE
if(!is.null(lower+upper)){
    P = P + coord_cartesian(xlim=c(lower, upper))
} 
ggsave(snakemake@output$png, P, width=3.5, height=1.1, dpi=snakemake@params$dpi)
saveRDS(P, snakemake@output$rds)
