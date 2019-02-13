source("~/programs/admixfrog/plotting/comparison_plot.R")

bin_size = as.integer(snakemake@wildcards$bin_size)
infile = snakemake@input$bin
panel = snakemake@wildcards$panel
names = snakemake@config$panels[[panel]]
cutoff = as.numeric(snakemake@wildcards$cutoff)
n_samples = length(names)

TRACK = strsplit(snakemake@wildcards$TRACK, "_")[[1]]

data = load_data(infile, names)
data$TRACK = rowSums(data[,TRACK]) > cutoff

R = rle_fit_pars(data)
P = rle_fit_plot(data, R) + ggtitle(snakemake@output$rleplot)
P + facet_wrap(~sample, scale="free", strip.position="left")
ggsave(snakemake@output$rleplot, width=20, height=10)

P3= P + scale_y_log10()
ggsave(snakemake@output$rlelogplot, width=20, height=10)

P2 =  plot_m_gamma(R, bin_size, 29) + ggtitle(snakemake@output$rleplot)
ggsave(snakemake@output$gammaplot, width=10, height=1 + (n_samples)  * .8, limitsize=F)
