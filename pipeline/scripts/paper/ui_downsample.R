source("scripts/plotting/lib.R")

YSCALE = scale_y_continuous(name="Prob.", breaks=c(.5, 1), expand=expand_scale(0,0))
XSCALE = scale_x_continuous(name="Position (cM)", 
			    breaks=seq(0, 283, 20), expand=expand_scale(0,0))
THEME2 = theme(strip.text.y = element_text(size = 6))
STUFF = list(XSCALE, YSCALE, THEME2)


# DOWNSAMPLE PLOT
ds = c(0.002, 0.005, 0.01, 1)
ds_names =sprintf("%sx", ds * 40)
fname = sprintf("admixfrog/ds%s/20000/NEA_DEN/UstIshim_archaicadmixture.bin.xz", ds)
print(fname)

a = load_bin_data(fname, ds_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	filter(value>0.01, variable != "AFR") %>%
    mutate(sample=fct_rev(sample))
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE+ XSCALE 
ggsave("figures/paper/ui_ds.pdf", width=8, height=1.8)
ggsave("figures/paper/ui_ds.png", width=8, height=1.8)

# DOWNSAMPLE PLOT 2
ds = c(0.00025, 0.0025, 0.01, 1)
bs = c(20000, 50000, 5000, 5000)
states=c("AFR_NEA_DEN", rep("AFR_NEA_DEN", 3))
asc=c("archaicadmixture", rep("archaicadmixture", 3))
ds_names =sprintf("%sx", ds * 40)
fname = sprintf("admixfrog/ds%s/%s/%s/UstIshim_%s.bin.xz", ds, bs, states, asc)

a = load_bin_data(fname, ds_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	filter(value>0.01, variable != "AFR") %>%
    mutate(sample=fct_rev(sample))
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE+ XSCALE
ggsave("figures/paper/ui_ds2.png", width=8, height=1.8)
ggsave("figures/paper/ui_ds2.pdf", width=8, height=1.8)


# BIN SIZE PLOT
bs = c(2000, 5000, 10000, 20000)
bs_names = bs / 1e6
fname = sprintf("admixfrog/basic/%s/AFR_NEA_DEN/UstIshim_archaicadmixture.bin.xz", bs)
cat(fname)

a = load_bin_data(fname, bs_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=bs_names)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE
ggsave("figures/paper/ui_bs.pdf", width=8, height=1.8)
ggsave("figures/paper/ui_bs.png", width=8, height=1.8)

# REC PLOT
recs = c("AA_Map", "deCODE", "YRI_LD", "CEU_LD")
fname = sprintf("admixfrog/rec%s/5000/AFR_NEA_DEN/UstIshim_archaicadmixture.bin.xz", recs)
fname = c(fname, "admixfrog/posmode/5000/AFR_NEA_DEN/UstIshim_archaicadmixture.bin.xz")
rec_names = c(recs, "physical")
cat(fname)
a = load_bin_data(fname, rec_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=rec_names)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_pos(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none", 
          strip.text.x = element_text(size = 6)) +
    YSCALE + XSCALE
ggsave("figures/paper/ui_rec.png", width=8, height=1.8)
ggsave("figures/paper/ui_rec.pdf", width=8, height=1.8)

# GT PLOT
bs = c(2000, 5000, 10000, 20000)
bs_names = bs / 1e6
fname = sprintf("admixfrog/gtmode/%s/AFR_NEA_DEN/Ust_Ishim_archaicadmixture.bin.xz", bs)
cat(fname)
a = load_bin_data(fname, bs_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=bs_names)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE
ggsave("figures/paper/ui_gt.pdf", width=8, height=1.8)
ggsave("figures/paper/ui_gt.png", width=8, height=1.8)

# ASCERTAINMENT PLOT
panel = c("archaicadmixture", "A3700k", "A1240k", "thirdallele")
panel_names = c("AA", "3.7M", "1240", "NEA")
fname = sprintf("admixfrog/basic/5000/AFR_NEA_DEN/UstIshim_%s.bin.xz", panel)
a = load_bin_data(fname, panel_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=panel_names)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE
ggsave("figures/paper/ui_panel.pdf", P, width=8, height=1.8)
ggsave("figures/paper/ui_panel.png", P, width=8, height=1.8)

L = read_tsv("laurits/Ancient_Archaicsegments.txt") %>% 
	filter(name=="Ust_ishim", chrom==1) %>% select(start, end)
rec = read_csv("ref/ref_archaicadmixture.csv.xz") %>% filter(chrom==1)
L$map = approx(rec$pos, rec$map, L$start)$y
L$map_end = approx(rec$pos, rec$map, L$end)$y 
L$src = "Skov"
P= L %>% ggplot(aes(xmin=map, xmax=map_end, ymin=0, ymax=1)) + 
	geom_rect() +
	facet_wrap(~src, strip='left') + 
	coord_cartesian(xlim=range(0,285.8846)) + THEME +
	YSCALE  + XSCALE
ggsave("figures/paper/ui_skov.pdf", P, width=8, height=.8)
ggsave("figures/paper/ui_skov.png", P, width=8, height=.8)
