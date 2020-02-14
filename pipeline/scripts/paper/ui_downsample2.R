library(glue)
library(scales)
source("scripts/plotting/lib.R")

options(scipen=999)
YSCALE = scale_y_continuous(name="Prob.", breaks=c(.5, 1), expand=expand_scale(0,0))
XSCALE = scale_x_continuous(name="Position (cM)", 
			    breaks=seq(0, 283, 40), expand=expand_scale(0,0),
                labels=comma)
THEME2 = theme(strip.text.y = element_text(size = 7, angle=180))
STUFF = list(XSCALE, YSCALE, THEME2)


# DOWNSAMPLE PLOT  // DONE
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
	theme(legend.position="none") + YSCALE+ XSCALE + THEME2
ggsave("figures/paper/ui_ds2.png", width=3.5, height=1.5)
ggsave("figures/paper/ui_ds2.pdf", width=3.5, height=1.5)



# REC PLOT
recs = c("AA_Map", "deCODE", "YRI_LD")
fname = sprintf("admixfrog/rec%s/5000/AFR_NEA_DEN/UstIshim_archaicadmixture.bin.xz", recs)
fname = c(fname, "admixfrog/posmode/5000/AFR_NEA_DEN/UstIshim_archaicadmixture.bin.xz")
rec_names = c("AA", "dCode", "YRI", "phys")
cat(fname)
a = load_bin_data(fname, rec_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=rec_names)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_pos(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none", 
          strip.text.x = element_text(size = 6)) +
    YSCALE  + THEME2
options(scipen = 999)
ggsave("figures/paper/ui_rec.png", width=3.5, height=1.5)
ggsave("figures/paper/ui_rec.pdf", width=3.5, height=1.5)


#ala fu2014

L = read_tsv("laurits/Ancient_Archaicsegments.txt") %>% 
	filter(name=="Ust_ishim", chrom==1) %>% select(start, end)
rec = read_csv("ref/ref_archaicadmixture.csv.xz") %>% filter(chrom==1)
s = read_csv("samples/Ust_Ishim_archaicadmixture.in.xz")
fu = rec %>% filter(ALT_alt>0, AFR_alt==0) %>% left_join(s) %>% filter(!is.na(talt), talt>0)
L$map = approx(rec$pos, rec$map, L$start)$y
L$map_end = approx(rec$pos, rec$map, L$end)$y 
L$src = "Skov"
L2 = fu %>% mutate(map_end = map+0.05, src='Fu') %>% select(map, map_end, src) %>% bind_rows(L, .)
P= L2 %>% ggplot(aes(xmin=map, xmax=map_end, ymin=0, ymax=1)) + 
	geom_rect() +
	facet_wrap(~src, strip='left', ncol=1) + 
	coord_cartesian(xlim=range(0,285.8846))                         + 
    THEME +
	YSCALE  + XSCALE+ THEME2
ggsave("figures/paper/ui_skov.pdf", P, width=3.5, height=1.0)
ggsave("figures/paper/ui_skov.png", P, width=3.5, height=1.0)

#mode
mode = c("basic", "error", "posmode", "inbreeding", "errorpos", "nohyper")
mode = c("basic", "error", "nohyper")
panel_names = c("basic", "error", "fix", "GTs")
fname = glue("admixfrog/{mode}/5000/AFR_NEA_DEN/UstIshim_archaicadmixture.bin.xz")
fname = c(fname, "admixfrog/gtmode/5000/AFR_NEA_DEN/Ust_Ishim_archaicadmixture.bin.xz")
a = load_bin_data(fname, panel_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=panel_names)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE + THEME2
ggsave("figures/paper/ui_mode.pdf", P, width=3.5, height=1.5)
ggsave("figures/paper/ui_mode.png", P, width=3.5, height=1.5)


# ASCERTAINMENT PLOT 
panel = c("archaicadmixture", "A3700k", "A1240k", "hcneaden")#, 'thirdallele')
panel_names = c("AA", "3.7M", "1240", "NEA")# 'thirdallele')
fname = sprintf("admixfrog/basic/5000/AFR_NEA_DEN/UstIshim_%s.bin.xz", panel)
a = load_bin_data(fname, panel_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=panel_names)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/ui_panel.pdf", P, width=3.5, height=1.5)
ggsave("figures/paper/ui_panel.png", P, width=3.5, height=1.5)

# GT PLOT 
bs = c(2000, 5000, 10000, 20000)
bs_names = bs / 1e6
fname = sprintf("admixfrog/gtmode/%s/AFR_NEA_DEN/Ust_Ishim_archaicadmixture.bin.xz", bs)
a = load_bin_data(fname, bs_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=bs_names)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/ui_gt.pdf", width=3.5, height=1.5)
ggsave("figures/paper/ui_gt.png", width=3.5, height=1.5)
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
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/ui_bs.pdf", width=3.5, height=1.5)
ggsave("figures/paper/ui_bs.png", width=3.5, height=1.5)

# href PLOT 
refs=c("AFR_NEA", "AFK_NEA", "EUR_NEA", "EUR_VIN")
fname = sprintf("admixfrog/nohyper/5000/%s_DEN/Ust_Ishim_archaicadmixture.bin.xz", refs)
cat(fname)
href_names = c("AFR", "AFK", "EUR", "VIN")
a = load_bin_data(fname, href_names) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=href_names), 
           variable=str_replace(variable, 'VIN', 'NEA')) %>%
	filter(value>0.01, !variable %in% c("AFR", "AFK", "EUR"))
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/ui_href.pdf", width=3.5, height=1.5)
ggsave("figures/paper/ui_href.png", width=3.5, height=1.5)

ds = .1
fake_cont = c(0.05, 0.2, 0.5, 0.8)
fname = sprintf("admixfrog/cont%s_ds%s/5000/AFR_NEA_DEN/UstIshim_hcneaden.bin.xz", 
                fake_cont, ds)
cat(fname)

a = load_bin_data(fname, fake_cont) %>% filter(chrom==1) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=fake_cont)) %>%
	filter(value>0.01, variable != "AFR")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/ui_cont.pdf", width=3.5, height=1.5)
ggsave("figures/paper/ui_cont.png", width=3.5, height=1.5)
