library(glue)
library(scales)
source("scripts/plotting/lib.R")
options(scipen=999)


YSCALE = scale_y_continuous(name="Prob.", breaks=c(.5, 1), expand=expand_scale(0,0))
XSCALE = scale_x_continuous(name="Position (cM)", 
			    breaks=seq(0, 283, 20), expand=expand_scale(0,0),
                labels=comma)
THEME2 = theme(strip.text.y = element_text(size = 8, angle=180))
STUFF = list(XSCALE, YSCALE, THEME2)

ALL_NAMES = c()

# DOWNSAMPLE PLOT  
ds = c(0.0002, 0.01, 0.02, 1)
bs = c(50000, 20000, 20000, 5000)
states=c("VIN_DEN", rep("NEA_DEN", 3))
asc=c("hcneaden")
ds_names =sprintf("%sx", ds * 50)
fname = sprintf("admixfrog/ds%s/%s/%s/altai_%s.bin.xz", ds, bs, states, asc)

cat(fname)
ALL_NAMES = c(ALL_NAMES, fname)

a = load_bin_data(fname, ds_names) %>% filter(chrom==9) %>%
	bin_to_long %>% 
    mutate(sample=fct_rev(sample),
           variable=ifelse(variable=='N', 'NEA', variable),
           variable=str_replace(variable, 'VIN', 'NEA')
           ) %>%
	filter(value>0.01, variable != "NEA") 
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE+ XSCALE + THEME2
ggsave("figures/paper/altai_ds2.png", width=3.5, height=1.5)
ggsave("figures/paper/altai_ds2.pdf", width=3.5, height=1.5)



# REC PLOT
recs = c("AA_Map", "deCODE", "YRI_LD")
fname = sprintf("admixfrog/rec%s/5000/NEA_DEN/altai_hcneaden.bin.xz", recs)
fname = c(fname, "admixfrog/posmode/5000/NEA_DEN/altai_hcneaden.bin.xz")
rec_names = c("AA", "decode", "YRI", "phys")
cat(fname)

ALL_NAMES = c(ALL_NAMES, fname)
a = load_bin_data(fname, rec_names) %>% filter(chrom==9) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=rec_names)) %>%
	filter(value>0.01, variable != "NEA")
P = bin_colplot_pos(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none", 
          strip.text.x = element_text(size = 6)) +
    YSCALE  + THEME2
options(scipen = 999)
ggsave("figures/paper/altai_rec.png", width=3.5, height=1.5)
ggsave("figures/paper/altai_rec.pdf", width=3.5, height=1.5)



#mode
mode = c("error", "posmode", "inbreeding", "errorpos", "nohyper")
mode = c("basic", "nohyper", "inbreeding")
panel_names = c("no err", "fix", "IB", "GTs")
fname = glue("admixfrog/{mode}/5000/NEA_DEN/altai_hcneaden.bin.xz")
fname = c(fname, "admixfrog/gtmode/5000/NEA_DEN/Altai_hcneaden.bin.xz")
cat(fname)

ALL_NAMES = c(ALL_NAMES, fname)

a = load_bin_data(fname, panel_names) %>% filter(chrom==9) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=panel_names)) %>%
	filter(value>0.01, variable != "NEA")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE + THEME2
ggsave("figures/paper/altai_mode.png", P, width=3.5, height=1.5)
ggsave("figures/paper/altai_mode.pdf", P, width=3.5, height=1.5)


# ASCERTAINMENT PLOT 
panel = c("archaicadmixture", "A1240k", "hcneaden", "hcneaden2")#, 'thirdallele')
panel_names = c("AA", "1240", "NEA", "NEA2")# 'thirdallele')
fname = sprintf("admixfrog/basic/5000/NEA_DEN/altai_%s.bin.xz", panel)
ALL_NAMES = c(ALL_NAMES, fname)

a = load_bin_data(fname, panel_names) %>% filter(chrom==9) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=panel_names)) %>%
	filter(value>0.01, variable != "NEA")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/altai_panel.pdf", P, width=3.5, height=1.5)
ggsave("figures/paper/altai_panel.png", P, width=3.5, height=1.5)

# BIN SIZE PLOT
bs = c(2000, 5000, 10000, 20000)
bs_names = bs / 1e6
fname = sprintf("admixfrog/basic/%s/NEA_DEN/altai_hcneaden.bin.xz", bs)
cat(fname)
ALL_NAMES = c(ALL_NAMES, fname)

a = load_bin_data(fname, bs_names) %>% filter(chrom==9) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=bs_names)) %>%
	filter(value>0.01, variable != "NEA")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/altai_bs.pdf", width=3.5, height=1.5)
ggsave("figures/paper/altai_bs.png", width=3.5, height=1.5)

# href PLOT 
refs=c("VIN_DEN", "AFR_NEA_DEN",  "CHA_DEN", "N=ALT+D11_D=DEN+D12")
fname = sprintf("admixfrog/basic/5000/%s/altai_hcneaden.bin.xz", refs)
cat(fname)
ALL_NAMES = c(ALL_NAMES, fname)
href_names = c("VIN", "AFR", 'CHA', "D11")
a = load_bin_data(fname, href_names) %>% filter(chrom==9) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=href_names), 
           variable=ifelse(variable=='N', 'NEA', variable),
           variable=ifelse(variable=='ND', 'NEADEN', variable),
           variable=ifelse(variable=='D', 'DEN', variable),
           variable=str_replace(variable, 'VIN', 'NEA'),
           variable=str_replace(variable, 'CHA', 'NEA'),
           variable=str_replace(variable, 'AFK', 'AFR'),
           variable=str_replace(variable, 'EUR', 'AFR'),
           ) %>%
	filter(value>0.01, !variable %in% c("NEA"))
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/altai_href.png", width=3.5, height=1.5)
ggsave("figures/paper/altai_href.pdf", width=3.5, height=1.5)

ds = 0.02
fake_cont = c(0.05, 0.2, 0.5, 0.8)
fname = sprintf("admixfrog/cont%s_ds%s/20000/NEA_DEN/altai_hcneaden.bin.xz", 
                fake_cont, ds)
cat(fname)
ALL_NAMES = c(ALL_NAMES, fname)

a = load_bin_data(fname, fake_cont) %>% filter(chrom==9) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=fake_cont)) %>%
	filter(value>0.01, variable != "NEA")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/altai_cont.pdf", width=3.5, height=1.5)
ggsave("figures/paper/altai_cont.png", width=3.5, height=1.5)

ds = 0.005
fake_cont = c(0.05, 0.2, 0.5, 0.8)
fname = sprintf("admixfrog/cont%s_ds%s/50000/NEA_DEN/altai_hcneaden.bin.xz", 
                fake_cont, ds)
cat(fname)
ALL_NAMES = c(ALL_NAMES, fname)

a = load_bin_data(fname, fake_cont) %>% filter(chrom==9) %>%
	bin_to_long %>% 
	mutate(sample=factor(sample, levels=fake_cont)) %>%
	filter(value>0.01, variable != "NEA")
P = bin_colplot_map(a) + facet_wrap(~sample, ncol=1, strip='left') +
	theme(legend.position="none") + YSCALE + XSCALE+ THEME2
ggsave("figures/paper/altai_cont2.pdf", width=3.5, height=1.5)
ggsave("figures/paper/altai_cont2.png", width=3.5, height=1.5)
