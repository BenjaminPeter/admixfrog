library(tidyverse)
source("scripts/plotting/lib.R")
#save.image("salkhit.rdebug")

#ref = read_csv("ref/ref_archaicadmixture.csv.xz")
base="admixfrog/error/5000/AFR_NEA_DEN/Salkhit_archaicadmixture"



#bin = read_csv(sprintf("%s.bin.xz", base))
#runs = read_csv(sprintf("%s.rle.xz", base))
#snp = read_csv(sprintf("%s.snp.xz", base))

longest_run = runs %>% filter(target==TARGET) %>% arrange(-len) %>% head(1) %>% tail(1)
longest_runs = runs %>% filter(target==TARGET) %>% arrange(-len) %>% head(10)

plot_run <- function(snp, longest_run, TARGET='AFRDEN', ext=c(4e5, 5e4)){
    longest_run_snp = snp %>% 
            filter(chrom==longest_run$chrom, pos>=longest_run$pos - ext[1], 
                   pos < longest_run$pos_end + ext[2]) %>%
            mutate(in_region = pos >= longest_run$pos & pos < longest_run$pos_end) %>%
            left_join(ref %>% select(-map)) 

    lrs0 = longest_run_snp %>% mutate(pafr= AFR_alt / (AFR_alt + AFR_ref),
                                     pnea= VIN_alt / (VIN_alt + VIN_ref),
                                     #peur= EUR_ref / (EUR_alt + EUR_ref),
                                     pden= DEN_alt / (DEN_alt + DEN_ref)) 

    #TRIALS
    x = lrs0 %>%  
        mutate(fixed_den = abs(pden - pafr) == 1,
               fixed_nea = abs(pnea - pafr) == 1)
    x$fixed = 'none'
    x$fixed[x$fixed_nea]  = 'fixed in Neandertal'
    x$fixed[x$fixed_den]  = 'fixed in Denisova'
    x$fixed[x$fixed_den & x$fixed_nea]  = 'both'

    x= x %>% 
        mutate(p_=talt/(tref+talt)) %>% 
        select(snp_id, bin, fixed, pos, p, p_, pnea, pden, pafr) %>% 
        gather('k', 'v', p_:pafr) %>% 
        mutate(target=pos>longest_run$pos & pos < longest_run$pos_end,
               pos=scales::comma(pos)) #%>%
        #filter(!fixed %in% c('both', 'none'))

    x$k = as.factor(x$k)
    levels(x$k) = c('read ratio', 'p[African]', 'p[Denisova]', 'p[Neandertal]')
    x$k = fct_relevel(x$k, levels(x$k[c(2,1,5,3,4)]))


    cols = c("fixed in Neandertal"='#3759a2', "fixed in Denisova"="#ff4000",
             'none'='darkgrey', 'both'='darkgrey')


    P = x %>% ggplot() + 
        #geom_col(aes(x=pos, y=-.05, fill=target)) + 
        geom_col(aes(x=pos, y=v, alpha=target, fill=fixed))  + 
        facet_grid(k~bin, space='free_x', scales='free_x') + 
        theme_classic() + 
        theme(axis.text.x=element_text(angle=90, vjust=.5), strip.text.x=element_text(size=0)) +
        scale_alpha_discrete(range=c(0.3, 1)) +
        geom_area(aes(x=pos, y=p), lwd=1, color='black', data=x %>% filter(k=='read ratio')) +
        scale_x_discrete("SNP position") +
        scale_y_continuous("allele frequency", expand=c(0,0,0,0)) +
        scale_fill_manual(values=cols)
}

