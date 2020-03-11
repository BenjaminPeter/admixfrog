source("scripts/plotting/lib.R")
load("fig1.rdebug")

ref = read_csv("ref/ref_archaicadmixture.csv.xz", 
               col_type=cols(chrom=col_factor(), pos=col_integer()))
names = str_split(snakemake@input$snps, "[/_]", simplify=T)[,7]
rle = load_rle_data(snakemake@input$rle, names) %>% 
    filter(target=="AFRNEA", type=='het') 
S = load_snp_data(snakemake@input$snps, names)


rle = read_csv("stats/frags/mh/state_NEA_0.1/AFR_NEA_DEN/2000/paper_weu2c_archaicadmixture.mcmc.gz")

X = rle %>% 
    rowwise %>% 
    do(bin=seq(.$start, .$end)) %>% 
    bind_cols(rle, .) %>% 
    select(sample, bin, len, fid, M, age) %>% unnest
Z = X %>% inner_join(S) %>% inner_join(ref, by=c('chrom', 'pos'))

abba = function(VIN, CHA, p, PAN)(1-VIN) * CHA * p * (1-PAN) + VIN * (1-CHA) * (1-p) * PAN
baba = function(a,b,c,d)abba(b,a,c,d)

D = Z %>% mutate(
                 VIN = VIN_alt / (VIN_ref + VIN_alt),
                 CHA = CHA_alt / (CHA_ref + CHA_alt),
                 AFR = AFR_alt / (AFR_ref + AFR_alt),
                 ALT = ALT_alt / (ALT_ref + ALT_alt),
                 UST = UST_alt / (UST_ref + UST_alt),
                 PAN = PAN_alt / (PAN_ref + PAN_alt)) %>%
    mutate(ABBA = abba(VIN, CHA, p, PAN),
           BABA = baba(VIN, CHA, p, PAN)) %>%
    mutate(ABBA2 = abba(PAN, AFR, p, ALT),
           BABA2 = baba(PAN, AFR, p ,ALT)) %>%
    mutate(ABBA3 = abba(PAN, AFR, VIN, ALT),
           BABA3 = baba(PAN, AFR, VIN, ALT) 
           ) %>%
    mutate(F41 = ABBA3- BABA3,
           F42 = ABBA2 - BABA2) %>%
    mutate(l2 = floor(len/50))
    
V = D %>% 
    filter(len >= 00) %>%
    group_by(sample,) %>%
    summarize(
          ABBA = sum(ABBA, na.rm=T),
          BABA = sum(BABA, na.rm=T),
          ABBA2 = sum(ABBA2, na.rm=T),
          BABA2 = sum(BABA2, na.rm=T),
          F41 = mean(F41, na.rm=T),
          F42 = mean(F42, na.rm=T),
          D = (ABBA - BABA) / (ABBA + BABA),
          D2 = (ABBA2 - BABA2) / (ABBA2 + BABA2)
           ) %>% left_join(ages) %>% mutate(age=replace_na(age,0))


data = S %>% select(sample, chrom, pos, p, tref, talt) %>% left_join(ref)
    

D2 = data %>% mutate(
                    p = talt / (talt+tref),
                 VIN = VIN_alt / (VIN_ref + VIN_alt),
                 CHA = CHA_alt / (CHA_ref + CHA_alt),
                 AFR = AFR_alt / (AFR_ref + AFR_alt),
                 ALT = ALT_alt / (ALT_ref + ALT_alt),
                 UST = UST_alt / (UST_ref + UST_alt),
                 PAN = PAN_alt / (PAN_ref + PAN_alt)) %>%
    mutate(ABBA = abba(VIN, CHA, p, PAN),
           BABA = baba(VIN, CHA, p, PAN),
           ABBA4 = abba(VIN, CHA, p, PAN),
           BABA4 = baba(VIN, CHA, p, PAN),
           ABBA5 = abba(VIN, ALT, p, PAN),
           BABA5 = baba(VIN, ALT, p, PAN),
           ABBA2 = abba(VIN, p, CHA, PAN),
           BABA2 = baba(VIN, p, CHA, PAN),
           ABBA3 = abba(VIN, ALT, CHA, PAN),
           BABA3 = baba(VIN, ALT, CHA, PAN),
           ) %>%
    mutate(F41 = ABBA3- BABA3,
           F42 = ABBA2 - BABA2,
           F43 = ABBA4 - BABA4,
           F44 = ABBA5 - BABA5)

V2 = D2 %>%
    group_by(sample,) %>%
    summarize(
          ABBA = sum(ABBA, na.rm=T),
          BABA = sum(BABA, na.rm=T),
          ABBA2 = sum(ABBA2, na.rm=T),
          BABA2 = sum(BABA2, na.rm=T),
          ABBA3 = sum(ABBA3, na.rm=T),
          BABA3 = sum(BABA3, na.rm=T),
          sdF4r = sd(F42/F41, na.rm=T),
          F41 = mean(F41, na.rm=T),
          F42 = mean(F42, na.rm=T),
          F43 = mean(F43, na.rm=T),
          F44 = mean(F44, na.rm=T)
           )

V3 = V2 %>% left_join(ages) %>% mutate(age=replace_na(age,0))

           


