require(tidyverse)
sgdp= read_csv("config/sgdp.csv.gz") %>% 
    rename(sample=SGDP_ID, pop=Population_ID,
           region=Region,
           country=Country,
           lat=Latitude,
           long=Longitude,
           sex=Gender)
sgdp$age = 0
sgdp$type = 'modern'
sgdp$taxon="MH"
sgdp$source="Mallick2016"
sgdp$culter = 'modern'
sgdp$coverage = 30



ancient = read_csv("config/fu_ext.csv")
ancient$type='ancient'
ancient$region




data =  sgdp %>% bind_rows(ancient) %>%
    mutate(region = ifelse(30 < long & long < 56 & 25 < lat & lat < 45, 
                           "MiddleEastCaucasus", region)) 
data %>% write_csv("config/sgdp_ancient.csv")

