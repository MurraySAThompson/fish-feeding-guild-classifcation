# Author: Murray SA Thompson
# Contact: murray.thompson@cefas.gov.uk
# Version: 1
# January 2025

# script for: Thompson et al (2024). Modelled and observed fish feeding traits for the North Atlantic and Arctic Oceans (1836-2020) and population estimates of fish with different feeding traits from Northeast Atlantic scientific trawl surveys (1997-2020). Cefas, UK. V1. doi: https://doi.org/10.14466/CefasDataHub.149
# and example data from: https://github.com/MurraySAThompson/fish-feeding-guild-classifcation
# fguilds data from: https://doi.org/10.14466/CefasDataHub.149

# rm(list=ls()); gc()

pkgs = c("tidyverse", "data.table", "vroom", "vegan", "mapplots", "patchwork", "ggpubr") #
for(p in pkgs){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
} 

# set working directory where data are stored
#setwd()

# load example data and taxonomy
example_tx = read_csv('example_taxa.csv')
example_surv = read_csv('example_survey_data.csv')
fguilds = read_csv('feeding guilds.csv') 

####################################################################################################################
### Here is how to load and pre-process all of the fish survey data downloaded from https://data.cefas.co.uk/view/21421
## set directory where files are stored
# setwd()
# 
# # read all survey data txt files:
# list_of_files = list.files(
#   pattern = "\\.txt$", 
#   full.names = TRUE)
# 
# # Read all the files and create a FileName column to store filenames
# DT = rbindlist(sapply(list_of_files, fread, simplify = FALSE),
#                use.names = TRUE, fill=TRUE, idcol = "FileName")
#
# save(DT, file='combined_survey_data.R')
####################################################################################################################

# Gobiidae and Ammodytidae taxa resolved to species across different data:
go_am_tx = c("Pomatoschistus minutus",
             "Hyperoplus lanceolatus", 
             "Gobius niger", 
             "Gobiusculus flavescens", 
             "Gobius paganellus",
             "Ammodytes dubius",
             "Aphia minuta")

## append taxonomy and body mass bins
# replace example_surv with DT if using the full combined_survey_data
alldf = example_surv %>% as_tibble() %>%
  left_join(example_tx, c('SciName'='taxa')) %>%
  rename(taxa=revised_taxa, year=Year, latitude=ShootLat_degdec, longitude=ShootLong_degdec) %>%
  mutate(ices_statistical_rectangle=ices.rect2(longitude, latitude),
         # to iron out missing biomass info
         ind_weight_g = exp(log(LWRa)+LWRb*log(FishLength_cm)),
         DensBiom_kg_Sqkm = case_when(is.na(DensBiom_kg_Sqkm) ~ (ind_weight_g/1000)*DensAbund_N_Sqkm,
                                      TRUE ~ DensBiom_kg_Sqkm),
         # add body mass bins
         bin_num = case_when(ind_weight_g < 6.73e-8 ~ 1,
                             ind_weight_g < 4.38e-7 ~ 2,
                             ind_weight_g < 2.85e-6 ~ 3,
                             ind_weight_g < 1.85e-5 ~ 4,
                             ind_weight_g < 1.20e-4 ~ 5,
                             ind_weight_g < 7.83e-4 ~ 6,
                             ind_weight_g < 5.09e-3 ~ 7,
                             ind_weight_g < 3.31e-2 ~ 8,
                             ind_weight_g < 2.16e-1 ~ 9,
                             ind_weight_g < 1.40e+0 ~ 10,
                             ind_weight_g < 9.12e+0 ~ 11,
                             ind_weight_g < 5.93e+1 ~ 12,
                             ind_weight_g < 3.86e+2 ~ 13,
                             ind_weight_g < 2.51e+3 ~ 14,
                             ind_weight_g < 1.63e+4 ~ 15,
                             ind_weight_g < 1.06e+5 ~ 16,
                             ind_weight_g < 6.90e+5 ~ 17,
                             ind_weight_g < 4.49e+6 ~ 18,
                             ind_weight_g < 2.92e+7 ~ 19,
                             ind_weight_g < 1.90e+8 ~ 20)) %>%
  filter(Gear == 'GOV',
         year %in% c(1997:2020),
         !Survey_Acronym == 'GNSIntOT3',
         DensBiom_kg_Sqkm >0,
         !is.na(bin_num)) %>%
  distinct() %>%
  # resolve taxa, ready for integration with feeding guild information
  mutate(taxa = case_when(species %in% go_am_tx ~ species,
                          family == 'Gobiidae' ~ 'Gobiidae',
                          family == 'Ammodytidae' ~ 'Ammodytes',
                          TRUE ~ species)) %>%
  left_join(dplyr::select(fguilds, bin_number, taxa, f_guild), c('bin_num'='bin_number', 'taxa')) 

# haul-level intra-guild richness, total biomass and abundance
div_gld_haul = alldf %>% # sum at haul-level
  filter(functional_group=='fish') %>%
  group_by(HaulID, latitude, longitude, Survey_Acronym, ices_statistical_rectangle, year, f_guild, taxa) %>%
  summarise(biomass = sum(DensBiom_kg_Sqkm), # b if using 1000 individuals
            abundance = sum(DensAbund_N_Sqkm)) %>% # n if using 1000 individuals
  group_by(HaulID, latitude, longitude, Survey_Acronym, ices_statistical_rectangle, year, f_guild) %>%
  summarise(kgkm = sum(biomass),
            countkm = sum(abundance),
            spp=specnumber(abundance)) %>% 
  ungroup() %>%
  pivot_longer(cols = -(HaulID:f_guild),
               names_to = 'variable',
               values_to = 'value') %>%
  unite(temp, f_guild, variable, sep='#') %>% # introduce 0s
  pivot_wider(names_from = temp, 
              values_from = value,
              values_fill = list(value = 0)) %>%
  pivot_longer(cols = -(HaulID:year),
               names_to = 'variable',
               values_to = 'value') %>%
  separate(variable, c('f_guild', 'variable'), sep='#', remove=T) %>%
  pivot_wider(names_from = variable, 
              values_from = value,
              values_fill = list(value = 0)) %>%
  mutate(f_guild = recode(f_guild, 'NA'='No_guild'))

# classifies X% of biomass
n_clss = alldf %>%
  group_by(f_guild) %>%
  summarise(sum = sum(DensBiom_kg_Sqkm),
            rich = length(unique(taxa))) 
100*(1-(n_clss[nrow(n_clss),'sum']/ sum(n_clss$sum)))

# classifies X% of species
1-(n_clss[nrow(n_clss),'rich']/ sum(n_clss$rich)) 

# plotting
world_shp = sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))

# generate ices_statistical_rectangle-level means for plotting
mydata = div_gld_haul %>% 
  group_by(ices_statistical_rectangle, year, f_guild) %>%
  summarise(kgkm = mean(kgkm),
            countkm = mean(countkm),
            spp= mean(spp)) %>%
  mutate(f_guild = factor(f_guild, levels = c("Planktivore","Benthivore", "Bentho_piscivore", "Piscivore", "No_guild"))) 

# ices statistical rectangle centroids
coords = ices.rect(mydata$ices_statistical_rectangle)
mydata = mydata %>%
  bind_cols(coords)

# define plot area
xlims = range(mydata$lon)
ylims = range(mydata$lat)

# choose which to plot 
toplot = levels(droplevels(unique(mydata$f_guild)[!unique(mydata$f_guild) %in% 'No_guild']))
#toplot = unique(mydata$f_guild)

# generate plots for each feeding guild
plot_list_bmass = list()
for(gld in levels(mydata$f_guild)){ #gld = unique(mydata$f_guild)[1]
  
  pltdat = subset(mydata, f_guild==gld & !is.na(kgkm)) 
  
  Lvals = range(log10(mydata$kgkm+1), na.rm = T)
  breaks = seq(0, round(max(Lvals)))
  labels = as.character(10^breaks)
  lims = range(0, max(Lvals))
  
  p = pltdat %>% ggplot() +
    geom_tile(aes(x=lon, y=lat, fill = log10(kgkm+1)), color=NA) +
    scale_fill_viridis_c(name = 'kg km2',
                         limits = range(log10(pltdat$kgkm+1), na.rm=T),
                         breaks = breaks, 
                         labels=labels) +
    geom_sf(data = world_shp, 
            fill = 'black', 
            color = 'black',
            size = 0.1) +
    coord_sf(xlim = xlims, ylim = ylims)+
    labs(x='', y='', title=gsub('_', ' ', gld)) +
    guides(fill = guide_colourbar(barwidth = 1)) +
    theme(panel.background = element_rect(fill = 'grey80'),
          panel.border = element_rect(colour='black', fill=NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plot_list_bmass[[gld]] = p
}

# plot mean ices-level feeding guild biomass for 2020
do.call(ggarrange, c(plot_list_bmass[toplot], 
                                     ncol=2, nrow=2)) 

# plot biomass of unclassified fish
plot_list_bmass[["No_guild"]]
