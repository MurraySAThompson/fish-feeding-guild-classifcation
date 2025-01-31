# Author: Murray SA Thompson
# Contact: murray.thompson@cefas.gov.uk
# Version: 1
# January 2025

# data from: Thompson et al (2024). Modelled and observed fish feeding traits for the North Atlantic and Arctic Oceans (1836-2020) and population estimates of fish with different feeding traits from Northeast Atlantic scientific trawl surveys (1997-2020). Cefas, UK. V1. doi: https://doi.org/10.14466/CefasDataHub.149
# and example data from: https://github.com/MurraySAThompson/fish-feeding-guild-classifcation

# rm(list=ls()); gc()

pkgs = c("tidyverse", "vegan", "patchwork", "ggpubr", "betapart",
         "scales", "geosphere", "effects", "multcomp") 

for(p in pkgs){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# set working directory where data are stored
#setwd()

# load example data
example_df = read_csv('example data for cluster analyses.csv')
example_samps = read_csv('example samples for cluster analyses.csv')


# use only directly observed prey mass and count data ('Yes') or those including predictions from models ('No')
direct = 'Yes'
#direct = 'No'

if(direct == 'Yes'){
  
  my_samps = example_samps %>%
    filter(direct_observation == 'Y')
  
  my_df = example_df %>%
    filter(direct_observation == 'Y')

} else {
  
  my_samps = example_samps 
  
  my_df = example_df 
  
}
  
## stomach sample information for taxa and body mass bins (i.e. bin_number) 
# these are collectively identified and abbreviated in the pred_id column 
stomsum = my_samps %>%
  group_by(pred_taxa, bin_number, pred_id) %>%
  summarise(sum_n_stomachs = sum(n_stomachs),
            sum_unique_hauls = length(unique(haul_id)))


# calculate prey biomass at the family-level
gdiet_f = my_df %>%
  # estimate at sample level first
  group_by(sample_id, prey_family) %>%
  summarise (gprey = sum(prey_weight, na.rm = T)) %>%
  pivot_wider(names_from = prey_family, 
              values_from = gprey,
              values_fill = list(gprey = 0)) %>%
  pivot_longer(cols = -sample_id,
               names_to = 'prey_family',
               values_to = 'gprey') %>% 
  left_join(dplyr::select(my_samps, sample_id, pred_id), #bin_number, pred_taxa, data, pred_weight_g), 
                   'sample_id') %>%
  # then take means for pred_id
  group_by(pred_id, prey_family) %>%
  summarise(av_gprey = mean(gprey)) %>%
  mutate(gprop = 100*(av_gprey / sum(av_gprey))) %>%
  filter(av_gprey>0) %>%
  ungroup()

## select predators with sufficient observations
# 30 or more stomachs from 3 or more unique hauls
good = stomsum %>%
  filter(sum_n_stomachs >29 & sum_unique_hauls > 2,
         !is.na(pred_taxa),
         !is.na(bin_number)) %>%
  distinct(pred_id) %>%
  pull()

# identify predators which exploit many different prey taxa 
goodinf = subset(gdiet_f, pred_id %in% good) %>%
  group_by(pred_id) %>%
  summarise(prey_n = length(unique(prey_family)))

# select preds which exploit > median prey families
tt = table(gdiet_f$pred_id)
N = median(goodinf$prey_n)
subocc = subset(gdiet_f, pred_id %in% names(tt[tt >=N])) %>%
  arrange(desc(gprop))
subocc = Reduce(rbind, by(subocc, subocc["pred_id"], head, n=N)) 
#range(table(subocc$pred_id))

# remove rare prey from well-sampled predators - see https://doi.org/10.1111/1365-2664.13662 for justification
gdiet_f = rbind(subset(gdiet_f, !pred_id %in% subocc$pred_id), subocc) 

## pred_id matrices and averages
# occurrence data
ocdf = gdiet_f %>%
  filter(pred_id %in% good) %>%
  mutate(occ = 1) %>%
  dplyr::select(-c(av_gprey, gprop)) %>%
  pivot_wider(names_from = prey_family, 
              values_from = occ,
              values_fill = list(occ = 0)) %>%
  as.data.frame()
rownames(ocdf) = ocdf$pred_id
ocdf$pred_id = NULL
#summary(rowSums(ocdf[,-1]))

# iteratively remove preds and prey with only 1 entry since these hamper dissimilarity analyses
ocdf = ocdf[colSums(ocdf) > 1]
ocdf = subset(ocdf, rowSums(ocdf)>1)
ocdf = ocdf[colSums(ocdf) > 1]
ocdf = subset(ocdf, rowSums(ocdf)>1)
min(colSums(ocdf)) # both should sum >1
min(rowSums(ocdf))

# final list of predators to classify and prey taxa to use
preds = row.names(ocdf)
preys = colnames(ocdf)

# generate final prey biomass matrix to use in the dissimilarity analysis
gdf = my_df %>%
  filter(pred_id %in% preds,
         prey_family %in% preys) %>%
  # estimate at sample level first
  group_by(sample_id, prey_family) %>%
  summarise (gprey = sum(prey_weight, na.rm = T)) %>%
  pivot_wider(names_from = prey_family, 
              values_from = gprey,
              values_fill = list(gprey = 0)) %>%
  pivot_longer(cols = -sample_id,
               names_to = 'prey_family',
               values_to = 'gprey') %>% 
  left_join(dplyr::select(my_samps, sample_id, pred_id), #bin_number, pred_taxa, data, pred_weight_g), 
            'sample_id') %>%
  # then take means for pred_id
  group_by(pred_id, prey_family) %>%
  summarise(av_gprey = mean(gprey)) %>%
  mutate(gprop = 100*(av_gprey / sum(av_gprey))) %>%
  filter(av_gprey>0) %>%
  ungroup() %>%
  dplyr::select(-av_gprey) %>%
  pivot_wider(names_from = prey_family, 
              values_from = gprop,
              values_fill = list(gprop = 0)) %>%
  as.data.frame()
rownames(gdf) = gdf$pred_id
gdf$pred_id = NULL

# calculate proportions of selected main prey functional groups at the stomach sample-level
prey_func_biomass_by_sample = my_df %>% 
  filter(pred_id %in% preds,          
         prey_family %in% preys,
         prey_weight >0) %>%
  # at sample-level
  group_by(sample_id, pred_id, prey_func) %>%
  summarise(sum_g = sum(prey_weight, na.rm = T)) %>%
  mutate(gprop = 100*(sum_g / sum(sum_g))) %>%
  dplyr::select(-sum_g) %>%
  pivot_wider(names_from = prey_func, 
              values_from= gprop,
              values_fill = list(gprop = 0))

# calculate PPMR for each stomach sample
ppmr = my_df %>% 
  left_join(dplyr::select(my_samps, sample_id, pred_weight_g), 'sample_id') %>% 
  mutate(ppmr = pred_weight_g/ind_prey_weight) %>%
  group_by(sample_id) %>% 
  summarise(totb = sum(prey_weight),
            totn = sum(prey_count),
            av_ind_prey_weight = mean(sum(prey_weight)/sum(prey_count)),
            bw_ppmr = weighted.mean(ppmr, prey_weight)) %>%
  ungroup() %>%
  left_join(prey_func_biomass_by_sample, 'sample_id') 

# calculate average % prey functional group for pred_ids  
av_g = prey_func_biomass_by_sample %>%
  group_by(pred_id) %>%
  summarise(fish = mean(fish, na.rm=T),
            zooplankton = mean(zooplankton, na.rm=T),
            benthos = mean(benthos, na.rm=T),
            nekton = mean(nekton, na.rm=T),
            other = mean(other, na.rm=T))
#rowSums(dplyr::select(av_g, fish, zooplankton, benthos, nekton,  other))

# calculate log10 averages of ppmr and prey mass for pred_ids  
av_ppmr = ppmr %>%
  group_by(pred_id) %>%
  filter(pred_id %in% preds) %>%
  summarise(log10_ppmr = log10(mean(bw_ppmr, na.rm=T)),
            log10_prey_mass = log10(mean(av_ind_prey_weight, na.rm=T)))

# generate the trait matrix for dissimilarity analysis
sdf = dplyr::select(av_g, -other) %>%
  left_join(av_ppmr, 'pred_id') %>%
  filter(pred_id %in% preds) %>%
  mutate_at(vars(-1), scales::rescale) %>%
  as.data.frame()
rownames(sdf) = sdf$pred_id
sdf$pred_id = NULL


############################################################################################################
### CLUSTER ANALYSIS

# calculate pairwise dissimilarities (jaccard suitable for occurrence data)
set.seed(1)
distdf1 = vegdist(sdf, na.rm = T, method='bray')
set.seed(10)
distdf2 = vegdist(ocdf, na.rm = T, method='jaccard')
set.seed(100)
distdf3 = vegdist(gdf, na.rm = T, method='bray')

# hclust hierarchical clustering of dissimilarity matrices
set.seed(3)
fitdf1 = hclust(distdf1, method="ward.D2")# 
set.seed(30)
fitdf2 = hclust(distdf2, method="ward.D2")# 
set.seed(300)
fitdf3 = hclust(distdf3, method="ward.D2")# 

# define n guilds:
nguilds = 4 

# predator latitude and longitude centroids 
pred_info = my_samps %>%
  filter(pred_id %in% preds) %>%
  group_by(pred_id) %>%
  summarise(latitude = mean(latitude, na.rm=T),
            longitude = mean(longitude, na.rm=T))

# predator size information
bin_info =  my_samps %>%
  filter(pred_id %in% preds) %>%
  group_by(pred_id, bin_number) %>%
  summarise(min_g = min(pred_weight_g, na.rm=T),
            max_g = max(pred_weight_g, na.rm=T),
            min_cm = min(pred_length_cm, na.rm=T),
            max_cm = max(pred_length_cm, na.rm=T),
            mean_g = mean(pred_weight_g, na.rm=T),
            mean_cm = mean(pred_length_cm, na.rm=T)) %>%
  arrange(bin_number, pred_id)

# feeding guild table with all approaches and supporting information
fguilds = data.frame(pred_id=c(fitdf1$labels, 
                               fitdf2$labels, 
                               fitdf3$labels), 
                     method=c(rep('trait', length(fitdf1$labels)), 
                              rep('occurrence', length(fitdf2$labels)),
                              rep('biomass', length(fitdf3$labels))),
                     #guilds_3=c(cutree(fitdf1, k = 3), cutree(fitdf2, k = 3), cutree(fitdf3, k = 3))) %>%
                     guilds=c(cutree(fitdf1, k = nguilds), cutree(fitdf2, k = nguilds), cutree(fitdf3, k = nguilds))) %>%
  pivot_wider(names_from = method, 
              values_from = guilds,
              values_fill = list(guilds = NA_real_)) %>%
  left_join(dplyr::select(stomsum, bin_number, pred_taxa, pred_id), 'pred_id') %>%
  left_join(av_g, 'pred_id') %>%
  left_join(av_ppmr, 'pred_id') %>%
  mutate(trait = as.factor(paste('Guild_', trait, sep='')),
         occurrence = as.factor(paste('Guild_', occurrence, sep='')),
         biomass = as.factor(paste('Guild_', biomass, sep=''))) %>%
  arrange(trait, occurrence, biomass, pred_taxa, bin_number) %>%
  left_join(pred_info, 'pred_id') %>%
  left_join(bin_info, c('pred_id', 'bin_number'))

# generate latitude and longitude centroids for the entire dataset 
av_coords = fguilds %>%
  summarise(av_lat = mean(latitude),
            av_lon = mean(longitude))

# generate latitude and longitude centroids for each guild
dist_2_centroid = unique(dplyr::select(fguilds, pred_id, pred_taxa, bin_number, occurrence, trait, biomass, latitude, longitude)) %>% 
  pivot_longer(cols = -c(pred_id, pred_taxa, bin_number, latitude, longitude),
               names_to = 'method',
               values_to = 'guild')  %>%
  group_by(method, guild) %>%
  summarise(av_lat = mean(latitude),
            av_lon = mean(longitude)) %>%
  mutate(lat = av_coords$av_lat,
         lon = av_coords$av_lon) %>%
  as.data.frame()

# calculate distances to centroid
dist_2_centroid$km_2_centroid = distCosine(dist_2_centroid[,c('av_lon', 'av_lat')], dist_2_centroid[,c('lon', 'lat')])/1000

# best method = lowest sum of distances to centroid - i.e. preds from different ecosystems 
# share more common guilds rather than occupying unique ecosystem-specific guilds
dist_2_centroid %>%
  group_by(method) %>%
  summarise(sum_km_2_centroid = sum(km_2_centroid))

# world map
world_shp = sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))

# guild-specific data centroids
dat_cen = dist_2_centroid %>%
  distinct(lat, lon) 

# set map limits 
xlims = range(dist_2_centroid$av_lon)
ylims = range(dist_2_centroid$av_lat)

# generate the plot
dist_2_centroid %>% ggplot() +
  geom_sf(data = world_shp, 
          fill = 'black', 
          color = 'black',
          size = 0.01) +
  coord_sf(xlim = c(-75,5), ylim = c(35,60)) +
  labs(x='Latitude', y='Longitude') +
  geom_point(aes(x = av_lon, y = av_lat, colour = method), size=5) + 
  geom_point(data=dat_cen, aes(x = lon, y = lat), colour = 'white', shape = 3, size=5) + 
  guides(fill = guide_colourbar(barwidth = 1)) +
  theme(panel.background = element_rect(fill = 'grey80'),
        panel.border = element_rect(colour='black', fill=NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


############################################################################################################
### test different approaches ability to classifying guilds against one-another and randomly assigned guilds
# re-sample data repeatedly and re-cluster, comparing results

# prey biomass data
gdiet_fc = my_df %>%
  filter(pred_id %in% preds,          
         prey_family %in% preys)  %>% 
  group_by(prey_family, sample_id) %>%
  summarise (gprey = sum(prey_weight, na.rm = T)) %>%
  pivot_wider(names_from = prey_family, 
              values_from = gprey,
              values_fill = list(gprey = 0)) %>%
  ungroup() 

# feeding trait data combined with prey biomass data
s_df = dplyr::select(ppmr, pred_id, sample_id, bw_ppmr, av_ind_prey_weight, 
                     zooplankton, benthos, fish, nekton)  %>%  
  filter(pred_id %in% preds) %>%
  mutate(log10_ppmr = log10(bw_ppmr),
         log10_prey_mass = log10(av_ind_prey_weight+1)) %>%
  dplyr::select(-c(bw_ppmr, av_ind_prey_weight)) %>%
  left_join(gdiet_fc, 'sample_id') 

# identify and remove any 0s and NAs
s_df$rowsum = rowSums(dplyr::select(s_df, -c(sample_id, pred_id, #data,
                                      zooplankton, benthos, fish, nekton, log10_ppmr, log10_prey_mass)))
s_df = s_df %>% 
  filter(!is.na(rowsum),
         rowsum > 0) %>%
  dplyr::select(-rowsum)

# table to store results
samp_fguilds = data.frame(pred_id=NA, data=NA, guilds=NA, method=NA, r_o=NA,
                          zooplankton=NA, benthos=NA, fish=NA, nekton=NA, 
                          log10_ppmr=NA, log10_prey_mass=NA)

ndfs = 100 # number of times we re-sample data (we use n = 1000 in the main study)
nsamps = 30 # number of samples per predator per data re-sample

# data resampling:
for(ro in c('Random', 'Observed')) { # ro = 'Observed' # ro = 'Random'
  
  for(i in 1:ndfs) { # i = 1
    
    # sample data
    set.seed(i)
    sub_samp = s_df %>%
      group_by(pred_id) %>% 
      sample_n(nsamps, replace = T) %>% # 
      mutate(data=i) %>%
      pivot_longer(cols = -c(pred_id, sample_id, data),
                   names_to = 'variable',
                   values_to = 'value')
    
    # traits data
    fglds_props = sub_samp %>%
      filter(variable %in% c('zooplankton', 
                             'benthos', 'fish', 'nekton', 'log10_ppmr', 'log10_prey_mass')) %>%
      group_by(pred_id, variable) %>%
      summarise(av_vals = mean(value)) %>%
      pivot_wider(names_from = variable, 
                  values_from = av_vals,
                  values_fill = list(av_vals = 0))
    
    # biomass data
    df = sub_samp %>%
      filter(!is.na(pred_id)) %>%
      group_by(pred_id, data, variable) %>%
      summarise(value = mean(value, na.rm=T)) %>%
      filter(!is.na(value),
             value>0) %>%
      pivot_wider(names_from = variable, 
                  values_from = value,
                  values_fill = list(value = 0)) %>%
      ungroup() %>%
      dplyr::select(-data) %>%
      as.data.frame()
    
    if(ro == 'Random') {
      
      # randomly order preds
      set.seed(i)
      df$pred_id = sample(df$pred_id)
      
    } else {
      
      df$pred_id = df$pred_id
      
    }
    
    rownames(df) = df$pred_id
    df$pred_id = NULL
    
    # prey taxa occurrence 
    df2 = decostand(dplyr::select(df, -c(fish, zooplankton, nekton, 
                                         benthos, log10_ppmr, log10_prey_mass)), method="pa")
    
    # biomass simple prey functional groups
    df3 = dplyr::select(df, fish, zooplankton, nekton,  
                        benthos, log10_ppmr, log10_prey_mass)
    
    # biomass taxa (remove those simple columns from df)
    df = dplyr::select(df, -c(fish, zooplankton, 
                              nekton, benthos, log10_ppmr, log10_prey_mass))
    
    
    # create dissimilarities
    set.seed(i)
    distdf = vegdist(df, method='bray')
    distdf2 = vegdist(df2, method='jaccard')
    distdf3 = vegdist(df3, method='bray')
    
    # cluster analysis
    set.seed(i)
    fitdf = hclust(distdf, method="ward.D2")
    fitdf2 = hclust(distdf2, method="ward.D2")
    fitdf3 = hclust(distdf3, method="ward.D2")
    
    # store results
    fguild_i = data.frame(pred_id=fitdf$labels, 
                          data=rep(i, length(fitdf$labels)), 
                          guilds=cutree(fitdf, k = nguilds),
                          r_o=ro,
                          method = 'Biomass') %>%
      left_join(fglds_props, 'pred_id')
    
    fguild_i2 = data.frame(pred_id=fitdf2$labels, 
                           data=rep(i, length(fitdf2$labels)), 
                           guilds=cutree(fitdf2, k = nguilds),
                           r_o=ro,
                           method = 'Occurrence') %>%
      left_join(fglds_props, 'pred_id')
    
    fguild_i3 = data.frame(pred_id=fitdf3$labels, 
                           data=rep(i, length(fitdf3$labels)), 
                           guilds=cutree(fitdf3, k = nguilds),
                           r_o=ro,
                           method = 'Trait') %>%
      left_join(fglds_props, 'pred_id')
    
    samp_fguilds = rbind(samp_fguilds, fguild_i, fguild_i2, fguild_i3) %>%
      filter(!is.na(pred_id))
    
  }
  
}



###########  estimate within method guild predator species turnover based on ndfs & nsamps 

res_df = data.frame(guild=NA, dist_cen=NA, method=NA, r_o=NA) 
t_dat = data.frame(gld_dat=NA)

for(ro in c('Random', 'Observed')) { # ro = 'Observed' # ro = 'Random'

  df_ro = subset(samp_fguilds, r_o == ro) 
  
  for(j in unique(samp_fguilds$method)) { # j = 'Trait'
    
    # pred x group and data matrix
    df_j = ungroup(df_ro) %>%
      filter(method == j ) %>%
      # create variable which captures guild and data.frame number 
      # (i.e. identifier for each classification per re-sampled data-frame) 
      mutate(gld_dat = paste(guilds, data, sep='_')) %>%
      group_by(gld_dat) %>%
      reframe(pred_id = unique(pred_id)) %>%
      mutate(occ = 1) %>%
      pivot_wider(names_from = gld_dat, 
                  values_from = occ,
                  values_fill = list(occ = 0)) %>%
      arrange(pred_id) %>% ungroup() %>%
      as.data.frame()
    rownames(df_j) = df_j$pred_id
    
    tdat = as.data.frame(t(df_j[,-1]))
    
    # turnover between guilds
    set.seed(1)
    sim = beta.pair(tdat, index.family = "jaccard")$beta.jac
    set.seed(33)
    sim_fit = hclust(sim, method="ward.D2")
    
    # store results
    glds = data.frame(gld_dat=sim_fit$labels, 
                      gld=cutree(sim_fit, k = nguilds)) 
    
    # append info
    tdat$gld_dat = rownames(tdat)
    tdat = tdat %>%
      left_join(glds, 'gld_dat') %>%
      separate(gld_dat, into = c('group', 'data'), remove = F) %>%
      mutate(method = j,
             r_o = ro)
    
    # 
    set.seed(2)
    tdist = vegdist(tdat[ , -which(names(tdat) %in% c('gld_dat', 'group', 'data', 'gld', 'method', 'r_o'))], method='jaccard')
    set.seed(4)
    tbet = with(tdat, betadisper(tdist, gld, type = "median", bias.adjust=TRUE))
    
    res_df = res_df %>%
      bind_rows(data.frame(guild=tbet$group, dist_cen=tbet$distances, method=j, r_o=ro, data=tdat$data)) %>%
      filter(!is.na(guild))
    
    t_dat = t_dat %>%
      bind_rows(tdat) %>%
      filter(!is.na(gld)) 
    
  }
  
}


## test methods and guild vs one-another and null hypothesis (i.e. random)
# approach with lowest distance to centroid here produces most robust guilds 
mod_df = res_df %>%
  filter(r_o == 'Observed') %>%
  mutate(Guild = as.factor(guild),
         Method = as.factor(method),
         Data = as.factor(data),
         cntrst = as.factor(paste(Method, Guild, sep='_')))

# plot histograms, note sqrt() transformation
mod_df %>%
  ggplot(aes(sqrt(dist_cen))) +
  geom_histogram() +
  facet_wrap(vars(Method))


# biomass_model
b_df = subset(mod_df, Method=='Biomass')
b_mod = aov(sqrt(dist_cen) ~ Guild + Data, b_df) 
b_comp = drop1(b_mod, test = "F")
b_comp_tab = data.frame(Model='Biomass',
                        Predictor=rownames(b_comp),
                        Df=b_comp$Df,
                        AIC = round(b_comp$AIC),
                        F_value=round(b_comp$'F value', 2),
                        p = ifelse(b_comp$`Pr(>F)` <0.001, '<0.001', round(b_comp$`Pr(>F)`, 3)))
# occ_model
o_df = subset(mod_df, Method=='Occurrence')
o_mod = aov(sqrt(dist_cen) ~ Guild + Data, o_df) 
o_comp = drop1(o_mod, test = "F")
o_comp_tab = data.frame(Model='Occurrence',
                        Predictor=rownames(o_comp),
                        Df=o_comp$Df,
                        AIC = round(o_comp$AIC),
                        F_value=round(o_comp$'F value', 2),
                        p = ifelse(o_comp$`Pr(>F)` <0.001, '<0.001', round(o_comp$`Pr(>F)`, 3)))
# trait_model
t_df = subset(mod_df, Method=='Trait')
t_mod = aov(sqrt(dist_cen) ~ Guild + Data, t_df) 
t_comp = drop1(t_mod, test = "F")
t_comp_tab = data.frame(Model='Trait',
                         Predictor=rownames(t_comp),
                         Df=t_comp$Df,
                         AIC = round(t_comp$AIC),
                         F_value=round(t_comp$'F value', 2),
                         p = ifelse(t_comp$`Pr(>F)` <0.001, '<0.001', round(t_comp$`Pr(>F)`, 3)))

#  diagnostic plots
par(mfrow=c(6,2), oma=c(0,0,0,0), mar=c(2,2,2,1))
plot(b_mod); plot(o_mod); plot(t_mod)

# model comparing approaches
a_mod_full = aov(sqrt(dist_cen) ~ Guild + Method + Data, mod_df) 
mod_comp = drop1(a_mod_full, test = "F")
mod_comp_tab = data.frame(Model='Across methods',
                          Predictor=rownames(mod_comp),
                          Df=mod_comp$Df,
                          AIC = round(mod_comp$AIC),
                          F_value=round(mod_comp$'F value', 2),
                          p = ifelse(mod_comp$`Pr(>F)` <0.001, '<0.001', round(mod_comp$`Pr(>F)`, 3)))

# hypothesis: can detect feeding guilds; differences between methods
rbind(b_comp_tab, o_comp_tab, t_comp_tab, mod_comp_tab)

final_a_mod = a_mod_full
#summary(final_a_mod) 

par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(2,2,2,1))
plot(final_a_mod)

# test differences between approaches
set.seed(15)
a_mod_test = summary(glht(final_a_mod, linfct = mcp(Method = "Tukey")))
data.frame(test = names(a_mod_test$test$coefficients),
                       fit = round(a_mod_test$test$coefficients, 3),
                       SE = round(a_mod_test$test$sigma, 3),
                       t = round(a_mod_test$test$tstat, 3),
                       p = ifelse(a_mod_test$test$pvalues[1:3] < 0.001, '0.001', 
                                  round(a_mod_test$test$pvalues[1:3], 3)),
                       row.names = NULL)

# plot differences 
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(2,10,2,1))
plot(glht(final_a_mod, linfct = mcp(Method = "Tukey")))


# plotting data
clust_plot_df = res_df %>%
  filter(r_o == 'Observed') %>%
  mutate(cntrst = as.factor(paste(method, guild, sep='_')))

rand_plot_df = res_df %>%
  filter(r_o == 'Random') %>%
  mutate(cntrst = as.factor(paste(method, guild, sep='_')))

mod_plot_info = as.data.frame(effect('Method', final_a_mod)) %>%
  mutate(dist_cen = fit^2,
         lwr = lower^2,
         upr = upper^2)
mod_plot_means = mod_plot_info %>%
  group_by(Method) %>%
  summarise(av_dist = round(mean(dist_cen), 2)) %>%
  mutate(lbl = paste('mean =', av_dist, sep=' '))

rand_plot_means = rand_plot_df %>%
  group_by(method) %>%
  summarise(av_dist = round(mean(dist_cen), 2)) %>%
  mutate(lbl = paste('mean =', av_dist, sep=' '))

clust_mod_plot = ggplot(clust_plot_df, aes(x = method, y = dist_cen)) +
  geom_point(cex = 1, pch = 1.0, position = position_jitter(w = 0.1, h = 0)) +
  ylim(0, 1) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', 
               colour = "red", width=.5, linewidth=0.75) +
  labs(x='Method', y='Distance to centroid', title='Cluster-based feeding guilds') +
  geom_text(data=mod_plot_means, 
            mapping = aes(x = Method, y = 0.9, label = lbl, colour = 'red')) +
  theme(legend.position="none",
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour='black', fill=NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


rand_mod_plot = ggplot(rand_plot_df, aes(x = method, y = dist_cen)) +
  geom_point(cex = 1, pch = 1.0,position = position_jitter(w = 0.1, h = 0)) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', 
               colour = "red", width=.5, linewidth=0.75) +
  ylim(0, 1) +
  labs(x='Method', y='Distance to centroid', title='Randomly assigned feeding guilds') +
  geom_text(data=rand_plot_means, 
            mapping = aes(x = method, y = 0.9, label = lbl, colour = 'red')) +
  theme(legend.position="none",
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour='black', fill=NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


clust_mod_plot / rand_mod_plot

## plot fguilds - these can be used to help name guilds based on prey functional group contributions

# organise the data  
pltdat = fguilds %>%
  dplyr::select(trait, pred_taxa, bin_number, zooplankton, benthos, fish, 
                mean_cm, log10_prey_mass, log10_ppmr) %>%
  pivot_longer(cols = -c(trait, pred_taxa, bin_number),
               names_to = 'response',
               values_to = 'value') %>%
  rename(Guild = trait) %>%
  mutate(response = recode(response, 'mean_cm'='Predator size (cm)', 'log10_prey_mass'='Log10(prey mass [g])',
                           'log10_ppmr' = 'Log10(PPMR)', 'zooplankton'='Zooplankton', 'benthos'='Benthos', 
                           'fish'='Fish'),
         response = factor(response , levels =c('Predator size (cm)', 'Log10(prey mass [g])', 'Log10(PPMR)',
                                                'Zooplankton', 'Benthos', 
                                                'Fish'))) 

# subset to plot
toplot_main = levels(pltdat$response)[!levels(pltdat$response) %in% c("Other", "Nekton")]
#toplot_supp = levels(pltdat$response)[levels(pltdat$response) %in% c('nekton', 'other')]

# generate plots
S_plot_list = list()
for(res in levels(pltdat$response)){ 
  
  mydata = subset(pltdat, response==res)
  
  if(res %in% c("Predator size (cm)","Log10(prey mass [g])","Log10(PPMR)")) {
    
    p = ggplot(mydata, aes(x=Guild, y=value, fill = Guild)) + 
      geom_point(cex = 3, pch = 21,  position = position_jitter(w = 0.1, h = 0)) +
      stat_summary(fun.data = 'mean_se', geom = 'errorbar', 
                   colour = "black", width=.5, linewidth=0.75) +
      xlab("") + ylab(res) + ggtitle(res) + 
      theme(#legend.position="none", 
        legend.text = element_text(size=16),
        legend.key.size = unit(1.5, 'cm'),
        text = element_text(size=13), 
        axis.text.x=element_blank (),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
    
  } else {
    
    p = ggplot(mydata, aes(x=Guild, y=value, fill = Guild)) + 
      geom_point(cex = 3, pch = 21,  position = position_jitter(w = 0.1, h = 0)) +
      stat_summary(fun.data = 'mean_se', geom = 'errorbar', 
                   colour = "black", width=.5, linewidth=0.75) +
      coord_cartesian(ylim=c(0, 100)) +
      xlab("") + ylab('%') + ggtitle(res) + 
      theme(#legend.position="none", 
        legend.text = element_text(size=16),
        legend.key.size = unit(1.5, 'cm'),
        text = element_text(size=13), 
        axis.text.x=element_blank (),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
    
  }
  
  S_plot_list[[res]] = p
}

do.call(ggarrange, c(S_plot_list[toplot_main], align='h', ncol=2, nrow=3, common.legend = TRUE, legend = "bottom"))  
