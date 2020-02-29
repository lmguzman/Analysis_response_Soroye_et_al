library(tidyverse)
library(raster)


#### shape data for model ######################################################################################

continent <- "BOTH"

# import occupancy data
occu_maps <- readRDS(paste0('model_out/', continent, "/bombus_occumaps.RDS"))
spplist <- names(occu_maps)

# import exposure data

sppcontient <- read.csv('data/bombus_err_obs.csv')

TEIspp_bs <- read_rds("data/TEIspp_bs.RDS")
PEIspp_bs <- read_rds("data/PEIspp_bs.RDS")
TEIspp_delta <- read_rds("data/TEIspp_delta.RDS")
PEIspp_delta <- read_rds("data/PEIspp_delta.RDS")

#name layers
names(TEIspp_bs) <- names(PEIspp_bs) <- names(TEIspp_delta) <- names(PEIspp_delta) <- spplist

# import mean climate data
meantemp_bs <- read_rds("data/meantemp_bs.RDS")
meanpre_bs <- read_rds("data/meanprecip_bs.RDS")
meantemp_delta <- read_rds("data/meantemp_delta.RDS")
meanpre_delta <- read_rds("data/meanprecip_delta.RDS")

# import continent data
continent <- read_rds("data/bombus_continent.RDS")

# continent exclusion table (narrow spp occupancy predictions to continent they live)
#sppcontient <- read.csv(paste0(dir_sppcontinent, "bombus_err_obs.csv")) %>% dplyr::select(species, exclude_from_cont)

# name quadrats
quadID <- 1:ncell(continent)

# desired form: 2d matrix where rows= observations, cols= variables
occudat_mod <- 
  do.call(
    rbind, 
    lapply(spplist, function(x){ #iterate through species (and bind together at the end)
      
      # sum over first and last period (i.e. remove seasons)
      mat_spp <- occu_maps[[x]]
      
      # current probability of occupancy
      pr_current <- mat_spp[[2]]
      
      # change in probability of occupancy at sites where species were ever detected 
      # (not sensible extrapolate outside of any data and creates incredible amount of noise)
      pr_change <- pr_current
      
      values(pr_change) <- ifelse(values(mat_spp[[1]])==1 | values(mat_spp[[2]])==1, (values(mat_spp[[2]]) - values(mat_spp[[1]]))*100, NA )
      
      # spp name
      spp <- rep(x, NROW(getValues(pr_current)))
      
      sppdat <- data.frame(curr_occu= getValues(pr_current),
                           hist_occu= getValues(mat_spp[[1]]),
                           procc_change= getValues(pr_change),
                           species = as.factor(spp),
                           continent= as.factor(getValues(continent)),
                           site= as.factor(quadID),
                           TEI_bs= getValues(subset(TEIspp_bs, x)),
                           TEI_delta= getValues(subset(TEIspp_delta, x)),
                           PEI_bs= getValues(subset(PEIspp_bs, x)),
                           PEI_delta= getValues(subset(PEIspp_delta, x)),
                           avgtemp_bs= getValues(meantemp_bs),
                           avgtemp_delta= getValues(meantemp_delta),
                           avgprecip_bs= getValues(meanpre_bs),
                           avgprecip_delta= getValues(meanpre_delta),
                           
                           # standardize continuous variables
                           sc_TEI_bs= (values(subset(TEIspp_bs, x))-mean(values(TEIspp_bs), na.rm=T))/sd(values(TEIspp_bs), na.rm=T),
                           sc_TEI_delta= (values(subset(TEIspp_delta, x))-mean(values(TEIspp_delta), na.rm=T))/sd(values(TEIspp_delta), na.rm=T),
                           sc_PEI_bs= (values(subset(PEIspp_bs, x))-mean(values(PEIspp_bs), na.rm=T))/sd(values(PEIspp_bs), na.rm=T),
                           sc_PEI_delta= (values(subset(PEIspp_delta, x))-mean(values(PEIspp_delta), na.rm=T))/sd(values(PEIspp_delta), na.rm=T),
                           sc_avgtemp_bs= scale(getValues(meantemp_bs)),
                           sc_avgtemp_delta= scale(getValues(meantemp_delta)),
                           sc_avgprecip_bs= scale(getValues(meanpre_bs)),
                           sc_avgprecip_delta= scale(getValues(meanpre_delta)))
      
      return(sppdat)
      
    } )
  )
str(occudat_mod)
summary(occudat_mod)


##### Figure 3 ########


data_occ_catdelta <- occudat_mod %>% 
  filter(!is.na(procc_change), !is.na(TEI_bs)) %>% 
  left_join(sppcontient, by = "species") %>% 
  # filter(procc_change < -0.01 | procc_change > 0.01) %>% # testing whether near-zero values influenced results (they dont seem to)
  filter((is.na(exclude_from_cont)) | (exclude_from_cont == 1 & continent == 2) | (exclude_from_cont == 2 & continent == 1)) %>% 
  mutate(catTEI_delta= factor(case_when(TEI_delta < -0.005 ~ "colder",
                                        TEI_delta > -0.005 & TEI_delta < 0.005 ~ "same",
                                        TEI_delta > 0.005 ~ "hotter"),
                              levels = c("colder", "same", "hotter")),
         catPEI_delta= factor(case_when(PEI_delta < -0.005 ~ "drier",
                                        PEI_delta > -0.005 & PEI_delta < 0.005 ~ "same",
                                        PEI_delta > 0.005 ~ "wetter"),
                              levels = c("drier", "same", "wetter")),
         continent= factor(ifelse(continent== 1, "North America", "Europe"),
                           levels = c("North America", "Europe")),
         catTEI_bs= factor(case_when(
           TEI_bs < quantile(filter(occudat_mod, !is.na(TEI_bs))$TEI_bs, 0.25) ~ "near cold limit",
           TEI_bs > quantile(filter(occudat_mod, !is.na(TEI_bs))$TEI_bs, 0.25) & 
             TEI_bs < quantile(filter(occudat_mod, !is.na(TEI_bs))$TEI_bs, 0.75) ~ "middle",
           TEI_bs > quantile(filter(occudat_mod, !is.na(TEI_bs))$TEI_bs, 0.75) ~ "near hot limit"),
           levels = c("near cold limit", "middle", "near hot limit")),
         catPEI_bs= factor(case_when(
           PEI_bs < quantile(filter(occudat_mod, !is.na(PEI_bs))$PEI_bs, 0.25) ~ "near dry limit",
           PEI_bs > quantile(filter(occudat_mod, !is.na(PEI_bs))$PEI_bs, 0.25) & 
             PEI_bs < quantile(filter(occudat_mod, !is.na(PEI_bs))$PEI_bs, 0.75) ~ "middle",
           PEI_bs > quantile(filter(occudat_mod, !is.na(PEI_bs))$PEI_bs, 0.75) ~ "near wet limit"),
           levels = c("near dry limit", "middle", "near wet limit")))
summary(data_occ_catdelta)

saveRDS(data_occ_catdelta, paste0(outdir9, "data_occupancy.RDS"))


#### visualize trends
# set category colors
TEIcatcolors <- c("colder" = "#2c7bb6", "same" = "#5e3c99", "hotter" = "#d7191c")
PEIcatcolors <- c("wetter" = "#2c7bb6", "same" = "#5e3c99", "drier" = "#d7191c")
TEIcatcolors_bs <- c("near cold limit" = "#2c7bb6", "middle" = "#5e3c99", "near hot limit" = "#d7191c")
PEIcatcolors_bs <- c("near wet limit" = "#2c7bb6", "middle" = "#5e3c99", "near dry limit" = "#d7191c")

original <- ggplot(data= data_occ_catdelta, aes(x= TEI_delta, y=procc_change, color= catTEI_bs)) + 
  geom_rug(aes(alpha= 0.1, size= 0.1)) + 
  geom_vline(xintercept = 0, linetype= "dashed") +  geom_hline(yintercept = 0, linetype= "dashed") + 
  stat_smooth(method="lm", formula=y~x, size=2, aes(fill= catTEI_bs)) + 
  scale_color_manual(values = TEIcatcolors_bs) + scale_fill_manual(values = TEIcatcolors_bs) +
  labs(x= "Change in thermal position", y= "Change in Occupancy (%)", color= " Baseline \n thermal position") + guides(fill= F, alpha= F, size=F) +
  theme_bw() + theme(text= element_text(size=20)) + facet_wrap(~continent) 



############ NA/ EU separated FIG 3 ##############


continent <- 'NA-constrained'

# import occupancy data
occu_maps <- readRDS(paste0('model_out/', continent, "/bombus_occumaps.RDS"))
spplist <- names(occu_maps)

dataBeeContinent <- read.csv('data/bombus_err_obs.csv')

# import exposure data
names(dataBeeContinent)

layer_rem <- dataBeeContinent[,c("species"), drop = FALSE] %>% 
  arrange(species) %>% 
  mutate(layer = 1:n()) %>% 
  filter(!species %in% spplist) %>% 
  dplyr::select(layer) %>% unlist()

TEIspp_bs <- read_rds("data/TEIspp_bs.RDS")
PEIspp_bs <- read_rds("data/PEIspp_bs.RDS")
TEIspp_delta <- read_rds("data/TEIspp_delta.RDS")
PEIspp_delta <- read_rds("data/PEIspp_delta.RDS")

## remove layer if not using all species
TEIspp_bs <- dropLayer(TEIspp_bs, layer_rem) 
PEIspp_bs <- dropLayer(PEIspp_bs, layer_rem) 
TEIspp_delta <- dropLayer(TEIspp_delta, layer_rem) 
PEIspp_delta <- dropLayer(PEIspp_delta, layer_rem) 

#name layers
names(TEIspp_bs) <- names(PEIspp_bs) <- names(TEIspp_delta) <- names(PEIspp_delta) <- spplist

# import mean climate data
meantemp_bs <- read_rds("data/meantemp_bs.RDS")
meanpre_bs <- read_rds("data/meanprecip_bs.RDS")
meantemp_delta <- read_rds("data/meantemp_delta.RDS")
meanpre_delta <- read_rds("data/meanprecip_delta.RDS")

#import continent data
continent <- read_rds("data/bombus_continent.RDS")

# continent exclusion table (narrow spp occupancy predictions to continent they live)
#sppcontient <- read.csv(paste0(dir_sppcontinent, "bombus_err_obs.csv")) %>% dplyr::select(species, exclude_from_cont)

# name quadrats
quadID <- 1:ncell(continent)

# desired form: 2d matrix where rows= observations, cols= variables
occudat_mod_NA <- 
  do.call(
    rbind, 
    lapply(spplist, function(x){ #iterate through species (and bind together at the end)
      
      # sum over first and last period (i.e. remove seasons)
      mat_spp <- occu_maps[[x]]
      
      # current probability of occupancy
      pr_current <- mat_spp[[2]]
      
      # change in probability of occupancy at sites where species were ever detected 
      # (not sensible extrapolate outside of any data and creates incredible amount of noise)
      pr_change <- pr_current
      
      values(pr_change) <- ifelse(values(mat_spp[[1]])==1 | values(mat_spp[[2]])==1, (values(mat_spp[[2]]) - values(mat_spp[[1]]))*100, NA )
      
      # spp name
      spp <- rep(x, NROW(getValues(pr_current)))
      
      sppdat <- data.frame(curr_occu= getValues(pr_current),
                           hist_occu= getValues(mat_spp[[1]]),
                           procc_change= getValues(pr_change),
                           species = as.factor(spp),
                           continent= as.factor(getValues(continent)),
                           site= as.factor(quadID),
                           TEI_bs= getValues(subset(TEIspp_bs, x)),
                           TEI_delta= getValues(subset(TEIspp_delta, x)),
                           PEI_bs= getValues(subset(PEIspp_bs, x)),
                           PEI_delta= getValues(subset(PEIspp_delta, x)),
                           avgtemp_bs= getValues(meantemp_bs),
                           avgtemp_delta= getValues(meantemp_delta),
                           avgprecip_bs= getValues(meanpre_bs),
                           avgprecip_delta= getValues(meanpre_delta),
                           
                           # standardize continuous variables
                           sc_TEI_bs= (values(subset(TEIspp_bs, x))-mean(values(TEIspp_bs), na.rm=T))/sd(values(TEIspp_bs), na.rm=T),
                           sc_TEI_delta= (values(subset(TEIspp_delta, x))-mean(values(TEIspp_delta), na.rm=T))/sd(values(TEIspp_delta), na.rm=T),
                           sc_PEI_bs= (values(subset(PEIspp_bs, x))-mean(values(PEIspp_bs), na.rm=T))/sd(values(PEIspp_bs), na.rm=T),
                           sc_PEI_delta= (values(subset(PEIspp_delta, x))-mean(values(PEIspp_delta), na.rm=T))/sd(values(PEIspp_delta), na.rm=T),
                           sc_avgtemp_bs= scale(getValues(meantemp_bs)),
                           sc_avgtemp_delta= scale(getValues(meantemp_delta)),
                           sc_avgprecip_bs= scale(getValues(meanpre_bs)),
                           sc_avgprecip_delta= scale(getValues(meanpre_delta)))
      
      return(sppdat)
      
    } )
  )


continent <- 'EU-constrained'


# import occupancy data
occu_maps <- readRDS(paste0('model_out/', continent, "/bombus_occumaps.RDS"))
spplist <- names(occu_maps)

# import exposure data
names(dataBeeContinent)

layer_rem <- dataBeeContinent[,c("species"), drop = FALSE] %>% 
  arrange(species) %>% 
  mutate(layer = 1:n()) %>% 
  filter(!species %in% spplist) %>% 
  dplyr::select(layer) %>% unlist()

TEIspp_bs <- read_rds("data/TEIspp_bs.RDS")
PEIspp_bs <- read_rds("data/PEIspp_bs.RDS")
TEIspp_delta <- read_rds("data/TEIspp_delta.RDS")
PEIspp_delta <- read_rds("data/PEIspp_delta.RDS")

## remove layer if not using all species
TEIspp_bs <- dropLayer(TEIspp_bs, layer_rem) 
PEIspp_bs <- dropLayer(PEIspp_bs, layer_rem) 
TEIspp_delta <- dropLayer(TEIspp_delta, layer_rem) 
PEIspp_delta <- dropLayer(PEIspp_delta, layer_rem) 

#name layers
names(TEIspp_bs) <- names(PEIspp_bs) <- names(TEIspp_delta) <- names(PEIspp_delta) <- spplist

# import mean climate data
meantemp_bs <- read_rds("data/meantemp_bs.RDS")
meanpre_bs <- read_rds("data/meanprecip_bs.RDS")
meantemp_delta <- read_rds("data/meantemp_delta.RDS")
meanpre_delta <- read_rds("data/meanprecip_delta.RDS")

#import continent data
continent <- read_rds("data/bombus_continent.RDS")

# continent exclusion table (narrow spp occupancy predictions to continent they live)
#sppcontient <- read.csv(paste0(dir_sppcontinent, "bombus_err_obs.csv")) %>% dplyr::select(species, exclude_from_cont)

# name quadrats
quadID <- 1:ncell(continent)

# desired form: 2d matrix where rows= observations, cols= variables
occudat_mod_EU <- 
  do.call(
    rbind, 
    lapply(spplist, function(x){ #iterate through species (and bind together at the end)
      
      # sum over first and last period (i.e. remove seasons)
      mat_spp <- occu_maps[[x]]
      
      # current probability of occupancy
      pr_current <- mat_spp[[2]]
      
      # change in probability of occupancy at sites where species were ever detected 
      # (not sensible extrapolate outside of any data and creates incredible amount of noise)
      pr_change <- pr_current
      
      values(pr_change) <- ifelse(values(mat_spp[[1]])==1 | values(mat_spp[[2]])==1, (values(mat_spp[[2]]) - values(mat_spp[[1]]))*100, NA )
      
      # spp name
      spp <- rep(x, NROW(getValues(pr_current)))
      
      sppdat <- data.frame(curr_occu= getValues(pr_current),
                           hist_occu= getValues(mat_spp[[1]]),
                           procc_change= getValues(pr_change),
                           species = as.factor(spp),
                           continent= as.factor(getValues(continent)),
                           site= as.factor(quadID),
                           TEI_bs= getValues(subset(TEIspp_bs, x)),
                           TEI_delta= getValues(subset(TEIspp_delta, x)),
                           PEI_bs= getValues(subset(PEIspp_bs, x)),
                           PEI_delta= getValues(subset(PEIspp_delta, x)),
                           avgtemp_bs= getValues(meantemp_bs),
                           avgtemp_delta= getValues(meantemp_delta),
                           avgprecip_bs= getValues(meanpre_bs),
                           avgprecip_delta= getValues(meanpre_delta),
                           
                           # standardize continuous variables
                           sc_TEI_bs= (values(subset(TEIspp_bs, x))-mean(values(TEIspp_bs), na.rm=T))/sd(values(TEIspp_bs), na.rm=T),
                           sc_TEI_delta= (values(subset(TEIspp_delta, x))-mean(values(TEIspp_delta), na.rm=T))/sd(values(TEIspp_delta), na.rm=T),
                           sc_PEI_bs= (values(subset(PEIspp_bs, x))-mean(values(PEIspp_bs), na.rm=T))/sd(values(PEIspp_bs), na.rm=T),
                           sc_PEI_delta= (values(subset(PEIspp_delta, x))-mean(values(PEIspp_delta), na.rm=T))/sd(values(PEIspp_delta), na.rm=T),
                           sc_avgtemp_bs= scale(getValues(meantemp_bs)),
                           sc_avgtemp_delta= scale(getValues(meantemp_delta)),
                           sc_avgprecip_bs= scale(getValues(meanpre_bs)),
                           sc_avgprecip_delta= scale(getValues(meanpre_delta)))
      
      return(sppdat)
      
    } )
  )



########### figure 3 correct #######


od_eu <- occudat_mod_EU %>% 
  filter(!is.na(procc_change), !is.na(TEI_bs))

od_na <- occudat_mod_NA %>% 
  filter(!is.na(procc_change), !is.na(TEI_bs)) 

occu_dat_binded <- bind_rows(od_eu, od_na) 

data_occ_catdelta_3 <- occu_dat_binded %>% 
  mutate(catTEI_delta= factor(case_when(TEI_delta < -0.005 ~ "colder",
                                        TEI_delta > -0.005 & TEI_delta < 0.005 ~ "same",
                                        TEI_delta > 0.005 ~ "hotter"),
                              levels = c("colder", "same", "hotter")),
         catPEI_delta= factor(case_when(PEI_delta < -0.005 ~ "drier",
                                        PEI_delta > -0.005 & PEI_delta < 0.005 ~ "same",
                                        PEI_delta > 0.005 ~ "wetter"),
                              levels = c("drier", "same", "wetter")),
         continent= factor(ifelse(continent== 1, "North America", "Europe"),
                           levels = c("North America", "Europe")),
         catTEI_bs= factor(case_when(
           TEI_bs < quantile(filter(occudat_mod, !is.na(TEI_bs))$TEI_bs, 0.25) ~ "near cold limit",
           TEI_bs > quantile(filter(occudat_mod, !is.na(TEI_bs))$TEI_bs, 0.25) & 
             TEI_bs < quantile(filter(occudat_mod, !is.na(TEI_bs))$TEI_bs, 0.75) ~ "middle",
           TEI_bs > quantile(filter(occudat_mod, !is.na(TEI_bs))$TEI_bs, 0.75) ~ "near hot limit"),
           levels = c("near cold limit", "middle", "near hot limit")),
         catPEI_bs= factor(case_when(
           PEI_bs < quantile(filter(occudat_mod, !is.na(PEI_bs))$PEI_bs, 0.25) ~ "near dry limit",
           PEI_bs > quantile(filter(occudat_mod, !is.na(PEI_bs))$PEI_bs, 0.25) & 
             PEI_bs < quantile(filter(occudat_mod, !is.na(PEI_bs))$PEI_bs, 0.75) ~ "middle",
           PEI_bs > quantile(filter(occudat_mod, !is.na(PEI_bs))$PEI_bs, 0.75) ~ "near wet limit"),
           levels = c("near dry limit", "middle", "near wet limit")))



#### visualize trends

removing <- ggplot(data= data_occ_catdelta_3, aes(x= TEI_delta, y=procc_change, color= catTEI_bs)) + 
  geom_rug(aes(alpha= 0.1, size= 0.1)) + 
  geom_vline(xintercept = 0, linetype= "dashed") +  geom_hline(yintercept = 0, linetype= "dashed") + 
  stat_smooth(method="lm", formula=y~x, size=2, aes(fill= catTEI_bs)) + 
  scale_color_manual(values = TEIcatcolors_bs) + scale_fill_manual(values = TEIcatcolors_bs) +
  labs(x= "Change in thermal position", y= "Change in Occupancy (%)", color= " Baseline \n thermal position") + guides(fill= F, alpha= F, size=F) +
  theme_bw() + theme(text= element_text(size=20)) + facet_wrap(~continent) 


fig3_comparison <- plot_grid(original, removing, labels = c("A", "B"), ncol = 1)

ggsave(fig3_comparison, filename = "Fig3_comparison.jpeg", width = 10, height = 10)

