library(tidyverse)
library(raster)
library(cowplot)

############# Figure 2 ##############

continent <- 'NA-constrained'

occu_maps <- readRDS(paste0('model_out/', continent,  "/bombus_occumaps.RDS"))
spplist <- names(occu_maps)
continent_map <- read_rds("data/bombus_continent.RDS")

sppcontient <- read.csv('data/bombus_err_obs.csv')

occudat_mod <- 
  do.call(
    rbind, 
    lapply(spplist, function(x){ #iterate through species (and bind together at the end)
      
      # sum over first and last period (i.e. remove seasons)
      mat_spp <- occu_maps[[x]]
      
      # current probability of occupancy
      pr_current <- mat_spp[[2]]
      # spp name
      
      pr_change <- pr_current
      
      spp <- rep(x, NROW(getValues(pr_current)))
      
      values(pr_change) <- ifelse(values(mat_spp[[1]])==1 | values(mat_spp[[2]])==1, (values(mat_spp[[2]]) - values(mat_spp[[1]]))*100, NA )
      
      sppdat <- data.frame(curr_occu= getValues(pr_current),
                           hist_occu= getValues(mat_spp[[1]]),
                           procc_change = getValues(pr_change),
                           species = as.factor(spp),
                           continent= as.factor(getValues(continent_map)))
    }))

f2_NA <- occudat_mod %>% 
  filter(!is.na(procc_change)) %>%
  filter(!(species == "distinguendus" & continent == "1")) %>% 
  group_by(species, continent) %>% 
  summarise(siteoccu_hist= mean(hist_occu, na.rm = T), siteoccu_pres = mean(curr_occu, na.rm = T)) %>% 
  mutate(siteoccu_ch= siteoccu_pres-siteoccu_hist,
         siteoccu_percch= ((siteoccu_pres-siteoccu_hist)/siteoccu_hist)*100) 

ungroup() %>% 
  dplyr::select(siteoccu_percch) %>% 
  summarise(mean = mean(siteoccu_percch), se = sd(siteoccu_percch)/sqrt(n()))

#### EU constrained ####

continent <- 'EU-constrained'

occu_maps <- readRDS(paste0('model_out/', continent,  "/bombus_occumaps.RDS"))
spplist <- names(occu_maps)
continent_map <- read_rds("data/bombus_continent.RDS")

sppcontient <- read.csv('data/bombus_err_obs.csv')

occudat_mod <- 
  do.call(
    rbind, 
    lapply(spplist, function(x){ #iterate through species (and bind together at the end)
      
      # sum over first and last period (i.e. remove seasons)
      mat_spp <- occu_maps[[x]]
      
      # current probability of occupancy
      pr_current <- mat_spp[[2]]
      # spp name
      
      pr_change <- pr_current
      
      spp <- rep(x, NROW(getValues(pr_current)))
      
      values(pr_change) <- ifelse(values(mat_spp[[1]])==1 | values(mat_spp[[2]])==1, (values(mat_spp[[2]]) - values(mat_spp[[1]]))*100, NA )
      
      sppdat <- data.frame(curr_occu= getValues(pr_current),
                           hist_occu= getValues(mat_spp[[1]]),
                           procc_change = getValues(pr_change),
                           species = as.factor(spp),
                           continent= as.factor(getValues(continent_map)))
    }))


f2_EU <- occudat_mod %>% 
  filter(!is.na(procc_change)) %>%
  group_by(species, continent) %>% 
  summarise(siteoccu_hist= mean(hist_occu, na.rm = T), siteoccu_pres = mean(curr_occu, na.rm = T)) %>% 
  mutate(siteoccu_ch= siteoccu_pres-siteoccu_hist,
         siteoccu_percch= ((siteoccu_pres-siteoccu_hist)/siteoccu_hist)*100)

f2_EU %>% ungroup() %>% 
  dplyr::select(siteoccu_percch) %>% 
  summarise(mean = mean(siteoccu_percch), se = sd(siteoccu_percch)/sqrt(n()))


#### BOTH #####
  
continent <- 'BOTH'

occu_maps <- readRDS(paste0('model_out/', continent,  "/bombus_occumaps.RDS"))
spplist <- names(occu_maps)
continent_map <- read_rds("data/bombus_continent.RDS")

sppcontient <- read.csv('data/bombus_err_obs.csv')

occudat_mod <- 
  do.call(
    rbind, 
    lapply(spplist, function(x){ #iterate through species (and bind together at the end)
      
      # sum over first and last period (i.e. remove seasons)
      mat_spp <- occu_maps[[x]]
      
      # current probability of occupancy
      pr_current <- mat_spp[[2]]
      # spp name
      
      pr_change <- pr_current
      
      spp <- rep(x, NROW(getValues(pr_current)))
      
      values(pr_change) <- ifelse(values(mat_spp[[1]])==1 | values(mat_spp[[2]])==1, (values(mat_spp[[2]]) - values(mat_spp[[1]]))*100, NA )
      
      sppdat <- data.frame(curr_occu= getValues(pr_current),
                           hist_occu= getValues(mat_spp[[1]]),
                           procc_change = getValues(pr_change),
                           species = as.factor(spp),
                           continent= as.factor(getValues(continent_map)))
    }))

f2_both <- occudat_mod %>% 
  filter(!is.na(procc_change)) %>% 
  left_join(sppcontient, by = "species") %>% 
  filter((is.na(exclude_from_cont)) | (exclude_from_cont == 1 & continent == 2) | (exclude_from_cont == 2 & continent == 1)) %>% 
  group_by(species, continent) %>% 
  summarise(siteoccu_hist= mean(hist_occu, na.rm = T), siteoccu_pres = mean(curr_occu, na.rm = T)) %>% 
  mutate(siteoccu_ch= siteoccu_pres-siteoccu_hist,
         siteoccu_percch= ((siteoccu_pres-siteoccu_hist)/siteoccu_hist)*100)

f2_both %>%
  filter(!(species == "distinguendus" & continent == "1")) %>% 
  group_by(continent) %>% 
  summarise(siteoccu_percchSE= sd(siteoccu_percch)/sqrt(length(siteoccu_percch)), mean = mean(siteoccu_percch))


#### Figure 1 -  response ######

f3_both <- f2_both %>% 
  filter(!(species == "distinguendus" & continent == "1")) %>% 
  mutate(Analysis ="Original\n(with known\nabsences)")

f3_na <- f2_NA %>% 
  mutate(Analysis ="Removing\nknown\nabsences")

f3_eu <- f2_EU  %>% 
  mutate(Analysis ="Removing\nknown\nabsences")

analysis_fig2 <- bind_rows(f3_both, f3_na, f3_eu) %>% 
  mutate(continent = ifelse(continent == 1, "North America", "Europe")) %>%
  mutate(continent = factor(continent, levels = c("North America", "Europe"))) %>% 
  ggplot(aes(x = siteoccu_percch, fill = Analysis)) + geom_histogram(alpha = 0.7, position="identity") + facet_wrap(~continent) +
  xlab("Change in probability of site occupancy \n relative to historical site baseline (%)") + ylab("Number of species") +
  theme_cowplot() + 
  theme(strip.background = element_blank(),
        legend.position = 'bottom', 
        strip.text = element_text(size = 25),
        axis.title  = element_text(size = 25), 
        legend.text = element_text(size = 25)) + geom_vline(xintercept = 0, colour = 'red', size = 1) + 
  scale_fill_manual(name = "", values = c("grey","black")) 

analysis_violin <- bind_rows(f3_both, f3_na, f3_eu) %>% 
  mutate(continent = ifelse(continent == 1, "North America", "Europe")) %>%
  mutate(continent = factor(continent, levels = c("North America", "Europe"))) %>% 
  ggplot(aes(x = Analysis, y = siteoccu_percch)) + facet_wrap(~continent) + geom_violin() + theme_cowplot() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 25),
        axis.title  = element_text(size = 25),
        axis.text.x = element_text(size = 20)) + ylab("Change in probability of site occupancy \n relative to historical site baseline (%)") +
  xlab("") + geom_hline(yintercept = 0, colour = 'red', size = 1)

two_plots <- plot_grid(analysis_fig2, analysis_violin, labels = c("A","B"), label_size = 25)

ggsave(analysis_fig2, filename = "Fig_2_re_done.jpeg", width = 10, height = 6)

ggsave(two_plots, filename = "Fig_2_re_done_violin.jpeg", width = 15, height = 9)

ggsave(analysis_violin, filename = "violin.jpeg", width = 9, height = 7)



