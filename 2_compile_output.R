library(tidyverse)
library(parallel)
library(raster)
library(R2jags)
library(R2WinBUGS)

### check Rhats####

### run for each fo the continent options

continent <- 'BOTH' ###BOTH recreates original results 

continent <- 'NA-constrained'  ##### removes all sites where species have never been observed
continent <- 'EU-constrained' ##### removes all sites where species have never been observed

dir <- file.path('model_out', continent, 'occufits')

completed <- list.files(dir)
  
check_rhat <- function(file){
  load(sprintf('%s/%s', dir, file))
  summ <- sppout_cei$BUGSoutput$summary
  vars <- c(grep('alpha_p', rownames(summ)),
            grep('b_lp_samp', rownames(summ)))
  rhat <- summ[vars,'Rhat']
  rhat
}
rhat_out <- sapply(completed, check_rhat)
t(rhat_out) %>% 
  as.data.frame() %>% 
  rownames_to_column('species') %>% 
  summarise_at(vars(`alpha_p[1]`:`b_lp_samp`), max)

t(rhat_out) %>% 
  as.data.frame() %>% 
  rownames_to_column('species') %>% 
  filter(`alpha_psi[1]`>1.1 | `alpha_psi[2]` > 1.1)

## plot chains for worst performing parameter
worst <- which(rhat_out==max(rhat_out),arr.ind=TRUE)
load(sprintf('%s/%s', dir, completed[worst[,2]]))
sims.arr <- sppout_cei$BUGSoutput$sims.array
matplot(sims.arr[,,rownames(worst)], type='l', las=1,
        xlab='Iteration', ylab=rownames(worst))
## for NA (which have run 10^5), this is very ugly... and needs to be
## run longer




###### build the rasters based on the summary ####

# import data

### run for each fo the continent options

continent <- 'NA-constrained'
continent <- 'EU-constrained'

continent <- 'BOTH'



dir <- file.path('model_out', continent, 'occufits')

dataOccu <- read_rds('data/dataOccu.RDS')

completed <- list.files(dir)

map <- read_rds('data/bombus_continent.RDS')
values(map) <- 1:ncell(map)

outdir8 <- "model_out"

summary_results <- function(file){
  
  #sppout_cei <- readRDS(sprintf('%s/%s', dir, file))
  
  load(paste0(dir,"/", file), verbose=T)
  ## extract occu probability and save
  
  sp_name <- gsub('_ss.RData', '', file)
  
  quad_id_df <- data.frame(rowid = 1:length(win.data$quadID), quadID= win.data$quadID)
  
  ## extract occu probability and save
  spp_occu_p1 <- as.data.frame(sppout_cei$BUGSoutput$summary) %>% rownames_to_column("coef") %>% 
    separate(coef, c("coef", "period"), sep="\\[", extra="merge", fill="right") %>% filter(coef=="z") %>% 
    separate(period, c("rowid", "period"), sep=",") %>% 
    mutate(period= as.numeric(gsub(pattern="\\]", x=period, "")),
           rowid= as.numeric(rowid)) %>% 
    filter(period == 1) %>% 
    left_join(quad_id_df) %>% 
    dplyr::select(-'2.5%', -'25%', -'50%', -'75%', -'97.5%', -'Rhat', -'n.eff', -'coef', -rowid, -period)    
  
  spp_occu_p2 <- as.data.frame(sppout_cei$BUGSoutput$summary) %>% rownames_to_column("coef") %>% 
    separate(coef, c("coef", "period"), sep="\\[", extra="merge", fill="right") %>% filter(coef=="z") %>% 
    separate(period, c("rowid", "period"), sep=",") %>% 
    mutate(period= as.numeric(gsub(pattern="\\]", x=period, "")),
           rowid= as.numeric(rowid)) %>% 
    filter(period == 2) %>% 
    left_join(quad_id_df) %>% 
    dplyr::select(-'2.5%', -'25%', -'50%', -'75%', -'97.5%', -'Rhat', -'n.eff', -'coef', -rowid, -period)
  
  spp_occu <- spp_occu_p1 %>% left_join(spp_occu_p2, by= "quadID", suffix= c("_p1","_p2"))
  
  ## set raster templates
  rasmean_p1 <- rasmean_p2 <- rassd_p1 <- rassd_p2 <- map
  
  ##set raster values
  values(rasmean_p1) <- unlist(dplyr::select(left_join(data.frame(quadID = values(map)), spp_occu, by = "quadID"), mean_p1))
  values(rasmean_p2) <- unlist(dplyr::select(left_join(data.frame(quadID = values(map)), spp_occu, by = "quadID"), mean_p2))
  
  saveRDS(rasmean_p1, paste0('model_out/', continent, "/occu_p1/", sp_name, ".RDS"))
  writeRaster(rasmean_p1, paste0('model_out/', continent, "/occu_p1/", sp_name), format= "GTiff", overwrite= TRUE)
  
  saveRDS(rasmean_p2, paste0('model_out/', continent, "/occu_p2/", sp_name, ".RDS"))
  writeRaster(rasmean_p2, paste0('model_out/', continent,"/occu_p2/", sp_name), format= "GTiff", overwrite= TRUE)
}

sapply(completed, summary_results)



completed <- list.files(paste0('model_out/', continent, '/occu_p1'))
completed <- completed[grep('RDS', completed)]

occu_maps <- lapply(completed, function(x){

  rasmean_p1 <- readRDS(file.path('model_out',continent, "occu_p1", x))
  rasmean_p2 <- readRDS(file.path('model_out',continent, "occu_p2", x))
  
  spp_occumap <- stack(rasmean_p1, rasmean_p2)
  
  return(spp_occumap)
  
})
names(occu_maps) <- gsub('\\.RDS', '', completed)

saveRDS(occu_maps, paste0('model_out/', continent, '/bombus_occumaps.RDS'))
