## occupancy models to correct imperfect detection in bumblebee data, and model CEI/avg clim to ext/col

library(parallel)
library(raster)
library(R2jags)
library(R2WinBUGS)


## INCLUDE NUMBER OF DETECTIONS AS A FIXED EFFECT ON EFFORT

## JAGS model 
model.jags <- function() {
  for (k in 1:nyears){ ## repeat for each year (with distinct priors)
    
    ## specify priors
    alpha_psi[k] ~ dnorm(0, 0.01)
    alpha_p[k]   ~ dnorm(0, 0.01) 
    
    ## covar effects
    b_lp_samp[k] ~ dnorm(0, 0.1)
    
    ## Ecological submodel: Define state conditional on parameters 
    for (i in 1:nsites) { 
      
      z[i,k] ~ dbern(psi[i,k]) ## true occurrence z at site i 
      psi[i,k] <- 1 / (1 + exp(-lpsi.lim[i,k])) 
      lpsi.lim[i,k] <- min(999, max(-999, lpsi[i,k])) ## stabilize logit to avoid numerical under or overflow
      lpsi[i,k] <- alpha_psi[k]
      
      for (j in 1:nsurveys) {   ## Observation model
        
        y[i,j,k] ~ dbern(mu_p[i,j,k]) ## Detection-nondetection at i and j
        
        mu_p[i,j,k] <- z[i,k] * p[i,j,k] 
        p[i,j,k] <- 1 / (1 + exp(-lp.lim[i,j,k])) 
        lp.lim[i,j,k] <- min(999, max(-999, lp[i,j,k])) 
        lp[i,j,k] <- alpha_p[k] + b_lp_samp[k] * sampmat[i,j,k]
        
      } ## j
    } ## i
    ## extra tings
    ## Derived quantities 
    n_occ[k] <- sum(z[,k]) ## Number of occupied sites
    mean_p[k] <- exp(alpha_p[k]) / (1 + exp(alpha_p[k])) ## Average detection
  } ## k
}

## run for Europe, North America, or all together
run.case <- function(case) {

  ## run a single species
  run.for.single.spp <- function(sp) {
    cat(sprintf('Running %s in %s\n', sp, case))

    ## bundle data 
    win.data <- dataOccu[[sp]]
    ##
    ## subset if only running one continent
    if(case!='BOTH') {

      win.data$y       <- win.data$y[win.data$continent==cc,,]
      win.data$sampmat <- win.data$sampmat[win.data$continent==cc,,]
      win.data$quadID  <- win.data$quadID[win.data$continent==cc]

      ## constrain to sites with at least once detection
      if(case %in% c('NA-constrained', 'EU-constrained')) {

        sites.keep <- which(apply(win.data$y, 1, sum)>1)
        win.data$y       <- win.data$y[sites.keep,,]
        win.data$sampmat <- win.data$sampmat[sites.keep,,]
        win.data$quadID  <- win.data$quadID[sites.keep]
      }
      
      win.data$nsites <- dim(win.data$y)[1]
    }

    ## break out if no sites ...
    if(win.data$nsites==0)
      return(NA)
    
    ## Initial values 
    Zst <- apply(win.data$y, c(1, 3), max)
    inits <- function(){
      list(z=Zst)
    }
    
    ## Parameters monitored 
    params <- c('alpha_psi',
                'alpha_p',
                'b_lp_samp',
                'n_occ',
                'mean_p' ,
                'z')

    sppout_cei <- jags(data=win.data, 
                       inits=inits, 
                       parameters.to.save=params,
                       model.file=model.jags,
                       n.chains=3,
                       n.thin=nt,
                       n.iter=ni,
                       n.burnin=nb)

    fn <- sprintf('%s_ss.RData', sp)
    save(win.data, sppout_cei,
         file=file.path('model_out', case, 'occufits', fn))

    NULL
  }
  
  ## limit only to species on the relevant continent
  ##
  ## North America
  ## continent
  if(case %in% c('NA', 'NA-constrained')) {
    cc <- 0
    keep <-
      is.na(dataBeeContinent[,'exclude_from_cont']) |
      dataBeeContinent[,'exclude_from_cont']==2
    dataOccu <- dataOccu[keep]
  }
  ##
  ## Europe
  if(case %in% c('EU', 'EU-constrained')) {
    cc <- 1
    keep <-
      is.na(dataBeeContinent[,'exclude_from_cont']) |
      dataBeeContinent[,'exclude_from_cont']==1
    dataOccu <- dataOccu[keep]
  }

  ## run models
  if(ncores==1)
    lapply(names(dataOccu), run.for.single.spp)

  if(ncores>1)
    mclapply(names(dataOccu), run.for.single.spp,
             mc.cores=ncores,
             mc.preschedule=FALSE)
  NULL
}

## import data
##
## dataOccu <- readRDS('data/dataOccu.RDS')
## save(dataOccu, file='data/dataOccu.RData')
load('data/dataOccu.RData', verbose=TRUE)
##
dataBeeContinent <- read.csv('data/bombus_err_obs.csv')
## re-order the latter
dataBeeContinent <-
  dataBeeContinent[match(names(dataOccu), dataBeeContinent[,'species']),]

## MCMC settings 
ni <- 1e5
nt <- 1e3
nb <- 5e4
ncores <- 20

### run for each fo the continent options

##### removes all sites where species have never been observed
run.case(case='EU-constrained')
run.case(case='NA-constrained')

### case BOTH re-creates the original results
run.case(case='BOTH')



