
## change 'fast' to TRUE for a quick run,
## to see if updates to the model run correctly
fast <- FALSE

if(fast) {
    niter <- 500
    nburn <- 0
} else {
    niter <- 600000
    nburn <- 100000
}    

nChains <- 3


## this is for my testing, to install the older (0.5-1) version
## of nimble, to make sure it runs correctly (and identically)
## using either version.
## just leave this commented out.
##
##remove.packages('nimble')
##install.packages('nimble', repos = 'http://r-nimble.org', type = 'source')


## NEW
## next line just loads the NIMBLE library
## it's different because I load it different when running on a remote cluster
if(Sys.info()['nodename'] == 'gandalf') library(nimble, lib.loc = '~/Documents/') else library(nimble)


##library(plyr)   ## NEW not needed ???
library(VGAM)
library(coda)


## NEW
## input custom samplers and distributions for your model
source('defs.R')


## Import data
##mdata <- read.csv('mortalities.30km.dword.1996.2013.csv')  ## original
mdata <- read.csv('mortalities.simulated.April.pDet.const.csv')  ## DT adding to try??
## FL regions
regions <- c('ATL', 'USJ', 'NW', 'SW')
nRegions <- length(regions)
## Winter habitat quality
habitats <- c('Low', 'Medium', 'High')
nHabitats <- length(habitats)
## Causes of death
causes <- c('watercraft', 'wcs', 'debris', 'cold', 'redtide', 'other')
nCauses <- length(causes)
full.causes <- c(causes, 'undetermined')
## Reorder with red tide last for ease of programming USJ having no red tide
#causes2 <- c('watercraft', 'wcs', 'debris', 'cold', 'other', 'redtide')
#full.causes2 <- c(causes2, 'undetermined')
## Age classes
classes = c('Calves', 'Subadults', 'Adults')
##classes2 <- c('Calves1', 'Calves2', 'Subadults', 'Adults')  ## NEW nowhere used
nClasses <- length(classes)
severities <- c('Normal', 'Cold', 'Severe')
nSeverities <- length(severities)
intensities <- c("Baseline", "Moderate", "Intense")
nIntensities <- length(intensities)

adult.base.mort <- 1 - c(0.972883, 0.978949, 0.977913, 0.977058)
calf1.mort.ratio <- 0.190/0.031
calf2.mort.ratio <- 0.085/0.031
base.mort <- array(c(1 - sqrt((1 - calf1.mort.ratio * adult.base.mort) * (1 - calf2.mort.ratio * adult.base.mort)), adult.base.mort, adult.base.mort),
                   dim = c(nRegions, nClasses),
                   dimnames = list(regions, classes))
##
data <- transform(mdata,
                  area    = factor(area,    levels = regions),
                  age     = factor(age,     levels = classes),
                  habitat = factor(habitat, levels = habitats))
##
## Constants
## CALF <- 1       ## NEW not used anywhere
## SUBADULT <- 2   ##
## ADULT <- 3      ##
STARTYEAR <- 1996
ENDYEAR <- 2013
years <- STARTYEAR:ENDYEAR
nYears <- ENDYEAR - STARTYEAR + 1
REDTIDE <- 5
##
## Moderate and intense red tide years
mod_tide_years <- c(2002, 2003, 2005, 2006, 2012)
int_tide_years <- c(1996, 2013)
modTideYears <- intTideYears <- rep(0, nYears)
tideYears <- rep(1, nYears)
names(modTideYears) <- names(intTideYears) <- names(tideYears) <- STARTYEAR:ENDYEAR
modTideYears[as.character(mod_tide_years)] <- 1
intTideYears[as.character(int_tide_years)] <- 1
tideYears[as.character(mod_tide_years)] <- 2
tideYears[as.character(int_tide_years)] <- 3

coldDesignations <- matrix(
    c("Cold", rep("Normal", 4), "Cold", "Normal", "Cold", rep("Normal", 6), "Severe", "Cold", rep("Normal", 2),
      "Cold", rep("Normal", 4), "Cold", rep("Normal", 7), "Cold", "Severe", "Cold", rep("Normal", 2),
      "Cold", "Normal", "Cold", "Normal", "Normal", "Cold", "Normal", "Severe", rep("Normal", 6), "Severe", "Cold", rep("Normal", 2),
      "Cold", rep("Normal", 4), "Cold", rep("Normal", 8), "Severe", rep("Normal", 3)),
    nYears, nRegions,
    dimnames = list(STARTYEAR:ENDYEAR, regions))
##coldDesignations
coldYears <- matrix(as.integer(factor(coldDesignations, levels = severities)), 
                    nYears, nRegions,
                    dimnames = list(STARTYEAR:ENDYEAR, regions))

wcsProtectYear <- 2001
## Uninformative priors for proportions and pDet
prior1 <- prior2 <- prior3 <- array(1, c(nClasses, nRegions, nCauses),
                                    dimnames = list(classes, regions, causes))
##
## I haven't figured out how to fix this at zero in NIMBLE, but giving it a prior close to zero might work okay
prior1[,'USJ','redtide'] <- prior2[,'USJ','redtide'] <- 0.001
prior3[,'USJ','redtide'] <- 100
##
## Move from data frame to array
data.array <- array(NA, c(nYears, nClasses, nRegions, nHabitats, nCauses + 1))
dimnames(data.array) <- list(STARTYEAR:ENDYEAR, classes, regions, habitats, full.causes)
for (class in classes) {
    for (region in regions) {
        for (qual in habitats) 
            data.array[ , class, region, qual, ] <- as.matrix(subset(data, age==class & area==region & habitat == qual)[,full.causes])
    }
}

##########################################
##### Full model #########################
##########################################

fraction.code.calf <- nimbleCode({
    ## Red tide effect factors and additional mortality
    tide_mort[1] <- 0
    tide_mort[2] ~ dunif(0, 1-baseMort[SW])
    tide_mort[3] ~ dunif(0, 1-baseMort[SW])
    ## I originally defined the UME factors defined separately (see commented
    ## code).  I removed this to simplify the structure of the model slightly, but
    ## it didn't seem to help anything.
    ##   tideFactor[1] <- 0
    ##   tideFactor[2] <- tide_mort[2] / baseMort[SW]
    ##   tideFactor[3] <- tide_mort[3] / baseMort[SW]
    ## 
    ## Cold effect additional mortality
    cold_mort[1, 1] ~ dunif(0, 1-baseMort[1]) # Low normal
    cold_mort[2, 1] ~ dunif(0, 1-baseMort[1]) # Medium normal
    cold_mort[3, 1] <- 0                      # High normal
    cold_mort[1, 2] ~ dunif(0, 1-baseMort[1]) # Low cold
    cold_mort[2, 2] ~ dunif(0, 1-baseMort[1]) # Medium cold
    cold_mort[3, 2] <- 0                      # High cold
    cold_mort[1, 3] ~ dunif(0, 1-baseMort[1]) # Low severe
    cold_mort[2, 3] ~ dunif(0, 1-baseMort[1]) # Medium severe
    cold_mort[3, 3] ~ dunif(0, 1-baseMort[1]) # High severe
    ##
    ## Loop over regions
    for (area in 1:nRegions) {
        ##
        ## Proportions of mortality for region area
        for (cause in 1:nCauses) {
            pi0[cause, area] ~ dgamma(prior1[area,cause], 1.0)
            pDet[cause, area] ~ dbeta(prior2[area,cause], prior3[area,cause])
        }
        ## Total
        pi[1:nCauses,area] <- pi0[1:nCauses,area] / sum(pi0[1:nCauses,area])
        ## Loop over years
        for (year in 1:nYears) {
            ## REMOVED pDet2 from here
            ##pDet2[year,area,1:nCauses] <- (pDet[1:nCauses, area] + (area==SW) * (tideYears[year] > 1) * tideVector[1:nCauses] * tide_mort[tideYears[year]] / baseMort[SW]) /
            ##    (1 + (area==SW) * (tideYears[year] > 1) * tideVector[1:nCauses] * tide_mort[tideYears[year]] / baseMort[SW])
          ## Loop over habitat qualities
          for (habitat in 1:nHabitats) {
            theta[year,area,habitat,1:nCauses] <- (pi[1:nCauses, area] + (area==SW) * (tideYears[year] > 1) * tideVector[1:nCauses] * tide_mort[tideYears[year]] / baseMort[SW] +
                                                     (1 - (habitat == 3) * (coldYears[year, area] < 3)) * coldVector[1:nCauses] * cold_mort[habitat, coldYears[year, area]] / baseMort[area]) /
              (1 + (area==SW) * (tideYears[year] > 1) * tide_mort[tideYears[year]] / baseMort[SW] + 
                 (1 - (habitat == 3) * (coldYears[year, area] < 3)) * cold_mort[habitat, coldYears[year, area]] / baseMort[area])
            ## Total carcasses due to each cause
            X[year,area,habitat,1:nCauses] ~ dmulti(prob = theta[year,area,habitat,1:nCauses],
                                                    size = totals[year,area,habitat])
            ## REMOVED data1 declaration from here
            ##for (cause in 1:nCauses)
            ##    data1[year,area,habitat,cause] ~ dbin(prob = pDet2[year,area,cause],
            ##                                          size = X[year,area,habitat,cause])
          }
      }
    }
    ##
    ## area=1,2,3, only use pDet, gives conjugate sampling
    for (area in 1:(nRegions-1)) {
        for (year in 1:nYears) {
            for (habitat in 1:nHabitats) {
                for (cause in 1:nCauses)
                    data1[year,area,habitat,cause] ~ dbin(prob = pDet[cause,area], size = X[year,area,habitat,cause])
            }
        }
    }
    ##
    ## area=4=SW, use pDet2, results in RW sampling
    for (year in 1:nYears) {
        pDet2[year,1:nCauses] <- (pDet[1:nCauses, 4] + (tideYears[year] > 1) * tideVector[1:nCauses] * tide_mort[tideYears[year]] / baseMort[SW]) /
            (1 + (tideYears[year] > 1) * tideVector[1:nCauses] * tide_mort[tideYears[year]] / baseMort[SW])
        ## Loop over habitat qualities
        for (habitat in 1:nHabitats) {
            for (cause in 1:nCauses)
                data1[year,4,habitat,cause] ~ dbin(prob = pDet2[year,cause], size = X[year,4,habitat,cause])
        }
    }
})

constants.fraction.calf <- list(
    data1 = data.array[ , 'Calves', , , 1:nCauses],
    nCauses = nCauses,
    nRegions = nRegions, 
    nYears = nYears,
    nSeverities = nSeverities,
    nHabitats = nHabitats,
    SW = 4,
    tideVector = c(0, 0, 0, 0, 1, 0),
    coldVector = c(0, 0, 0, 1, 0, 0),
    prior1 = prior1['Calves', , ], 
    prior2 = prior2['Calves', , ], 
    prior3 = prior3['Calves', , ], 
    tideYears = tideYears,
    coldYears = coldYears,
    baseMort = as.vector(base.mort[,'Calves'])
)
#     undet = data.array[ , 'Calves', , , 1 + nCauses] # Not needed for this version
# )
##
data.fraction.calf <- list(totals = apply(data.array[ , 'Calves', , , ], 1:3, sum)) 

inits.calf <- function(data.array) {
  pi <- p <- pDet <- matrix(0, nCauses, nRegions, dimnames = list(causes, regions))
  U <- array(NA, c(nYears, nRegions, nHabitats, nCauses))
  dimnames(U) <- list(STARTYEAR:ENDYEAR, regions, habitats, causes)
  baseMort <- base.mort[,'Calves']
  for (region in regions) {
    if (region=='SW') {
      p[ , region] <- rdiric(1, rep(1, nCauses))
      for (habitat in habitats) {
        for (yr in STARTYEAR:ENDYEAR) {
          year <- as.character(yr)
          U[year, region, habitat, ] <- rmultinom(1, data.array[year,'Calves',
                                                                region, habitat, nCauses+1], p[ , region])
        }
      }
      pi[,region] <- rdiric(1, apply(data.array[(1 - modTideYears) * 
                                                  (1 - intTideYears) * coldYears[,region]==1,
                                                'Calves',region, , 1:nCauses], 3, sum)+1)
      for (cause in causes) 
        pDet[cause, region] <- rbeta(1, sum(data.array[(1 - modTideYears) * 
                                                         (1 - intTideYears),'Calves', region, , cause])+1,  
                                     sum(U[(1 - modTideYears) * (1 - intTideYears),region, , cause])+1)
    }
    else if (region=='USJ') {
      p[-REDTIDE, region] <- rdiric(1, rep(1, nCauses - 1))
      p[REDTIDE, region] <- 1e-100
      p[REDTIDE - 1, region] <- p[REDTIDE - 1, region] - 1e-100
      for (habitat in habitats) {
        for (yr in STARTYEAR:ENDYEAR) {
          year <- as.character(yr)
          U[year, region, habitat, -REDTIDE] <- rmultinom(1, data.array[year,'Calves',
                                                                        region, habitat, nCauses+1], p[-REDTIDE, region])
          U[year, region, habitat, REDTIDE] <- 0
        }
      }
      pi[-REDTIDE,region] <- rdiric(1, apply(data.array[coldYears[,region]==1,'Calves',
                                                        region, , -c(REDTIDE, nCauses + 1)], 3, sum)+1)
      pi[REDTIDE, region] <- 1e-100
      pi[REDTIDE - 1, region] <- pi[REDTIDE - 1, region] - 1e-100
      for (cause in causes) 
        pDet[cause, region] <- rbeta(1, sum(data.array[ ,'Calves', region, , cause])+1,  
                                     sum(U[,region, , cause])+1)
      pDet[REDTIDE, region] <- 1e-100
    } else {
      p[,region] <- rdiric(1, rep(1, nCauses))
      for (habitat in habitats) {
        for (yr in STARTYEAR:ENDYEAR) {
          year <- as.character(yr)
          U[year, region, habitat, ] <- rmultinom(1, data.array[year, 'Calves',
                                                                region, habitat, nCauses+1], p[ , region])
        }
      }
      pi[,region] <- rdiric(1, apply(data.array[coldYears[,region]==1,'Calves',region, , 1:nCauses], 3, sum)+1) 
      for (cause in causes) 
        pDet[cause, region] <- rbeta(1, sum(data.array[ ,'Calves', region, , cause])+1,  
                                     sum(U[,region, , cause])+1)
    }
  }
  cold_mort <- matrix(0, nHabitats, nSeverities, dimnames = list(habitats, severities))
  cold_mort[1:2, 'Normal'] <- runif(2, 0, 0.1)
  cold_mort[1:2, 'Cold'] <- runif(2, 0, 0.2)
  cold_mort[, 'Severe'] <- runif(3, 0, 1 - max(baseMort))
  tide_mort <- c(0, runif(2, 0, c(0.1, 0.3)))
  list(pi0 = pi * nCauses, pDet = pDet, X = U + data.array[,'Calves', , , 1:nCauses], 
       tide_mort = tide_mort, cold_mort = cold_mort) 
}
##
set.seed(0)
inits0 <- inits.calf(data.array)
inits0$pDet

## Nimble model
fraction.model.calf <- nimbleModel(fraction.code.calf, constants = constants.fraction.calf, data = data.fraction.calf, inits = inits0)
fraction.comp.calf <- compileNimble(fraction.model.calf)  ## DT taking out to skip compile

fraction.model.calf$calculate()
## [1] -1456.4
fraction.comp.calf$calculate()
## [1] -1456.4

## DT just configure for pDet[1, 1]
nodes <- c('pDet[1, 1]')  ## DT just print for nodes
fraction.mcmcConf.calf <- configureMCMC(fraction.model.calf, nodes = nodes)
fraction.mcmcConf.calf$printSamplers()  ## DT just print for nodes

## Configure, set up, and compile MCMC
fraction.mcmcConf.calf <- configureMCMC(fraction.model.calf)
fraction.mcmcConf.calf$printSamplers()

##
## NEW
## this is necessary, too, to assign the RW_multinomial samplers:
# fraction.mcmcConf.calf$removeSamplers('X')
# for(node in fraction.model.calf$expandNodeNames('X'))
#   fraction.mcmcConf.calf$addSampler(target = node, type = 'RW_multinomial')
##
##fraction.mcmcConf5$printSamplers()
##
fraction.mcmcConf.calf$getMonitors()
fraction.mcmcConf.calf$resetMonitors()
fraction.mcmcConf.calf$addMonitors('pi')
fraction.mcmcConf.calf$addMonitors('pDet')
fraction.mcmcConf.calf$addMonitors('tide_mort')
fraction.mcmcConf.calf$addMonitors('cold_mort')
##
fractionMCMC.calf <- buildMCMC(fraction.mcmcConf.calf)
CfractionMCMC.calf <- compileNimble(fractionMCMC.calf, project = fraction.model.calf, resetFunctions = TRUE)


## nodes pegged at (or near) zero
nodesToExclude <- c('cold_mort[3, 1]', 'cold_mort[3, 2]', 'pDet[5, 2]', 'pi0[5, 2]', 'tide_mort[1]', 'pi[5, 2]')

runNIMBLE <- function(seed) {
  message('beginning chain ', i, '...')
  set.seed(seed)
  inits1 <- inits.calf(data.array = data.array)
  fraction.comp.calf$setInits(inits1)
  fraction.comp.calf$calculate()
  CfractionMCMC.calf$run(niter)
  samples <- as.matrix(CfractionMCMC.calf$mvSamples)
  nodeIndToExclude <- which(colnames(samples) %in% nodesToExclude)
  samples <- samples[, -nodeIndToExclude]  ## remove nodes pegged at (or near) zero
  samples <- samples[(nburn+1):nrow(samples), ]  ## remove burnin
  return(samples)
}

samplesList <- vector('list', nChains)

options(scipen = 999)   ## used for printing output file names
saveFileName <- paste0('real.calf.pDet.niter', niter, '.RData')

for(i in 1:nChains) {  
  samplesList[[i]] <- runNIMBLE(i)
  save(list = c('samplesList'), file = saveFileName)
}







