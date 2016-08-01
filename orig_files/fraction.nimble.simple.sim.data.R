
## change 'fast' to TRUE for a quick run,
## to see if updates to the model run correctly
fast <- FALSE

if(fast) {
    niter <- 5000
    nburn <- 0
} else {
    niter <- 250000
    nburn <- 50000
}    

nChains <- 1


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
mdata <- read.csv('mortalities.simulated.July.simple.csv')
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
ages <- classes[-1] # For subadult/adult analysis
nAges <- length(ages)
severities <- c('Normal', 'Cold', 'Severe')
nSeverities <- length(severities)
##
## Baseline mortality 
##base.mort <- array(c(0.123000, 0.093000, 0.100000, 0.103000, 
##                     0.027417, 0.020739, 0.022355, 0.023073,
##                     0.027417, 0.020739, 0.022355, 0.023073),
##                   dim = c(nRegions, nClasses),
##                   dimnames = list(regions, classes))
## NEW in fraction.nimble.calf6.R
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
##modTideYears
##intTideYears
##tideYears
##
## Cold and severe years
##coldDesignations <- matrix(c(
##    'Cold', rep('Normal', 4), 'Cold', 'Normal', 'Cold', rep('Normal', 6), 'Severe', rep('Normal', 3),
##    'Cold', rep('Normal', 4), 'Cold', rep('Normal', 7), 'Cold', 'Severe', 'Cold', rep('Normal', 2),
##    'Cold', 'Normal', 'Cold', 'Normal', 'Normal', 'Cold', 'Normal', 'Severe', rep('Normal', 6), 'Severe', 'Cold', rep('Normal', 2),
##    'Cold', rep('Normal', 4), 'Cold', rep('Normal', 8), 'Severe', rep('Normal', 3)),
##                           nYears, nRegions,
##                           dimnames = list(STARTYEAR:ENDYEAR, regions))
####coldDesignations
##coldYears <- matrix(as.integer(factor(coldDesignations, levels = severities)), 
##                    nYears, nRegions,
##                    dimnames = list(STARTYEAR:ENDYEAR, regions))
## NEW in fraction.nimble.calf6.R
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
##
## Uninformative priors for proportions
prior1 <- array(1, c(nClasses, nRegions, nCauses))
dimnames(prior1) <- list(classes, regions, causes)
##
## I haven't figured out how to fix this at zero in NIMBLE, but giving it a prior close to zero might work okay
prior1[,'USJ','redtide'] <- 0.001
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
#     ## Red tide effect factors and additional mortality
#     tide_mort[1] <- 0
#     tide_mort[2] ~ dunif(0, 1-baseMort[SW])
#     tide_mort[3] ~ dunif(0, 1-baseMort[SW])
#     ## I originally defined the UME factors defined separately (see commented
#     ## code).  I removed this to simplify the structure of the model slightly, but
#     ## it didn't seem to help anything.
#     ##   tideFactor[1] <- 0
#     ##   tideFactor[2] <- tide_mort[2] / baseMort[SW]
#     ##   tideFactor[3] <- tide_mort[3] / baseMort[SW]
#     ## 
#     ## Cold effect additional mortality
#     cold_mort[1, 1] ~ dunif(0, 1-baseMort[1]) # Low normal
#     cold_mort[2, 1] ~ dunif(0, 1-baseMort[1]) # Medium normal
#     cold_mort[3, 1] <- 0                      # High normal
#     cold_mort[1, 2] ~ dunif(0, 1-baseMort[1]) # Low cold
#     cold_mort[2, 2] ~ dunif(0, 1-baseMort[1]) # Medium cold
#     cold_mort[3, 2] <- 0                      # High cold
#     cold_mort[1, 3] ~ dunif(0, 1-baseMort[1]) # Low severe
#     cold_mort[2, 3] ~ dunif(0, 1-baseMort[1]) # Medium severe
#     cold_mort[3, 3] ~ dunif(0, 1-baseMort[1]) # High severe
#     ##
    ## Loop over regions
    for (area in 1:nRegions) {
        ##
        ## Proportions of mortality for region area
        for (cause in 1:nCauses) {
            pi0[cause, area] ~ dgamma(prior1[area,cause], 1.0)
            p0[cause, area] ~ dgamma(prior1[area,cause], 1.0)
        }
        ## Determined
        pi[1:nCauses,area] <- pi0[1:nCauses,area] / sum(pi0[1:nCauses,area])
        ## Undetermined
        p[1:nCauses, area] <- p0[1:nCauses,area] / sum(p0[1:nCauses,area])
        ##
        for (habitat in 1:nHabitats) {
            ## Loop over years
            for (year in 1:nYears) {
#                 theta[year,area,habitat,1:nCauses] <- (pi[1:nCauses, area] + (area==SW) * (tideYears[year] > 1) * tideVector[1:nCauses] * tide_mort[tideYears[year]] / baseMort[SW] +
#                                                            (1 - (habitat == 3) * (coldYears[year, area] < 3)) * coldVector[1:nCauses] * cold_mort[habitat, coldYears[year, area]] / baseMort[area]) /
#                                                                (1 + (area==SW) * (tideYears[year] > 1) * tide_mort[tideYears[year]] / baseMort[SW] + 
#                                                                     (1 - (habitat == 3) * (coldYears[year, area] < 3)) * cold_mort[habitat, coldYears[year, area]] / baseMort[area])
#                 eta[year,area,habitat,1:nCauses] <- (p[1:nCauses, area] + (1 - (habitat == 3) * (coldYears[year, area] < 3)) * coldVector[1:nCauses] * cold_mort[habitat, coldYears[year, area]] / baseMort[area]) /
#                     (1 + (1 - (habitat == 3) * (coldYears[year, area] < 3)) * cold_mort[habitat, coldYears[year, area]] / baseMort[area])
                ## Determined carcasses due to each cause
                ## Ultimately we want this statement for X (total carcasses) instead of data1, where 
                ## X[year,area,habitat,1:nCauses] <- data1[year,area,habitat,1:nCauses] + U[year,area,habitat,1:nCauses] AND
                ## X[year,area,habitat,1:nCauses] ~ dmulti(prob = theta[year,area,habitat,1:nCauses], size = totals[year,area,habitat])
                ## with totals redefined to include undetermined carcasses
                ##data1[year,area,habitat,1:nCauses] ~ dmulti(prob = theta[year,area,habitat,1:nCauses],
                ##                                            size = totals[year,area,habitat])
                ## Undetermined carcasses due to each cause
                U[year,area,habitat,1:nCauses] ~ dmulti(prob = p[1:nCauses, area],
                                                        size = undet[year,area,habitat])
                ## NEW
                ## use of custom distribution:
                zeros[year,area,habitat] ~ dmultiSum(prob = pi[1:nCauses, area],
                                                     x1 = data1[year,area,habitat,1:nCauses],
                                                     x2 = U[year,area,habitat,1:nCauses])
            }
        }
    }
})

constants.fraction.calf6 <- list(
    data1 = data.array[ , 'Calves', , , 1:nCauses],
    nCauses = nCauses,
    nRegions = nRegions, 
    nYears = nYears,
#    nSeverities = nSeverities,
    nHabitats = nHabitats,
#    SW = 4,
#    tideVector = c(0, 0, 0, 0, 1, 0),
#    coldVector = c(0, 0, 0, 1, 0, 0),
    prior1 = prior1['Calves', , ], 
#    tideYears = tideYears,
#    coldYears = coldYears,
    ##totals = apply(data.array[ , 'Calves', , , 1:nCauses], 1:3, sum),  ## NEW no longer necessary
#    baseMort = as.vector(base.mort[,'Calves']),
    undet = data.array[ , 'Calves', , , 1 + nCauses]
)
##
## NEW
data6 <- list(
    zeros = array(0, c(nYears,nRegions,nHabitats))
)
##
inits.calf6 <- function() {
    p <- pi <- matrix(0, nCauses, nRegions)
    dimnames(p) <- dimnames(pi) <- list(causes, regions)
    U <- array(NA, c(nYears, nRegions, nHabitats, nCauses))
    dimnames(U) <- list(STARTYEAR:ENDYEAR, regions, habitats, causes)
    baseMort <- base.mort[,'Calves']
    for (region in regions) {
        if (region=='USJ') {
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
            pi[-REDTIDE,region] <- rdiric(1, apply(data.array[,'Calves',
                                                                   region, , -c(REDTIDE, nCauses + 1)], 3, sum)+1)
            pi[REDTIDE, region] <- 1e-100
            pi[REDTIDE - 1, region] <- pi[REDTIDE - 1, region] - 1e-100
        } else {
            p[,region] <- rdiric(1, rep(1, nCauses))
            for (habitat in habitats) {
                for (yr in STARTYEAR:ENDYEAR) {
                    year <- as.character(yr)
                    U[year, region, habitat, ] <- rmultinom(1, data.array[year, 'Calves',
                                                                          region, habitat, nCauses+1], p[ , region])
                }
            }
            pi[,region] <- rdiric(1, apply(data.array[,'Calves',region, , 1:nCauses], 3, sum)+1) 
        }
    }
#     cold_mort <- matrix(0, nHabitats, nSeverities, dimnames = list(habitats, severities))
#     cold_mort[1:2, 'Normal'] <- runif(2, 0, 0.1)
#     cold_mort[1:2, 'Cold'] <- runif(2, 0, 0.2)
#     cold_mort[, 'Severe'] <- runif(3, 0, 1 - max(baseMort))
#     tide_mort <- c(0, runif(2, 0, c(0.1, 0.3)))
    list(pi0 = pi * nCauses, p0 = p * nCauses, U = U) #, 
#         tide_mort = tide_mort, cold_mort = cold_mort) 
}
##
set.seed(0)
inits6 <- inits.calf6()

inits.list <- vector('list', nChains)
for (chain in 1:nChains) {
  set.seed(chain)
  inits.list[[chain]] <- inits.calf6()
}

inits.list[[1]]$U[6, 4, 3, ]
## NEW now uses 'data' also
## Nimble model
fraction.model.calf <- nimbleModel(fraction.code.calf, constants = constants.fraction.calf6, data = data6, inits = inits6)
fraction.comp.calf <- compileNimble(fraction.model.calf)

fraction.model.calf$calculate()
## [1] -3992.413
fraction.comp.calf$calculate()
## [1] -3992.413

## Configure, set up, and compile MCMC
fraction.mcmcConf.calf <- configureMCMC(fraction.model.calf)
fraction.mcmcConf.calf$printSamplers()
##
## NEW
## this is necessary, too, to assign the RW_multinomial samplers:
fraction.mcmcConf.calf$removeSamplers('U')
for(node in fraction.model.calf$expandNodeNames('U'))
    fraction.mcmcConf.calf$addSampler(target = node, type = 'RW_multinomial')
##
##fraction.mcmcConf5$printSamplers()
##
fraction.mcmcConf.calf$getMonitors()
fraction.mcmcConf.calf$resetMonitors()
#fraction.mcmcConf.calf$addMonitors('tide_mort')
#fraction.mcmcConf.calf$addMonitors('cold_mort')
fraction.mcmcConf.calf$addMonitors('pi')
fraction.mcmcConf.calf$addMonitors('p')
fraction.mcmcConf.calf$addMonitors('U')
#fraction.mcmcConf.calf$addMonitors('eta')
##
fractionMCMC.calf <- buildMCMC(fraction.mcmcConf.calf)
CfractionMCMC.calf <- compileNimble(fractionMCMC.calf, project = fraction.model.calf, resetFunctions = TRUE)


## nodes pegged at (or near) zero
nodesToExclude <- c('pi[5, 2]', 'p[5, 2]',
#                    paste('eta[', rep(1:nYears, nHabitats), ', 2, ', rep(1:nHabitats, each = nYears), ', 5]', sep = ''),
                    paste('U[', rep(1:nYears, nHabitats), ', 2, ', rep(1:nHabitats, each = nYears), ', 5]', sep = ''))
nodesToExclude

runNIMBLE <- function(seed) {
    message('beginning chain ', i, '...')
    set.seed(seed)
    inits6 <- inits.calf6()
    fraction.comp.calf$setInits(inits6)
    fraction.comp.calf$calculate()
    CfractionMCMC.calf$run(niter)
    samples <- as.matrix(CfractionMCMC.calf$mvSamples)
    nodeIndToExclude <- which(colnames(samples) %in% nodesToExclude)
    samples <- samples[, -nodeIndToExclude]  ## remove nodes pegged at (or near) zero
    samples <- samples[(nburn+1):nrow(samples), ]  ## remove burnin
    return(samples)
}

samplesList <- vector('list', nChains)

for(i in 1:nChains)   samplesList[[i]] <- runNIMBLE(i)

options(scipen = 999)   ## used for printing output file names
saveFileName <- paste0('sim3.niter', niter, '.RData')

save(list = c('samplesList'), file = saveFileName)





