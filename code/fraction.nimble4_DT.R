## Packages
library(plyr)
library(VGAM)
library(nimble)
library(igraph)

## Import data
mdata <- read.csv('mortalities.simulated.April.pDet.const.csv')
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
causes2 <- c('watercraft', 'wcs', 'debris', 'cold', 'other', 'redtide')
full.causes2 <- c(causes2, 'undetermined')
## Age classes
classes = c('Calves', 'Subadults', 'Adults')
classes2 <- c('Calves1', 'Calves2', 'Subadults', 'Adults')
nClasses <- length(classes)
severities <- c('Normal', 'Cold', 'Severe')
nSeverities <- length(severities)

## Baseline mortality 
base.mort <- array(c(0.123000, 0.093000, 0.100000, 0.103000, 
                     0.027417, 0.020739, 0.022355, 0.023073,
                     0.027417, 0.020739, 0.022355, 0.023073),
                   dim = c(nRegions, nClasses),
                   dimnames = list(regions, classes))

data <- transform(mdata,
                  area    = factor(area,    levels = regions),
                  age     = factor(age,     levels = classes),
                  habitat = factor(habitat, levels = habitats))


## Constants
CALF <- 1
SUBADULT <- 2
ADULT <- 3
STARTYEAR <- 1996
ENDYEAR <- 2013
nYears <- ENDYEAR - STARTYEAR + 1

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
modTideYears
intTideYears
tideYears

## Cold and severe years
coldDesignations <- matrix(c('Cold', rep('Normal', 4), 'Cold', 'Normal', 'Cold', rep('Normal', 6), 'Severe', rep('Normal', 3),
                             'Cold', rep('Normal', 4), 'Cold', rep('Normal', 7), 'Cold', 'Severe', 'Cold', rep('Normal', 2),
                             'Cold', 'Normal', 'Cold', 'Normal', 'Normal', 'Cold', 'Normal', 'Severe', rep('Normal', 6), 'Severe', 'Cold', rep('Normal', 2),
                             'Cold', rep('Normal', 4), 'Cold', rep('Normal', 8), 'Severe', rep('Normal', 3)),
                           nYears, nRegions, dimnames = list(STARTYEAR:ENDYEAR, regions))
coldDesignations
coldYears <- matrix(as.integer(factor(coldDesignations, levels = severities)), 
                    nYears, nRegions,
                    dimnames = list(STARTYEAR:ENDYEAR, regions))

## Uninformative priors for proportions
prior1 <- array(1, c(nClasses, nRegions, nCauses))
dimnames(prior1) <- list(classes, regions, causes2)
## I haven't figured out how to fix this at zero in NIMBLE, but giving it a prior close to zero might work okay
prior1[,'USJ','redtide'] <- 0.001
rdirch(1, prior1['Calves','USJ',])

## Move from data frame to array
data.array <- array(NA, c(nYears, nClasses, nRegions, nHabitats, nCauses + 1))
dimnames(data.array) <- list(STARTYEAR:ENDYEAR, classes, regions, habitats, full.causes2)
for (class in classes) {
    for (region in regions) {
        for (qual in habitats) 
            data.array[ , class, region, qual, ] <- as.matrix(subset(data, age==class & area==region & habitat == qual)[,full.causes2])
    }
}

fraction.code4 <- nimbleCode({
    
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
                                        # 
    ## Cold effect additional mortality
    cold_mort[1, 1] ~ dunif(0, 1-baseMort[1]) # Low normal
    cold_mort[2, 1] ~ dunif(0, 1-baseMort[1]) # Medium normal
    cold_mort[3, 1] <- 0                        # High normal
    cold_mort[1, 2] ~ dunif(0, 1-baseMort[1]) # Low cold
    cold_mort[2, 2] ~ dunif(0, 1-baseMort[1]) # Medium cold
    cold_mort[3, 2] <- 0                        # High cold
    cold_mort[1, 3] ~ dunif(0, 1-baseMort[1]) # Low severe
    cold_mort[2, 3] ~ dunif(0, 1-baseMort[1]) # Medium severe
    cold_mort[3, 3] ~ dunif(0, 1-baseMort[1]) # High severe
    
    
    ## Loop over regions
    for (area in 1:nRegions) {
        
        ## Proportions of mortality for region area
        for (cause in 1:nCauses) {
            pi0[cause, area] ~ dgamma(prior1[area,cause], 1.0)
            p0[cause, area] ~ dgamma(prior1[area,cause], 1.0)
        }
        ## Determined
        pi[1:nCauses,area] <- pi0[1:nCauses,area] / sum(pi0[1:nCauses,area])
        ## Undetermined
        p[1:nCauses, area] <- p0[1:nCauses,area] / sum(p0[1:nCauses,area])
        
        for (habitat in 1:nHabitats) {
            ## Loop over years
            for (year in 1:nYears) {
                theta[year,area,habitat,1:nCauses] <- (pi[1:nCauses, area] + (area==SW) * (tideYears[year] > 1) * tideVector[1:nCauses] * tide_mort[tideYears[year]] / baseMort[SW] +
                                                           (1 - (habitat == 3) * (coldYears[year, area] < 3)) * coldVector[1:nCauses] * cold_mort[habitat, coldYears[year, area]] / baseMort[area]) /
                                                               (1 + (area==SW) * (tideYears[year] > 1) * tide_mort[tideYears[year]] / baseMort[SW] + 
                                                                    (1 - (habitat == 3) * (coldYears[year, area] < 3)) * cold_mort[habitat, coldYears[year, area]] / baseMort[area])
                ## Determined carcasses due to each cause
                ## Ultimately we want this statement for X (total carcasses) instead of data1, where 
                ## X[year,area,habitat,1:nCauses] <- data1[year,area,habitat,1:nCauses] + U[year,area,habitat,1:nCauses] AND
                ## X[year,area,habitat,1:nCauses] ~ dmulti(prob = theta[year,area,habitat,1:nCauses], size = totals[year,area,habitat])
                ## with totals redefined to include undetermined carcasses
                data1[year,area,habitat,1:nCauses] ~ dmulti(prob = theta[year,area,habitat,1:nCauses],
                                                            size = totals[year,area,habitat])
                ## Undetermined carcasses due to each cause
                U[year,area,habitat,1:nCauses] ~ dmulti(prob = p[1:nCauses, area],
                                                        size = undet[year,area,habitat])
            }
        }
    }
})


## This very simplified version of the model (and other ones where pi isn't
## modified into theta) seem to work fine.
det.code1 <- nimbleCode({
    ## Loop over regions
    for (area in 1:nRegions) {
        
        ## Proportions of mortality for region area
        pi[1:nCauses,area] ~ ddirch(alpha = prior1[age_class,area,1:nCauses])
        
        for (habitat in 1:nHabitats) {
            ## Loop over years
            for (year in 1:nYears) {
                ## Multinomial number of carcasses due to each cause
                data1[year,area,habitat,1:nCauses] ~ dmulti(prob = pi[1:nCauses,area],
                                                            size = totals[year,area,habitat])
            }
        }
    }
})


## Data structures for Nimble
data.fraction.calf4 <- list(data1=data.array[ , 'Calves', , , 1:nCauses], nCauses=nCauses, nRegions=nRegions, 
                            nYears = nYears, nSeverities = nSeverities, nHabitats = nHabitats,
                            SW = 4, tideVector = c(0, 0, 0, 0, 0, 1), coldVector = c(0, 0, 0, 1, 0, 0),
                            prior1 = prior1['Calves', , ], 
                            tideYears = tideYears, coldYears = coldYears,
                            totals = apply(data.array[ , 'Calves', , , 1:nCauses], 1:3, sum), 
                            baseMort = as.vector(base.mort[,'Calves']),
                            undet = data.array[ , 'Calves', , , 1 + nCauses])
data.det.calf <- list(data1=data.array[, 'Calves', , , 1:nCauses], nCauses=nCauses, nRegions=nRegions, 
                      nYears = nYears, nSeverities = nSeverities, nHabitats = nHabitats,
                      age_class = CALF, prior1 = prior1, 
                      totals = apply(data.array[ , 'Calves', , , 1:nCauses], 1:3, sum)) 


inits.calf4 <- function(){
    p <- pi <- matrix(0, nCauses, nRegions)
    dimnames(p) <- dimnames(pi) <- list(causes2, regions)
    U <- array(NA, c(nYears, nRegions, nHabitats, nCauses))
    dimnames(U) <- list(STARTYEAR:ENDYEAR, regions, habitats, causes2)
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
        }
        else if (region=='USJ') {
            p[1:(nCauses-1), region] <- rdiric(1, rep(1, nCauses - 1))
            p[nCauses, region] <- 1e-100
            p[nCauses - 1, region] <- p[nCauses - 1, region] - 1e-100
            for (habitat in habitats) {
                for (yr in STARTYEAR:ENDYEAR) {
                    year <- as.character(yr)
                    U[year, region, habitat, 1:(nCauses-1)] <- rmultinom(1, data.array[year,'Calves',
                                                                                       region, habitat, nCauses+1], p[1:(nCauses-1), region])
                    U[year, region, habitat, nCauses] <- 0
                }
            }
            pi[1:(nCauses-1),region] <- rdiric(1, apply(data.array[coldYears[,region]==1,'Calves',
                                                                   region, , 1:(nCauses-1)], 3, sum)+1)
            pi[nCauses, region] <- 1e-100
            pi[nCauses - 1, region] <- pi[nCauses - 1, region] - 1e-100
        } else  {
            p[,region] <- rdiric(1, rep(1, nCauses))
            for (habitat in habitats) {
                for (yr in STARTYEAR:ENDYEAR) {
                    year <- as.character(yr)
                    U[year, region, habitat, ] <- rmultinom(1, data.array[year, 'Calves',
                                                                          region, habitat, nCauses+1], p[ , region])
                }
            }
            pi[,region] <- rdiric(1, apply(data.array[coldYears[,region]==1,'Calves',region, , 1:nCauses], 3, sum)+1) 
        }
    }
    cold_mort <- matrix(0, nHabitats, nSeverities, dimnames = list(habitats, severities))
    cold_mort[1:2, 'Normal'] <- runif(2, 0, 0.1)
    cold_mort[1:2, 'Cold'] <- runif(2, 0, 0.2)
    cold_mort[, 'Severe'] <- runif(3, 0, 1 - max(baseMort))
    tide_mort <- c(0, runif(2, 0, c(0.1, 0.3)))
    list(pi0 = pi * nCauses, p0 = p * nCauses, U = U, 
         tide_mort = tide_mort, cold_mort = cold_mort) 
}

inits.det.calf <- function(){
    pi <- matrix(0, nCauses, nRegions)
    dimnames(pi) <- list(causes2, regions)
    for (region in regions) {
        pi[,region] <- rdiric(1, apply(data.array[ ,'Calves',region, , 1:nCauses], 3, sum)+1)
    }
    list(pi = pi) 
}

## Nimble model
inits4 <- inits.calf4()
fraction.model4 <- nimbleModel(fraction.code4, constants = data.fraction.calf4,
                               inits = inits4)

fraction.comp4 <- compileNimble(fraction.model4)

## Configure, set up, and compile MCMC
fraction.mcmcConf4 <- configureMCMC(fraction.model4)
fraction.mcmcConf4$printSamplers()
## I've tried some different options for sampling.  So far the most I can do is
## spread the problem from pi to other variables (see commented code below).
## fraction.mcmcConf3$removeSamplers('tide_mort')
## fraction.mcmcConf3$removeSamplers('pi')
## fraction.mcmcConf3$removeSamplers('cold_mort')
## fraction.model3$expandNodeNames('tide_mort')
## fraction.mcmcConf3$addSampler(type = 'RW_block', 
##                               target = c(fraction.model3$expandNodeNames('tide_mort'),
##                                          fraction.model3$expandNodeNames('pi'),
##                                          fraction.model3$expandNodeNames('cold_mort')))
## fraction.mcmcConf3$printSamplers()
fraction.mcmcConf4$getMonitors()
fraction.mcmcConf4$addMonitors('pi')
fraction.mcmcConf4$addMonitors('p')

fractionMCMC4 <- buildMCMC(fraction.mcmcConf4)
CfractionMCMC4 <- compileNimble(fractionMCMC4, project = fraction.model4, resetFunctions = T)

## Run the model
print(system.time(CfractionMCMC4$run(10000)))

## Examine results
sample.mat <- as.matrix(CfractionMCMC4$mvSamples)
summary(sample.mat)

## Yes, these still sum to 1
summary(rowSums(sample.mat[,c('pi[1, 1]', 'pi[2, 1]', 'pi[3, 1]', 'pi[4, 1]', 'pi[5, 1]', 'pi[6, 1]')]))
summary(rowSums(sample.mat[,c('pi[1, 2]', 'pi[2, 2]', 'pi[3, 2]', 'pi[4, 2]', 'pi[5, 2]', 'pi[6, 2]')]))

## Much simpler model
det.model1 <- nimbleModel(det.code1, constants = data.det.calf,
                          inits = inits.det.calf(),
                          dimensions = list(pi = c(nCauses, nRegions)))

det.comp1 <- compileNimble(det.model1)

mcmcConf <- configureMCMC(det.model1)
mcmcConf$printSamplers()
detMCMC <- buildMCMC(mcmcConf)
CdetMCMC <- compileNimble(detMCMC, project = det.model1)
CdetMCMC$run(50000)
summary(as.matrix(CdetMCMC$mvSamples))

