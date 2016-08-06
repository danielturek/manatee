# Read in and process individual-based recovery data

# Constants and categories
# FL regions



regions <- c('ATL', 'USJ', 'NW', 'SW')
nRegions <- length(regions)
# Winter habitat quality
habitats <- c('Low','Medium','High')
nHabitats <- length(habitats)
# Causes of death
causes <- c('watercraft', 'wcs', 'debris', 'cold', 'redtide', 'other')
nCauses <- length(causes)
full.causes <- c(causes, 'undetermined')
# Age classes
classes = c('Calves', 'Subadults', 'Adults')
#classes <- c('Calves', 'Adults')
nClasses <- length(classes)
STARTYEAR <- 1996; ENDYEAR <- 2013
years <- STARTYEAR:ENDYEAR
nYears <- length(years)
startMonth <-c(12,12,12,12); startDay <- c(1,1,1,1) # KB's version, based on survival analysis interval starts
nWCSYears <- 5
nProtectYears <- nYears - nWCSYears
wcsProtectYear <- 2001
names(startDay) <- names(startMonth) <- regions
#cold.years.NW <- c(1996, 1998, 2001, 2011)
#severe.years.NW <- c(2003, 2010)
##coldDesignations
severities <- c('Normal', 'Cold', 'Severe')
nSeverities <- length(severities)
coldDesignations <- matrix(
  c("Cold", rep("Normal", 4), "Cold", "Normal", "Cold", rep("Normal", 6), "Severe", "Cold", rep("Normal", 2),
    "Cold", rep("Normal", 4), "Cold", rep("Normal", 7), "Cold", "Severe", "Cold", rep("Normal", 2),
    "Cold", "Normal", "Cold", "Normal", "Normal", "Cold", "Normal", "Severe", rep("Normal", 6), "Severe", "Cold", rep("Normal", 2),
    "Cold", rep("Normal", 4), "Cold", rep("Normal", 8), "Severe", rep("Normal", 3)),
  nYears, nRegions,
  dimnames = list(STARTYEAR:ENDYEAR, regions))
coldDesignations
coldYears <- matrix(as.integer(factor(coldDesignations, levels = severities)), 
                    nYears, nRegions,
                    dimnames = list(STARTYEAR:ENDYEAR, regions))

intensities <- c("Baseline", "Moderate", "Intense")
nIntensities <- length(intensities)
tideDesignations <- matrix(c(rep("Baseline", nYears * 3), 
                             "Intense", rep("Baseline", 5), rep("Moderate", 2), "Baseline", rep("Moderate", 2), rep("Baseline", 5), "Moderate", "Intense"),
                           nYears, nRegions, dimnames = list(years, regions))
tideYears <- matrix(as.integer(factor(tideDesignations, levels = intensities)), 
                    nYears, nRegions,
                    dimnames = list(STARTYEAR:ENDYEAR, regions))



# Simulation parameters
# base.mort <- array(c(0.123, 0.093, 0.1, 0.103, 
#                      0.027417, 0.020739, 0.022355, 0.023073,
#                      0.027417, 0.020739, 0.022355, 0.023073), c(nRegions, nClasses),
#                    dimnames = list(regions, classes))
adult.base.mort <- 1 - c(0.972883, 0.978949, 0.977913, 0.977058)
calf1.mort.ratio <- 0.190/0.031
calf2.mort.ratio <- 0.085/0.031
base.mort <- array(c(1 - sqrt((1 - calf1.mort.ratio * adult.base.mort) * (1 - calf2.mort.ratio * adult.base.mort)), adult.base.mort, adult.base.mort),
                   dim = c(nRegions, nClasses),
                   dimnames = list(regions, classes))
base.mort
base.pi <- array(c(0.357, 0.015, 0.008, 0.21, 0.01, 0.4,
                   0.469, 0.048, 0.048, 0.386, 0, 0.049,
                   0.496, 0.032, 0.016, 0.149, 0.016, 0.291,
                   0.405, 0.008, 0.008, 0.068, 0.166, 0.345,
                   rep(c(0.357, 0.015, 0.008, 0.21, 0.01, 0.4,
                   0.469, 0.048, 0.048, 0.386, 0, 0.049,
                   0.496, 0.032, 0.016, 0.149, 0.016, 0.291,
                   0.522, 0.03, 0.011, 0.032, 0.284, 0.121), 2)), 
                 c(nCauses, nRegions, nClasses), dimnames = list(causes, regions, classes))
apply(base.pi, 2:3, sum)
base.pi
tide.mort <- t(matrix(c(rep(0, nClasses), 0.05, 0.02, 0.02, 0.18, 0.1, 0.1),
                    nClasses, nIntensities, dimnames = list(classes, intensities)))
cold.mort <- array(c(0.056, 0.004, 0, 0.094, 0.042, 0, 0.645, 0.410, 0.169,
                     0.005, 0, 0, 0.008, 0.007, 0, 0.15, 0.107, 0.063, 
                     rep(0, 3), 0.002, 0, 0, 0.053, 0.033, 0.023), 
                   c(nHabitats, nSeverities, nClasses), dimnames = list(habitats, severities, classes))
c.mort <- array(c(rep(0, nRegions), rep(c(0.007, 0, 0, 0), 2)), c(nRegions, nClasses),
                dimnames = list(regions, classes))
total.mort <- array(rep(base.mort, each = nYears * nHabitats), 
                    c(nYears, nHabitats, nRegions, nClasses), 
                    dimnames = list(years, habitats, regions, classes))
csm.mort <- array(rep(base.mort, each = nCauses * nYears * nHabitats) *
                    rep(base.pi, each = nYears * nHabitats), 
                  c(nYears, nHabitats, nCauses, nRegions, nClasses), 
                  dimnames = list(years, habitats, causes, regions, classes))
for (year in years) {
  for (age in classes) {
    for (habitat in habitats) {
      for (region in regions) {
        if (year < wcsProtectYear) {
          total.mort[as.character(year), habitat, region, age] <- 
            total.mort[as.character(year), habitat, region, age] + c.mort[region, age]
          csm.mort[as.character(year), habitat, 'wcs', region, age] <-
            csm.mort[as.character(year), habitat, 'wcs', region, age] + c.mort[region, age]
        }
        total.mort[as.character(year), habitat, region, age] <- 
          total.mort[as.character(year), habitat, region, age] + 
          cold.mort[habitat, coldYears[as.character(year), region], age]
        csm.mort[as.character(year), habitat, 'cold', region, age] <- 
          csm.mort[as.character(year), habitat, 'cold', region, age] + 
          cold.mort[habitat, coldYears[as.character(year), region], age]
        total.mort[as.character(year), habitat, region, age] <- 
          total.mort[as.character(year), habitat, region, age] + 
          tide.mort[tideYears[as.character(year), region], age]
        csm.mort[as.character(year), habitat, 'redtide', region, age] <- 
          csm.mort[as.character(year), habitat, 'redtide', region, age] + 
          tide.mort[tideYears[as.character(year), region], age]
      }
    }
  }
}
total.mort
max(total.mort)
theta <- prop.table(csm.mort, c(1:2, 4:5))
theta
pDet0 <- 2/3
base.p <- array(c(0.357, 0.015, 0.008, 0.21, 0.01, 0.4,
                   0.4, 0.1, 0.1, 0.35, 0, 0.05,
                   0.4, 0.03, 0.02, 0.2, 0.01, 0.34,
                   0.48, 0.01, 0.01, 0.1, 0.05, 0.35,
                   rep(c(0.357, 0.015, 0.008, 0.21, 0.01, 0.4,
                         0.469, 0.048, 0.048, 0.386, 0, 0.049,
                         0.496, 0.032, 0.016, 0.149, 0.016, 0.291,
                         0.55, 0.05, 0.01, 0.05, 0.05, 0.29), 2)), 
                 c(nCauses, nRegions, nClasses), dimnames = list(causes, regions, classes))
apply(base.p, 2:3, sum)
base.p - base.pi
tide.frac <- tide.mort / rep(base.mort["SW", ], each = nIntensities)
cold.frac <- aperm(array(rep(aperm(cold.mort, c(3, 2, 1)), each = nRegions) / rep(base.mort, nSeverities * nHabitats),
                   c(nRegions, nClasses, nSeverities, nHabitats), 
                   dimnames = list(regions, classes, severities, habitats)), c(3, 4, 1, 2))
cold.frac
c.frac <- c.mort / base.mort

eta <- array(rep(base.p, each = nHabitats * nYears), c(nYears, nHabitats, nCauses, nRegions, nClasses), 
             dimnames = list(years, habitats, causes, regions, classes))
theta2 <- array(rep(base.pi, each = nHabitats * nYears), c(nYears, nHabitats, nCauses, nRegions, nClasses), 
             dimnames = list(years, habitats, causes, regions, classes))
for (year in years) {
  for (age in classes) {
    for (habitat in habitats) {
      for (region in regions) {
        if (year < wcsProtectYear) {
          if (region == "SW") {
            theta2[as.character(year), habitat, , region, age] <- 
              (theta2[as.character(year), habitat, , region, age] + c(0, c.frac[region, age], rep(0, 4)) +
                 c(rep(0, 3), cold.frac[coldYears[as.character(year), region], habitat, region, age], 0, 0) +
                 c(rep(0, 4), tide.frac[tideYears[as.character(year), region], age], 0)) /
              (1 + c.frac[region, age] + cold.frac[coldYears[as.character(year), region], habitat, region, age] +
                 tide.frac[tideYears[as.character(year), region], age])
          }
          else {
            theta2[as.character(year), habitat, , region, age] <- 
              (theta2[as.character(year), habitat, , region, age] + c(0, c.frac[region, age], rep(0, 4)) +
                 c(rep(0, 3), cold.frac[coldYears[as.character(year), region], habitat, region, age], 0, 0)) /
              (1 + c.frac[region, age] + cold.frac[coldYears[as.character(year), region], habitat, region, age])
          }
          eta[as.character(year), habitat, , region, age] <- 
            (eta[as.character(year), habitat, , region, age] + c(0, c.frac[region, age], rep(0, 4)) +
               c(rep(0, 3), cold.frac[coldYears[as.character(year), region], habitat, region, age], 0, 0)) /
            (1 + c.frac[region, age] + cold.frac[coldYears[as.character(year), region], habitat, region, age])
        }
        else {
          if (region == "SW") {
            theta2[as.character(year), habitat, , region, age] <- 
              (theta2[as.character(year), habitat, , region, age] +
                 c(rep(0, 3), cold.frac[coldYears[as.character(year), region], habitat, region, age], 0, 0) +
                 c(rep(0, 4), tide.frac[tideYears[as.character(year), region], age], 0)) /
              (1 + cold.frac[coldYears[as.character(year), region], habitat, region, age] +
                 tide.frac[tideYears[as.character(year), region], age])
          }
          else {
            theta2[as.character(year), habitat, , region, age] <- 
              (theta2[as.character(year), habitat, , region, age] + 
                 c(rep(0, 3), cold.frac[coldYears[as.character(year), region], habitat, region, age], 0, 0)) /
              (1 + cold.frac[coldYears[as.character(year), region], habitat, region, age])
          }
          eta[as.character(year), habitat, , region, age] <- 
            (eta[as.character(year), habitat, , region, age] + 
               c(rep(0, 3), cold.frac[coldYears[as.character(year), region], habitat, region, age], 0, 0)) /
            (1 + cold.frac[coldYears[as.character(year), region], habitat, region, age])
        }
      }
    }
  }
}
theta - theta2
min(theta - theta2)
max(theta - theta2)
apply(eta, c(1:2, 4:5), sum)

pDet0 <- 2/3
zeta <- (theta - (1 - pDet0) * eta) / pDet0
zeta[tideYears[,"SW"] > 1, , , 'SW', ] <- theta[tideYears[,"SW"] > 1, , , 'SW', ] 
zeta[tideYears[,"SW"] > 1,, 'redtide', 'SW', ] <- zeta[tideYears[,"SW"] > 1,, 'redtide', 'SW', ] + 0.1
zeta[tideYears[,"SW"] > 1,, , 'SW', ] <- zeta[tideYears[,"SW"] > 1,, , 'SW', ] / 1.1
zeta[,,,'SW', 'Calves']
apply(zeta, c(1:2, 4:5), sum)

# eta["2013", , , "SW", ] <- eta["2011", , , "SW", ]
# eta["1996", , , "SW", ] <- eta["2001", , , "SW", ]
# eta["2002", , , "SW", ] <- eta["2004", , , "SW", ]
# eta["2003", , , "SW", ] <- eta["2004", , , "SW", ]
# eta["2005", , , "SW", ] <- eta["2004", , , "SW", ]
# eta["2006", , , "SW", ] <- eta["2004", , , "SW", ]
# eta["2012", , , "SW", ] <- eta["2004", , , "SW", ]
# eta["2007", , , "ATL", ] <- eta["2006", , , "ATL", ]
# eta["2008", , , "ATL", ] <- eta["2006", , , "ATL", ]
recovery <- 0.5
N <- array(rpois(nYears * nHabitats * nRegions * nClasses, 
                 250), 
           c(nYears, nHabitats, nRegions, nClasses),
           dimnames = list(years, habitats, regions, classes))
N
apply(N, c(1, 3), sum)
#determine <- array(0.65, c(nYears, nHabitats, nCauses, nRegions, nClasses))
#undetermine <- 0.35 * eta / theta

# best.u <- function(u, theta, eta, mean.u) {
#   u.array <- array(u, c(nYears, nHabitats, nCauses, nRegions, nClasses), 
#                    dimnames = list(years, habitats, causes, regions, classes))
#   eta.new <- prop.table(theta * u.array, c(1:2, 4:5))
#   return(sum((eta - eta.new)^2) + 100 * (mean(u) - mean.u)^2)
# }
# 
# opt <- optim(rep(0.35, nYears * nHabitats * nCauses * nRegions * nClasses), best.u,
#              theta = theta, eta = eta, mean.u = 0.35, method = "L-BFGS-B",
#              lower = rep(0, nYears * nHabitats * nCauses * nRegions * nClasses),
#              upper = rep(1, nYears * nHabitats * nCauses * nRegions * nClasses))
# undetermine <- array(opt$par, c(nYears, nHabitats, nCauses, nRegions, nClasses), 
#                      dimnames = list(years, habitats, causes, regions, classes))
# max(undetermine)
# min(undetermine)
# undetermine[c('2007', '2008'), 'Low', 'redtide', 'ATL', 'Calves'] <- 0.01
# undetermine[,,'redtide',,'Adults']
# eta0[,,'redtide',,'Adults']
# which(eta != theta, arr.ind = T)

pDet <- ifelse(zeta == eta, pDet0, (theta - eta) / (zeta - eta))
dimnames(pDet)
min(pDet)
max(pDet)
pDet[pDet==0] <- pDet0
pDet[pDet==1] <- pDet0


# Run sim
set.seed(44)
data.sim <- undeterm.sim <- data.frame(area = rep(regions, nClasses * nHabitats * nYears),
                                       year = rep(years, nClasses * nHabitats, each = nRegions), 
                                        age = rep(classes, nHabitats, each = nYears * nRegions), 
                                        habitat = rep(habitats, each = nYears * nClasses * nRegions), 
                                       watercraft = 0, wcs = 0, debris = 0, cold = 0, redtide = 0, 
                                       other = 0, undetermined = 0, total = 0)
for (age in classes) {
  for (year in years) {
    for (habitat in habitats) {
      for (region in regions) {
        dead <- rbinom(1, N[as.character(year), habitat, region, age], 
                       total.mort[as.character(year), habitat, region, age] * recovery)
        tot.by.cause <- rmultinom(1, dead, theta[as.character(year), habitat, , region, age])
        undetermined.by.cause <- rbinom(nCauses, tot.by.cause, 1 - pDet[as.character(year), habitat, , region, age]) 
        data.sim[data.sim$area==region & data.sim$habitat==habitat & data.sim$age==age & data.sim$year==year, 
                 causes] <- tot.by.cause - undetermined.by.cause
        data.sim[data.sim$area==region & data.sim$habitat==habitat & data.sim$age==age & data.sim$year==year, 
                   'undetermined'] <- sum(undetermined.by.cause)
        undeterm.sim[undeterm.sim$area==region & undeterm.sim$habitat==habitat & undeterm.sim$age==age & undeterm.sim$year==year, 
                   causes] <- undetermined.by.cause
      }
    }
  }
}
data.sim$total <- rowSums(data.sim[,full.causes])
data.sim
undeterm.sim$total <- rowSums(undeterm.sim[,causes])
undeterm.sim

write.csv(data.sim, 'mortalities.simulated.July.csv', row.names=F)
write.csv(undeterm.sim, 'undeterm.simulated.July.csv', row.names=F)

# Run sim with zeta and eta
set.seed(44)
data.sim <- undeterm.sim <- data.frame(area = rep(regions, nClasses * nHabitats * nYears),
                                       year = rep(years, nClasses * nHabitats, each = nRegions), 
                                       age = rep(classes, nHabitats, each = nYears * nRegions), 
                                       habitat = rep(habitats, each = nYears * nClasses * nRegions), 
                                       watercraft = 0, wcs = 0, debris = 0, cold = 0, redtide = 0, 
                                       other = 0, undetermined = 0, total = 0)
for (age in classes) {
  for (year in years) {
    for (habitat in habitats) {
      for (region in regions) {
        dead1 <- rbinom(1, N[as.character(year), habitat, region, age], 
                       total.mort[as.character(year), habitat, region, age] * recovery * pDet0)
        if (tideYears[as.character(year), region]>1) {
          if (year == 1996)
            year0 <- "2001"
          else
            year0 <- "2000"
          dead2 <- rbinom(1, N[as.character(year), habitat, region, age], 
                          total.mort[year0, habitat, region, age] * recovery * (1 - pDet0))
        }
        else
          dead2 <- rbinom(1, N[as.character(year), habitat, region, age], 
                          total.mort[as.character(year), habitat, region, age] * recovery * (1 - pDet0))
        determined.by.cause <- rmultinom(1, dead1, zeta[as.character(year), habitat, , region, age])
        undetermined.by.cause <- rmultinom(1, dead2, eta[as.character(year), habitat, , region, age]) 
        data.sim[data.sim$area==region & data.sim$habitat==habitat & data.sim$age==age & data.sim$year==year, 
                 causes] <- determined.by.cause
        data.sim[data.sim$area==region & data.sim$habitat==habitat & data.sim$age==age & data.sim$year==year, 
                 'undetermined'] <- sum(undetermined.by.cause)
        undeterm.sim[undeterm.sim$area==region & undeterm.sim$habitat==habitat & undeterm.sim$age==age & undeterm.sim$year==year, 
                     causes] <- undetermined.by.cause
      }
    }
  }
}
data.sim$total <- rowSums(data.sim[,full.causes])
data.sim
undeterm.sim$total <- rowSums(undeterm.sim[,causes])
undeterm.sim$undetermined <- rowSums(undeterm.sim[,causes])
undeterm.sim

write.csv(data.sim, 'mortalities.simulated.July.zeta.eta.csv', row.names=F)
write.csv(undeterm.sim, 'undeterm.simulated.July.zeta.eta.csv', row.names=F)

# Run sim with p and pi
set.seed(44)
data.sim <- undeterm.sim <- data.frame(area = rep(regions, nClasses * nHabitats * nYears),
                                       year = rep(years, nClasses * nHabitats, each = nRegions), 
                                       age = rep(classes, nHabitats, each = nYears * nRegions), 
                                       habitat = rep(habitats, each = nYears * nClasses * nRegions), 
                                       watercraft = 0, wcs = 0, debris = 0, cold = 0, redtide = 0, 
                                       other = 0, undetermined = 0, total = 0)
for (age in classes) {
  for (year in years) {
    for (habitat in habitats) {
      for (region in regions) {
        dead1 <- rbinom(1, N[as.character(year), habitat, region, age], 
                        total.mort[as.character(year), habitat, region, age] * recovery * pDet0)
        dead2 <- rbinom(1, N[as.character(year), habitat, region, age], 
                        total.mort[as.character(year), habitat, region, age] * recovery * (1 - pDet0))
        determined.by.cause <- rmultinom(1, dead1, base.pi[, region, age])
        undetermined.by.cause <- rmultinom(1, dead2, base.p[ , region, age]) 
        data.sim[data.sim$area==region & data.sim$habitat==habitat & data.sim$age==age & data.sim$year==year, 
                 causes] <- determined.by.cause
        data.sim[data.sim$area==region & data.sim$habitat==habitat & data.sim$age==age & data.sim$year==year, 
                 'undetermined'] <- sum(undetermined.by.cause)
        undeterm.sim[undeterm.sim$area==region & undeterm.sim$habitat==habitat & undeterm.sim$age==age & undeterm.sim$year==year, 
                     causes] <- undetermined.by.cause
      }
    }
  }
}
data.sim$total <- rowSums(data.sim[,full.causes])
data.sim
undeterm.sim$total <- rowSums(undeterm.sim[,causes])
undeterm.sim$undetermined <- rowSums(undeterm.sim[,causes])
undeterm.sim

write.csv(data.sim, 'mortalities.simulated.July.simple.csv', row.names=F)
write.csv(undeterm.sim, 'undeterm.simulated.July.simple.csv', row.names=F)

