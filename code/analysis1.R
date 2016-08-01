
##niter <- '500'
##niter <- '100000'
niter <- '250000'
inputFileName <- paste0('niter', niter, '.RData')
load(inputFileName)
options(scipen = 999)    ## used for printing output file names
source('defs.R')      ## need this for plots below

lapply(samplesList, dim)

library(coda)
mcmcs <- as.mcmc.list(lapply(samplesList, as.mcmc))

gelman.diag(mcmcs, transform = TRUE)
## Multivariate psrf
## 1.01

samples <- samplesList[[1]]

sort(apply(samples, 2, effectiveSize))
## minimum ESS from 200,000 samples is about 900 (reasonable)


for(j in 1:3) {
    samples <- samplesList[[j]]
    ##
    ## make traceplots and posterior density plots
    samplesPlot(samples, c('cold_mort[1, 1]', 'cold_mort[2, 1]', 'cold_mort[1, 2]', 'cold_mort[2, 2]', 'cold_mort[1, 3]', 'cold_mort[2, 3]', 'cold_mort[3, 3]'))
    dev.copy2pdf(file = paste0('cold_mort_', niter, '.pdf'))
    samplesPlot(samples, c('tide_mort[2]', 'tide_mort[3]'))
    dev.copy2pdf(file = paste0('tide_mort_', niter, '.pdf'))
    for(node in c('p', 'pi')) {
        for(area in 1:4) {
            nodenames <- paste0(node, '[', 1:6, ', ', area, ']')
            nodenames <- setdiff(nodenames, c('p[6, 2]', 'pi[6, 2]'))
            samplesPlot(samples, nodenames)
            regions <- c('ATL', 'USJ', 'NW', 'SW')
            ##dev.copy2pdf(file = paste0(node, '_', regions[area], '_', niter, '.pdf'))
        }
    }
    ##
}

## not doing plots for p0 and pi0, since I think p and pi
## are the actual quantities of interest.

summary(samples)


## posterior mean and median and 95% credible intervals
a <- cbind(
    apply(samples, 2, function(x) quantile(x, 0.025)),
    apply(samples, 2, mean),
    apply(samples, 2, median),
    apply(samples, 2, function(x) quantile(x, 0.975))
)
colnames(a) <- c('2.5%', 'mean', 'median', '97.5%')
print(a)



