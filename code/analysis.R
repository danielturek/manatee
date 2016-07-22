
niter <- 5000
inputFileName <- paste0('niter', niter, '.RData')
load(inputFileName)

library(coda)
mcmcs <- as.mcmc.list(lapply(samplesList, as.mcmc))

gelman.diag(mcmcs, transform = TRUE)


## utility for plotting MCMC samples
samplesPlot <- function(samples, ind=1:ncol(samples), burnin=NULL, width=7, height=4, legend=TRUE, legend.location='topright') {
    dev.new(height=height, width=width)
    par(mfrow=c(1,2), cex=0.7, cex.main=1.5, lab=c(3,3,7), mgp=c(0,0.6,0), mar=c(2,1,2,1), oma=c(0,0,0,0), tcl=-0.3, yaxt='n', bty='l')
    samples <- samples[, ind, drop=FALSE]
    if(!is.null(burnin))
        samples <- samples[(burnin+1):dim(samples)[1], , drop=FALSE]
    nparam <- ncol(samples)
    rng <- range(samples)
    plot(1:nrow(samples), ylim=rng, type='n', main='Traceplots', xlab='', ylab='')
    for(i in 1:nparam)
        lines(samples[,i], col=rainbow(nparam, alpha=0.75)[i])
    xMin <- xMax <- yMax <- NULL
    for(i in 1:nparam) {
        d <- density(samples[,i])
        xMin <- min(xMin,d$x); xMax <- max(xMax,d$x); yMax <- max(yMax, d$y) }
    plot(1, xlim=c(xMin,xMax), ylim=c(0,yMax), type='n', main='Posterior Densities', xlab='', ylab='')
    alpha_density <- 0.2
    for(i in 1:nparam)
        polygon(density(samples[,i]), col=rainbow(nparam, alpha=alpha_density)[i], border=rainbow(nparam, alpha=alpha_density)[i])
    if(legend & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
        legend(legend=dimnames(samples)[[2]], fill=rainbow(nparam, alpha=0.5), bty='n', x=legend.location)
}


samplesPlot(samples, c('cold_mort[1, 1]', 'cold_mort[2, 1]', 'cold_mort[1, 2]', 'cold_mort[2, 2]', 'cold_mort[1, 3]', 'cold_mort[2, 3]', 'cold_mort[3, 3]'))
dev.copy2pdf(file = paste0('cold_mort', niter, '.pdf'))

samplesPlot(samples, c('tide_mort[2]', 'tide_mort[3]'))
dev.copy2pdf(file = paste0('tide_mort', niter, '.pdf'))

for(node in c('p', 'pi')) {
    for(area in 1:4) {
        nodenames <- paste0(node, '[', 1:6, ', ', area, ']')
        nodenames <- setdiff(nodenames, c('p[6, 2]', 'pi[6, 2]'))
        samplesPlot(samples, nodenames)
        regions <- c('ATL', 'USJ', 'NW', 'SW')
        dev.copy2pdf(file = paste0(node, '_', regions[area], niter, '.pdf'))
    }
}

## not doing plots for p0 and pi0, since I think p and pi
## are genuinely the quantities of interest, here.

samples <- samplesList[[1]]
sort(apply(samples, 2, effectiveSize))



