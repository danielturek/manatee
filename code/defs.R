

## NEW
## define a custom distribution for the determ/stoch multinomial
dmultiSum <- nimbleFunction(
    run = function(x = double(), prob = double(1), x1 = double(1), x2 = double(1), log.p = double()) {
        len <- dim(x1)[1]
        data <- x1 + x2
        size <- sum(data)
        logL <- dmulti(data[1:len], prob = prob[1:len], size = size, log = 1)
        returnType(double())
        return(logL)
    }
)


## NEW
## define a custom distribution for the determ/stoch multinomial
rmultiSum <- nimbleFunction(
    run = function(n = integer(), prob = double(1), x1 = double(1), x2 = double(1)) {
        print('not implemented')
        returnType(double())
        return(1)
    }
)


## NEW
## register the new distribution with NIMBLE
registerDistributions(list(
    dmultiSum = list(
        BUGSdist = 'dmultiSum(prob, x1, x2)',
        types = c('prob = double(1)', 'x1 = double(1)', 'x2 = double(1)'),
        discrete = TRUE
    )
))


## NEW
## we'll need this new sampler to sample the multnomial distribution for U[...].
## previously it worked out of the box, since U had no dependencies so
## a posterior_predictive sampler was used.
sampler_RW_multinomial <- nimbleFunction( 
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        targetAllNodes <- unique(model$expandNodeNames(target))
        calcNodes      <- model$getDependencies(target) 
        lTarget        <- length(targetAsScalar)
        Ntotal         <- sum(values(model,target))
        NOverL         <- Ntotal / lTarget
        Zeros             <- matrix(0, lTarget, lTarget)
        Ones              <- matrix(1, lTarget, lTarget)
        timesRan          <- Zeros
        AcceptRates       <- Zeros
        ScaleShifts       <- Zeros
        totalAdapted      <- Zeros
        timesAccepted     <- Zeros
        ENSwapMatrix      <- Ones
        ENSwapDeltaMatrix <- Ones
        RescaleThreshold  <- 0.2 * Ones
        lpProp  <- 0
        lpRev   <- 0
        Pi      <- pi 
        PiOver2 <- Pi / 2
        u       <- runif(1, 0, Pi)
        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        my_decideAndJump       <- decideAndJump(model, mvSaved, calcNodes)
        if(model$getNodeDistribution(target) != 'dmulti')   stop('can only use RW_multinomial sampler for multinomial distributions')
        if(length(targetAllNodes) > 1)                      stop('cannot use RW_multinomial sampler on more than one target')
        if(adaptive & adaptInterval < 100)                  stop('adaptInterval < 100 is not recommended for RW_multinomial sampler')
    },
    run = function() {
        for(iFROM in 1:lTarget) {            
            for(iTO in 1:(lTarget-1)) {
                if(u > PiOver2) {                
                    iFrom <- iFROM
                    iTo   <- iTO
                    if (iFrom == iTo)
                        iTo <- lTarget
                    u <<- 2 * (u - PiOver2)
                } else {
                    iFrom <- iTO
                    iTo   <- iFROM
                    if (iFrom == iTo)
                        iFrom <- lTarget
                    u <<- 2 * (PiOver2 - u)
                }
                propValueVector <- generateProposalVector(iFrom, iTo)
                lpMHR <- my_setAndCalculateDiff$run(propValueVector) + lpRev - lpProp 
                jump  <- my_decideAndJump$run(lpMHR, 0, 0, 0)
                if(adaptive)   adaptiveProcedure(jump=jump, iFrom=iFrom, iTo=iTo)
            }
        }
    },
    methods = list(
        generateProposalVector = function(iFrom = integer(), iTo = integer()) { 
            propVector <- values(model,target) 
            pSwap      <- min(1, max(1, ENSwapMatrix[iFrom,iTo]) / propVector[iFrom]) 
            nSwap      <- rbinom(n=1,   size=propVector[iFrom], prob=pSwap) 
            lpProp    <<- dbinom(nSwap, size=propVector[iFrom], prob=pSwap, log=TRUE) 
            propVector[iFrom] <- propVector[iFrom] - nSwap 
            propVector[iTo]   <- propVector[iTo]   + nSwap 
            pRevSwap   <- min(1, max(1, ENSwapMatrix[iTo,iFrom]) / (propVector[iTo] + nSwap)) 
            lpRev     <<- dbinom(nSwap, size=propVector[iTo], prob=pRevSwap, log=TRUE) 
            returnType(double(1)) 
            return(propVector) 
        },
        adaptiveProcedure = function(jump=logical(), iFrom=integer(), iTo=integer()) {
            NVector <- values(model,target) 
            timesRan[iFrom, iTo] <<- timesRan[iFrom, iTo] + 1
            if(jump)
                timesAccepted[iFrom, iTo] <<- timesAccepted[iFrom, iTo] + 1
            if (timesRan[iFrom, iTo] %% adaptInterval == 0) {
                totalAdapted[iFrom, iTo] <<- totalAdapted[iFrom, iTo] + 1
                accRate                   <- timesAccepted[iFrom, iTo] / timesRan[iFrom, iTo]
                AcceptRates[iFrom, iTo]  <<- accRate
                if (accRate > 0.5) {
                    ENSwapMatrix[iFrom, iTo] <<-
                        min(Ntotal,
                            ENSwapMatrix[iFrom,iTo] + ENSwapDeltaMatrix[iFrom, iTo] / totalAdapted[iFrom,iTo])
                } else {
                    ENSwapMatrix[iFrom, iTo] <<-
                        max(1,
                            ENSwapMatrix[iFrom,iTo] - ENSwapDeltaMatrix[iFrom,iTo] / totalAdapted[iFrom,iTo])
                } 
                if(accRate<RescaleThreshold[iFrom,iTo] | accRate>(1-RescaleThreshold[iFrom,iTo])) {
                    if (ENSwapMatrix[iFrom, iTo] > 1 & ENSwapMatrix[iFrom, iTo] < Ntotal) {
                        ScaleShifts[iFrom, iTo]       <<- ScaleShifts[iFrom, iTo] + 1 
                        ENSwapDeltaMatrix[iFrom, iTo] <<- min(NOverL, ENSwapDeltaMatrix[iFrom, iTo] * totalAdapted[iFrom,iTo] / 10)
                        ENSwapDeltaMatrix[iTo, iFrom] <<- ENSwapDeltaMatrix[iFrom, iTo] 
                        RescaleThreshold[iFrom,iTo]   <<- 0.2 * 0.95^ScaleShifts[iFrom, iTo]
                    }
                }
                if(ENSwapMatrix[iFrom, iTo] < 1)
                    ENSwapMatrix[iFrom, iTo] <<- 1                
                ENSwapMatrix[iTo,iFrom]   <<- ENSwapMatrix[iFrom,iTo]
                timesRan[iFrom, iTo]      <<- 0
                timesAccepted[iFrom, iTo] <<- 0
            }
        },
        reset = function() {
            timesRan          <<- Zeros
            AcceptRates       <<- Zeros
            ScaleShifts       <<- Zeros
            totalAdapted      <<- Zeros
            timesAccepted     <<- Zeros
            ENSwapMatrix      <<- Ones
            ENSwapDeltaMatrix <<- Ones
            RescaleThreshold  <<- 0.2 * Ones
        }
    ), where = getLoadingNamespace()
)

