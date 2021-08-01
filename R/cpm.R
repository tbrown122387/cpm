#' checks if a log-density evaluation is not a valid number
#' 
#' @param num evaluation of a log-density
#' @return TRUE or FALSE
#' @export
#' @examples
#' isBadNum(NaN)
isBadNum <- function(num){
  !is.numeric(num) | is.na(num) | is.nan(num)
}


#' correlated pseudo-marginal: generates functions that output a big vector
#' 
#' @param paramKernSamp function(theta) -> theta proposal
#' @param logParamKernEval function(oldTheta, newTheta) -> logDensity.
#' @param logPriorEval function(theta) -> logDensity.
#' @param logLikeApproxEval function(y, thetaProposal, uProposal) -> logApproxDensity.
#' @param yData the observed data
#' @param numU integer number of u samples
#' @param numIters integer number of MCMC iterations
#' @param rho correlation tuning parameter (-1,1)
#' @param storeEvery increase this integer if you want to use thinning
#' @return vector of theta samples
#' @export
#' @examples
#' 
#' # sim data
#' realTheta1 <- .2 + .3
#' realTheta2 <- .2
#' realParams <- c(realTheta1, realTheta2)
#' numObs <- 10
#' realX <- rnorm(numObs, mean = 0, sd = sqrt(realTheta2))
#' realY <- rnorm(numObs, mean = realX, sd = sqrt(realTheta1 - realTheta2))
#' # tuning params
#' numImportanceSamps <- 1000
#' numMCMCIters <- 1000
#' randomWalkScale <- 1.5
#' recordEveryTh <- 1
#' sampler <- makeCPMSampler(
#'  paramKernSamp = function(params){
#'    return(params + rnorm(2)*randomWalkScale)
#'  },
#'  logParamKernEval = function(oldTheta, newTheta){
#'    dnorm(newTheta[1], oldTheta[1], sd = randomWalkScale, log = TRUE)
#'    + dnorm(newTheta[2], oldTheta[2], sd = randomWalkScale, log = TRUE)
#'  },
#'  logPriorEval = function(theta){
#'    if( (theta[1] > theta[2]) & all(theta > 0)){
#'      0
#'    }else{
#'      -Inf
#'    }
#'  },
#'  logLikeApproxEval = function(y, thetaProposal, uProposal){
#'    if( (thetaProposal[1] > thetaProposal[2]) & (all(thetaProposal > 0))){
#'      xSamps <- uProposal*sqrt(thetaProposal[2])
#'      logCondLikes <- sapply(xSamps,
#'                            function(xsamp) {
#'                              sum(dnorm(y,
#'                                        xsamp,
#'                                        sqrt(thetaProposal[1] - thetaProposal[2]),
#'                                        log = TRUE)) })
#'      m <- max(logCondLikes)
#'      log(sum(exp(logCondLikes - m))) + m - log(length(y))
#'    }else{
#'      -Inf
#'    }
#'  },
#'  realY, numImportanceSamps, numMCMCIters, .99, recordEveryTh
#')
#'res <- sampler(realParams)
makeCPMSampler <- function(paramKernSamp, logParamKernEval, 
                           logPriorEval, logLikeApproxEval,
                           yData, numU, numIters, 
                           rho = .99, storeEvery = 1){
  # checks
  stopifnot(typeof(paramKernSamp) == "closure")
  stopifnot(typeof(logParamKernEval) == "closure")
  stopifnot(typeof(logPriorEval) == "closure")
  stopifnot(typeof(logLikeApproxEval) == "closure")
  stopifnot(-1 < rho & rho < 1)
  
  # state
  U <- stats::rnorm(numU)
  theta <- vector(mode = "numeric", length = 0L)
  y <- yData
  logLikeApprox <- -Inf # set later
  logPrior <- -Inf # set later
  
  # function to be returned
  function(initParams){
    
    theta <<- initParams
    thetaSamps <- vector(mode = "list", 
                         length = numIters %/% storeEvery)
    numAccepts <- 0
    for(i in 1:numIters){
      prevIndex <- (i - 2 + storeEvery) %/% storeEvery
      thisIndex <- ( (i-1)  %/% storeEvery + 1)

      if(i > 1){
        
        thetaProposal <- paramKernSamp(thetaSamps[[prevIndex]])
        uProposal <- rho * U + sqrt(1 - rho^2) * stats::rnorm(numU)
        
        # hastings ratio
        propLogLikeEval <- logLikeApproxEval(y, thetaProposal, uProposal)
        propLogPriorEval <- logPriorEval(thetaProposal)
        forwardLogKern <- logParamKernEval(thetaProposal, theta)
        backwardLogKern <- logParamKernEval(theta, thetaProposal)
        stopifnot(!isBadNum(propLogLikeEval))
        stopifnot(!isBadNum(propLogPriorEval))
        stopifnot(!isBadNum(forwardLogKern))
        stopifnot(!isBadNum(backwardLogKern))
        
        logRatio <- propLogLikeEval - logLikeApprox
                  + propLogPriorEval - logPrior
                  + backwardLogKern - forwardLogKern
        
        accept <- log(stats::runif(1)) < logRatio
        if(accept){
          U <<- uProposal
          theta <<- thetaProposal
          logLikeApprox <<- propLogLikeEval
          logPrior <<- propLogPriorEval
          numAccepts <- numAccepts + 1
        }
      }else{ 
        logLikeApprox <<- logLikeApproxEval(y, theta, U)
        logPrior <<- logPriorEval(theta)
      }
      
      # record
      if((i-1) %% storeEvery == 0)
        thetaSamps[[thisIndex]] <- theta
    }
    
    stuff <- list(samples = thetaSamps, numAccepts = numAccepts, numIters = numIters)
    class(stuff) <- "cpmResults"
    stuff
  }  
}


#' calculates the posterior mean point estimate
#' 
#' @param x a cpmResults object
#' @param ... arguments to be passed to or from methods.
#' @return a vector of parameter estimates (posterior mean)
#' @export
mean.cpmResults <- function(x, ...){
  colMeans(do.call(rbind, x$samples))
}


#' prints a cpmResults object
#' 
#' @param x a cpmResults object
#' @param ... arguments to be passed to or from methods.
#' @return the same cpmResults object
#' @export
print.cpmResults <- function(x, ...){
  cat(
    "CPM Results: at a glance...\n",
    x$numIters, " iterations performed\n",
    length(x$samples), " samples retained\n",
    "posterior means (in the supplied parameterization): ", mean(x), "\n",
    "acceptance rate: ", x$numAccepts / x$numIters)
  invisible(x)
}

#' plots a cpmResults object
#' 
#' @param x a cpmResults object
#' @param ... arguments to be passed to or from methods.
#' @export
plot.cpmResults <- function(x, ...){
  graphics::pairs(do.call(rbind, x$samples), 
        labels = paste0("theta", 1:length(x$samples[[1]])))
}

