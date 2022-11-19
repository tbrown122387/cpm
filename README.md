# cPseudoMaRg

<!-- badges: start -->
  

[![DOI](https://zenodo.org/badge/381811804.svg)](https://zenodo.org/badge/latestdoi/381811804)
[![R-CMD-check](https://github.com/tbrown122387/cpm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tbrown122387/cpm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

An implementation of the [Correlated Pseudo-Marginal Sampler](https://arxiv.org/abs/1511.04992).

## Install

Install from CRAN by typing 

```
install.packages("cPseudoMaRg")
```

in an R console. Alternatively, install from Github by typing 

```
devtools::install_github("tbrown122387/cpm")
```

## Example

Another Random Effects Model that mimics the example in the above paper. They estimate a mean parameter, whereas the unknown parameters here are variance parameters. Also, this model's likelihood is nonidentifiable. 

```
# y | x, theta ~ Normal(x, SSy)
# x | theta ~ Normal(0, SSx)
# theta = (SSy + SSx, SS_x)
# p(theta | y) propto p(y | theta)p(theta)
# approx p(y | theta) with mean( p(y | xi, theta)  ) where xi ~ p(xi | theta)

# real data
realxVar <- .2
realyVar <- .3
realTheta1 <- realxVar + realyVar
realTheta2 <- realxVar
realParams <- c(realTheta1, realTheta2)
numObs <- 10
realX <- rnorm(numObs, mean = 0, sd = sqrt(realxVar))
realY <- rnorm(numObs, mean = realX, sd = sqrt(realyVar))

# tuning params
numImportanceSamps <- 1000
numMCMCIters <- 10000
randomWalkScale <- 1.5
recordEveryTh <- 1
myLLApproxEval <- function(y, thetaProposal, uProposal){
  if( (thetaProposal[1] > thetaProposal[2]) & (all(thetaProposal > 0))){
    xSamps <- uProposal*sqrt(thetaProposal[2])
    logCondLikes <- sapply(xSamps,
                           function(xsamp) {
                             sum(dnorm(y,
                                       xsamp,
                                       sqrt(thetaProposal[1] - thetaProposal[2]),
                                       log = T)) })
    m <- max(logCondLikes)
    log(sum(exp(logCondLikes - m))) + m - log(length(y))
  }else{
    -Inf
  }
}
sampler <- makeCPMSampler(
  paramKernSamp = function(params){
    return(params + rnorm(2)*randomWalkScale)
  },
  logParamKernEval = function(oldTheta, newTheta){
    dnorm(newTheta[1], oldTheta[1], sd = randomWalkScale, log = TRUE)
         + dnorm(newTheta[2], oldTheta[2], sd = randomWalkScale, log = TRUE)
  },
  logPriorEval = function(theta){
    if( (theta[1] > theta[2]) & all(theta > 0)){
      0
    }else{
      -Inf
    }
  },
  logLikeApproxEval = myLLApproxEval,
  realY, numImportanceSamps, numMCMCIters, .99, recordEveryTh
)
res <- sampler(realParams)

# look at output
print(res)
plot(res)
```


