########################################################################
##                        Diffusion Simulation                        ##
##                  Multivariate Normal Distribution                  ##
##              Parameter Sampling & Single Stimulation               ##
##                           Between Subject                          ##
########################################################################

library("gmm")
library("tmvtnorm")
library("rtdists")

## ---------------------------------------------------------------------
## Description
## ---------------------------------------------------------------------

## Sampling the parameter values from $multivariate$ distribution.
## It needs 'optimization algorithm', but just uses the means and
## covariance matrix of parameters for Ratcliff diffusion model.

## DDM - drift diffusion model - Ratcliff diffusion model

## ---------------------------------------------------------------------
## Function for parameter sampling
## ---------------------------------------------------------------------

## par_from_val(mu, sigma)

## INPUT
##    mu: the vector of the mean of each parameter
##        c(v a z t0 sv st0)
## sigma: the covariance matrix which consists of
##        (1) the variance of each parameter along the diagonal
##        (2) the covariance between each pair of parameters

## OUTPUT
## A vector of parameter values sampled from parameter distributions

## Input mu and cor matrix consist of six parameters
##   v: drift rate/quality of information accumulation, 10 > ~ > -10
##   a: threshold separation, > 0
##   z: relative starting point, 0.5*a is unbiased decision, 1 > ~ > 0
##  t0: non-decision time, > 0
##  sv: intertrial variability of drift rate, 10 > ~ > 0
## st0: intertrial variability of non-decision time, > 0
##      (Date range is from observation of the original data)

par_from_val <- function(mu, sigma) {
  # generate random numbers from the truncated mv normal distribution
  # n: times of random sampling, here just return one set of parameters for DDM
  gen <- rtmvnorm(n = 1, mean = mu, sigma = sigma,
                  #       c(  v,   a,  z,  t0, sv, st0)
                  lower = c(-10,   0,  0,   0,  0,   0),
                  upper = c( 10, Inf,  1, Inf, 10, Inf),
                  H = NULL)
  gen <- as.data.frame(gen)
  gen <- gen %>% 
    rename(
      v = V1,
      a = V2,
      z = V3,
      t0 = V4,
      sv = V5,
      st0 = V6)
  gen
}
