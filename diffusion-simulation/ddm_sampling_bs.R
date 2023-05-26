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

## ---------------------------------------------------------------------
## Function for data simulation
## ---------------------------------------------------------------------

## Simulating One Single Run for Between-Subject Design ################

## simulate_dt_bs(mu, sigma, size, pp, group)

## INPUT
##    mu: the vector of the mean of each parameter
## sigma: the covariance matrix which consists of
##  size: number of trials per participant
##    pp: number of participants per group
## group: the group label

## OUTPUT: One single simulation
## Simulated data of pp participants in a group (size trials per pp)
## using DDM parameters sampled from a mvt normal dis with mu and sigma

simulate_dt_bs <- function(mu, sigma, size, pp, group) {
  # params(mu, sigma): simulate n trials (n = size) for one pp
  #                    using one set of parameters sampled from 
  #                    the truncated multivariate normal distribution
  params <- function(mu, sigma) {
    # values: one set of parameters for one DDM model (each pp)
    values <- par_from_val(mu, sigma)
    trials <- rdiffusion(
      n = size,
      v = values$v,
      a = values$a,
      z = values$z * values$a, # relative scale
      t0 = values$t0,
      sv = values$sv,
      st0 = values$st0,
      stop_on_error = FALSE
      # return 0 if parameters values are outside the allowed range
    )
  }
  
  # repeat the simulation for all pps (n = pp)
  repeat
  {result <- map(1:pp, ~ params(mu, sigma)) %>%
    rbindlist(., idcol = TRUE) %>% # add id column
    mutate(group = rep(group))  # add group number 
  
  # remove those participants whose rdiffusion produce an error
  result <- result %>% 
    rename("id" = ".id") %>% 
    mutate(id = paste0(group, "_", id)) %>%
    group_by(id) %>%
    filter(mean(rt) != 0) %>% # by using this filter
    ungroup()
  
  # break until the valid simulation reaches the number of pps
  if (as.numeric(as.character(summarise(result, n_distinct(id)))) == pp) break
  }
  return(result)
}

