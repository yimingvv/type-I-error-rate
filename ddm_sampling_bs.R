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

## Sampling the parameter values
## Fitting multivariate distribution needs no 'optimization algorithm'
## and the covariance matrix which consists of the variance of each parameter along the diagonal and the covariance between every possible pair of variables in other positions of the matrix

## ---------------------------------------------------------------------
## Function for parameter sampling
## ---------------------------------------------------------------------

## par_from_val(mu, sigma)

## INPUT
##    mu: the vector of the mean of each parameter
## sigma: the covariance matrix which consists of
##        (1) the variance of each parameter along the diagonal
##        (2) the covariance between each pair of parameters

## OUTPUT
## A vector of parameter values sampled from a best-fit distribution

## Input sigma matrix consist of six parameters
##   a: threshold separation, > 0
## st0: intertrial variability of non-decision time, > 0
##  sv: intertrial variability of drift rate, 10 > ~ > 0
##  t0: non-decision time, > 0
##   v: drift rate/quality of information accumulation, 10 > ~ > -10
##   z: starting point, 0.5*a is unbiased decision, 1 > ~ > 0

par_from_val <- function(mu, sigma) {
  # generate random numbers from the truncated mv normal distribution
  # n: times of random sampling, here just return one set of parameters for DDM
  gen <- rtmvnorm(n = 1, mean = mu, sigma = sigma,
                  # WHY sv ~ [0, 10]
                  # WHY v ~ [-10, 10]
                  # WHY z ~ [0, 1]
                  # from observation of the original data
                  lower = c(0,     0,  0,   0, -10, 0),
                  upper = c(Inf, Inf, 10, Inf,  10, 1),
                  H = NULL)
  gen <- as.data.frame(gen)
  gen <- gen %>% 
    rename(
      a = V1,
      st0 = V2,
      sv = V3,
      t0 = V4,
      v = V5,
      z = V6)
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
## using the DDM parameters sampled from a mvt nd with mu and sigma

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
      t0 = values$t0,
      sv = values$sv,
      st0 = values$st0,
      z = values$z,
      stop_on_error = FALSE
      # return 0 if parameters values are outside the allowed range
    )
  }
  
  # repeat the simulation for all pps (n = pp)
  repeat
  {result <- rerun(pp, params(mu, sigma)) %>%
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
