########################################################################
##                        Diffusion Simulation                        ##
##                  Multivariate Normal Distribution                  ##
##                           Between Subject                          ##
########################################################################

## ---------------------------------------------------------------------
## General Description
## ---------------------------------------------------------------------

## Parameter distributions: Multivariate Normal Distribution (truncated)
## Analysis: ANOVA & GLM & GLMM
## Design: Between-Subjects

## This set of functions first simulate data using Drift Diffusion Model.
## The parameters are sampled from the multivariate normal distribution
## which takes correlations between parameters into consideration.
## The sampled parameters are then passed on to rdiffusion to generate data.
## The data are then analysed using ANOVA, GLM, or GLMM.

## ---------------------------------------------------------------------
## Preparation
## ---------------------------------------------------------------------

source("ddm_sampling_bs.R")
library("tibble")
library("car")
library("lme4")

## ---------------------------------------------------------------------
## Function for Data Simulation: One Group
## ---------------------------------------------------------------------

## One Iteration for Between-Subject Design ############################

## simulate_dt_bs(mu, sigma, size, pp, group)

## INPUT
##    mu: the vector of the mean of each parameter
## sigma: the covariance matrix which consists of
##  size: number of trials per participant
##    pp: number of participants per group
## group: the group label

## OUTPUT: Stimulate one group of data
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

## ---------------------------------------------------------------------
## Function for Data Analysis: One Factor
## ---------------------------------------------------------------------

## data_analysis_bs(mu, sigma, size, pp)
## Analyse one set of data with 'pp' participants and 'size' number of trials 

## INPUT
##    mu: the vector of the mean of each parameter
## sigma: the covariance matrix which consists of v, a, z, t0, sv, st0
##  size: number of trials per participant
##    pp: number of participants [per group]

## OUTPUT
##          aov_p: p-value generated by ANOVA
##          glm_p: p-value generated by GLM
##         glmm_p: p-value generated by GLMM
##     diff_props: difference in upper response proportion between groups
## mean_prop_real: mean upper response proportion for both groups  
##                 * converge with [prob in logistic simulation] for large samples
##        g1_prop: upper response proportion for group 1
##        g2_prop: upper response proportion for group 2
##           size: number of trials per participant
##          pp_g1: number of participants in group1, same as pp
##          pp_g2: number of participants in group2, same as pp
## DDM parameters: mean and sd for a, st0, sv, t0, v, z
##             df: each simulation of data
##     glm_object: glm model
##    glmm_object: glmm model

data_analysis_bs <- function(mu, sigma, size, pp) {
  # assume no main effect: same parameters distributions for pps in both groups
  
  # ------------ Simulate Data for Group 1 -------------
  g1 <- simulate_dt_bs(
    mu = mu,
    sigma = sigma,
    size = size,
    pp = pp,
    group = 1 # add group identifier
  )
  # calculate the number of pps
  pp_g1 <- pp
  
  # ------------ Simulate Data for Group 2 -------------
  g2 <- simulate_dt_bs(
    mu = mu, 
    sigma = sigma, 
    size = size,
    pp = pp,
    group = 2 # add group identifier
  )
  pp_g2 <- pp
  # pp_g2 <- as.numeric(as.character(summarise(g2, n_distinct(id))))
  
  # ------------ Prepare for analysis: resp_prop & n_trials -------------
  gg <- rbind(g1, g2) %>%
    # transfer upper/lower response to 1/0
    mutate(response = ifelse(response == "upper", 1, 0))
  
  # calculate the proportion of upper response
  mean_prop <- gg %>%
    group_by(group) %>%
    summarise(upper_props = mean(response))
  
  # upper response proportion for group1
  g1_prop <- mean_prop$upper_props[1] 
  # upper response proportion for group2
  g2_prop <- mean_prop$upper_props[2] 
  # overall mean of upper response proportion
  mean_prop_real <- mean(gg$response)
  
  # calculate the difference of the upper response proportion between groups
  diff_props <- g1_prop - g2_prop
  
  # add columns of resp_prop and n_trials
  gg_final <- gg %>% 
    group_by(id, group) %>% 
    summarise(resp_prop = mean(response), 
              n_trials = size) %>% 
    ungroup %>% 
    mutate(group = factor(group),
           id = factor(id))
  
  # -------------- Check if all response are the same --------------
  # when all responses are the same, set p-value to 1
  # so that the analysis model does not produce an error
  # particularly important when the number of trials is low]
  if (mean_prop_real == 1 || mean_prop_real == 0) {
    data <-
      tibble(
        aov_p = 1, 
        glm_p = 1,
        glmm_p = 1,
        diff_props = diff_props,
        mean_prop_real = mean_prop_real,
        g1_prop = g1_prop,
        g2_prop = g2_prop,
        size = size,
        pp_g1 = pp_g1,
        pp_g2 = pp_g2,
        
        # v, a, z, t0, sv, st0
        v_mean = mu[1],
        v_sd = sqrt(sigma[1, 1]),
        
        a_mean = mu[2],
        a_sd = sqrt(sigma[2, 2]),
        
        z_mean = mu[3],
        z_sd = sqrt(sigma[3, 3]),
        
        t0_mean = mu[4],
        t0_sd = sqrt(sigma[4, 4]),
        
        sv_mean = mu[5],
        sv_sd = sqrt(sigma[5, 5]),
        
        st0_mean = mu[6],
        st0_sd = sqrt(sigma[6, 6]),
        
        # save simulated data
        df = list(gg_final),
        glm_obj = list(1), # save glm object
        glmm_obj = list(1) # save glmm object
      )
  } else {
  
  # ------------ Analyse with ANOVA --------------
  try(aov_ss <- aov_ez(
    id = "id",
    dv = "resp_prop",
    data = gg_final,
    between = "group",
    fun_aggregate = mean
  )
  )
  if (exists("aov_ss")) {
    # return the p-value
    aov_p <- summary(aov_ss)[["Pr(>F)"]][[1]]
    aov_p
  } else {
    ## set p-value to 1
    aov_p = 1
    aov_p
  }
  
  # ------------ Analyse with GLM ---------------
  try(glm_gg <- glm(
    resp_prop ~ group,
    data = gg_final,
    weights = n_trials,
    family = binomial)
  )
  if (exists("glm_gg")) {
    # calculate type III ANOVA tables for various statistical models
    try(glm_anova <- car::Anova(glm_gg, type = 3))
    if (exists("glm_anova")) {
      glm_p <- glm_anova$`Pr(>Chisq)`
      glm <- tibble(
        glm_p = glm_p,
        glm_obj = list(glm_gg))
    } else {
      glm <- tibble(
        glm_p = 1,
        glm_obj = list(glm_gg))
      }
  } else {
    ## set p-values and so forth
    glm <- tibble(
      glm_p = 1,
      glm_obj = list(1)
    )
    glm
  }
  
  # ------------ Analyse with GLMM ---------------
  try(glmm_gg <- glmer(
    resp_prop ~ group + (1|id),
    data = gg_final,
    weights = n_trials,
    family = binomial
#    ,control = glmerControl(maxit = 1000)
  ))
  if (exists("glmm_gg")) {
    # calculate type III ANOVA tables for various statistical models
    try(glmm_anova <- car::Anova(glmm_gg, type = 3))
    if (exists("glmm_anova")) {
      glmm_p <- glmm_anova$`Pr(>Chisq)`[2]
      glmm <- tibble(
        glmm_p = glmm_p,
        glmm_obj = list(glmm_gg))
    } else {
      glmm <- tibble(
        glmm_p = 1,
        glmm_obj = list(glmm_gg))
      glmm
    }
  } else {
    ## set p-values and so forth
    glmm <- tibble(
      glmm_p = 1,
      glmm_obj = list(1)
    )
    glmm
  }
  
  # ------------ Return results of analysis & summary statistics -------------
  data <-
    tibble(
      aov_p = aov_p,
      glm_p = glm$glm_p,
      glmm_p = glmm$glmm_p,
      
      diff_props = diff_props,
      mean_prop_real = mean_prop_real,
      g1_prop = g1_prop,
      g2_prop = g2_prop,
      size = size,
      pp_g1 = pp_g1,
      pp_g2 = pp_g2,
      
      # v, a, z, t0, sv, st0
      v_mean = mu[1],
      v_sd = sqrt(sigma[1, 1]),
      
      a_mean = mu[2],
      a_sd = sqrt(sigma[2, 2]),
      
      z_mean = mu[3],
      z_sd = sqrt(sigma[3, 3]),
      
      t0_mean = mu[4],
      t0_sd = sqrt(sigma[4, 4]),
      
      sv_mean = mu[5],
      sv_sd = sqrt(sigma[5, 5]),
      
      st0_mean = mu[6],
      st0_sd = sqrt(sigma[6, 6]),
      
      # save simulated data
      df = list(gg_final),
      glm_obj = glm$glm_obj, # save glm object 
      glmm_obj = glmm$glmm_obj # save glmm object 
    )
  }
  return(data)
}

## ReRun ###############################################################

## reruns_bs(mu, sigma, size, pp, reruns)
## run for 'rerun' times of simulation
## simulate data with DDM with one BS factor with two levels / two groups
## return p-values from different analysis for this simulation

##     mu: the vector of the mean of each parameter
##  sigma: the covariance matrix which consists of
##   size: number of trials per participant
##     pp: number of participants [per group]
## reruns: times of rerun

reruns_bs <- function(mu, sigma, size, pp, reruns) {
  result <- map(1: reruns, 
                ~ data_analysis_bs(mu, sigma, size, pp))
  return(result)
}





