########################################################################
##                        Diffusion Simulation                        ##
##                  Multivariate Normal Distribution                  ##
##                           Between Subject                          ##
########################################################################

## ---------------------------------------------------------------------
## General Description
## ---------------------------------------------------------------------

## Parameter distributions: Multivariate Normal Distribution (truncated or not)
## Analysis: ANOVA & GLM & GLMM
## Design: Between-Subjects with 2 Factors

## This set of functions first simulate data using Drift Diffusion Model
## The parameter values are sampled from the multivariate normal distribution
## to take correlations between parameters into consideration
## The sampled parameters are then passed on to rdiffusion to generate data 
## The data is then analysed using ANOVA, GLM, or GLMM

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

## ---------------------------------------------------------------------
## Function for Data Simulation & Analysis: Two Factor
## ---------------------------------------------------------------------

## data_analysis_bs_2f(mu, sigma, size, pp)
## Analyse one set of data w 'pp' participants and 'size' number of trials
## with two BETWEEN-SUBJECT factors with interaction for 4 groups

## INPUT
##    mu: the vector of the mean of each parameter
## sigma: the covariance matrix which consists of a, st0, sv, t0, v, z
##  size: number of trials per participant
##    pp: number of participants [per group]

data_analysis_bs_2f <- function(mu, sigma, size, pp) {
  # assume no main effect: DDM parameter distributions are the same for all groups
  
  # ------------ Simulate Data for Groups 1 & 2 & 3 & 4-----------------
  # group1: a = 1, b = 1
  g1 <- simulate_dt_bs(
    mu = mu,
    sigma = sigma,
    size = size,
    pp = pp,
    group = 1 # add group identifier
  )
  g1 <- as.data.frame(g1) %>% 
    mutate(group = 1, a = 1, b = 1)
  
  # group2: a = 1, b = 2
  g2 <- simulate_dt_bs(
    mu = mu,
    sigma = sigma,
    size = size,
    pp = pp,
    group = 2 # add group identifier
  )
  g2 <- as.data.frame(g2) %>% 
    mutate(group = 2, a = 1, b = 2)
  
  # group3: a = 2, b = 1
  g3 <- simulate_dt_bs(
    mu = mu,
    sigma = sigma,
    size = size,
    pp = pp,
    group = 3 # add group identifier
  )
  g3 <- as.data.frame(g3) %>% 
    mutate(group = 3, a = 2, b = 1)
  
  # group4: a = 2, b = 2
  g4 <- simulate_dt_bs(
    mu = mu,
    sigma = sigma,
    size = size,
    pp = pp,
    group = 4 # add group identifier
  )
  g4 <- as.data.frame(g4) %>% 
    mutate(group = 4, a = 2, b = 2)
  
  # ------------ Prepare for analysis -------------
  gg <- rbind(g1, g2, g3, g4) %>%
    # transfer upper/lower response to 1/0
    mutate(response = ifelse(response == "upper", 1, 0)) %>% 
    mutate(group = factor(group),
           a = factor(a),
           b = factor(b),
           id = factor(id))
  
  # add columns of resp_prop and n_trials
  gg_final <- gg %>% 
    group_by(id, group) %>% 
    summarise(resp_prop = mean(response), 
              n_trials = size) %>% 
    ungroup() %>% 
    mutate(a = )
  # ???????????? how to set the group label
  
  # calculate the proportion of upper response
  mean_prop <- gg_final %>%
    group_by(group) %>%
    summarise(upper_props = mean(response))
  
  # upper response proportion for each group
  g1_prop <- mean_prop$upper_props[1]
  g2_prop <- mean_prop$upper_props[2]
  g3_prop <- mean_prop$upper_props[3]
  g4_prop <- mean_prop$upper_props[4]
  
  # overall mean of upper response proportion
  mean_prop_real <- mean(gg_final$response)
  
  # calculate the difference of conditional means
  # a: 1 vs 2
  diff_a <- (g1_prop+g2_prop)/2 - (g3_prop+g4_prop)/2
  # b: 1 vs 2
  diff_a <- (g1_prop+g3_prop)/2 - (g2_prop+g4_prop)/2
  
  # calculate the difference of the upper response proportion between groups
  diff_props <- g1_prop - g2_prop
  
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
        pp = pp,
        
        a_mean = mu[1],
        a_sd = sqrt(sigma[1, 1]),
        
        st0_mean = mu[2],
        st0_sd = sqrt(sigma[2, 2]),
        
        sv_mean = mu[3],
        sv_sd = sqrt(sigma[3, 3]),
        
        t0_mean = mu[4],
        t0_sd = sqrt(sigma[4, 4]),
        
        v_mean = mu[5],
        v_sd = sqrt(sigma[5, 5]),
        
        z_mean = mu[6],
        z_sd = sqrt(sigma[6, 6]),
        
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
      pp = pp,
      
      a_mean = mu[1],
      a_sd = sqrt(sigma[1, 1]),
      
      st0_mean = mu[2],
      st0_sd = sqrt(sigma[2, 2]),
      
      sv_mean = mu[3],
      sv_sd = sqrt(sigma[3, 3]),
      
      t0_mean = mu[4],
      t0_sd = sqrt(sigma[4, 4]),
      
      v_mean = mu[5],
      v_sd = sqrt(sigma[5, 5]),
      
      z_mean = mu[6],
      z_sd = sqrt(sigma[6, 6]),
      
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
  result <- rerun(reruns, 
                  data_analysis_bs(mu, sigma, size, pp))
  return(result)
}





