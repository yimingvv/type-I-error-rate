########################################################################
##                         Logistic Simulation                        ##
##                     with Individual Differences                    ##
##                              One Factor                            ##
########################################################################

library("dplyr")
library("tidyr")
library("tidyverse")
library("afex")
library("emmeans")
library("purrr")
library("stats")
library("lme4")

## -----------------------------------------------------------------------
## BETWEEN-SUBJECT DEISGN
## -----------------------------------------------------------------------

########################################################################

## logistic_bs_id_one_run(n, size, prob, sd = 0.5)
## one single run
## simulate data with logistic regression
## with one BETWEEN-SUBJECT factor with two levels / two groups
## and return p-values from two analysis for this simulation

## n: number of participants per group
## size: number of trials per participant
## prob: probability of success on each trial
## sd: standard deviation of individual variances around prob, default is 0.5

logistic_bs_id_one_run <- function(n, size, prob, sd = 0.5) {
  # assume no main effect: base prob is the same for all pps in both groups
  # assume individual variances around the base prob for each pp
  
  # ------------ Simulate Data for Groups 1 & 2 -----------------
  prob_g1 <- plogis(rnorm(n, qlogis(prob), sd))
  g1 <- rbinom(n, size, prob_g1)
  prob_g2 <- plogis(rnorm(n, qlogis(prob), sd))
  g2 <- rbinom(n, size, prob_g2)
  
  # ------------ Prepare for analysis: resp_prop & n_trials -----------------
  # combine by columns and rename the column [g1 g2]
  gg <- cbind(g1, g2)
  gg <- as.data.frame(gg)
  
  # calculate the mean response rate for each group [conditional mean]
  g1_prop <- mean(gg$g1)/size
  g2_prop <- mean(gg$g2)/size
  
  # calculate the difference of the mean response rate
  diff_props <- g1_prop - g2_prop
  # calculate the mean response rate for all participants [overall mean]
  mean_prop_real <- (g1_prop + g2_prop)/2
  
  gg <- data.frame(
    response = c(g1, g2),
    group = rep(c("g1", "g2"), each = n),
    id = seq.int(n*2)
  )

  # calculate response proportion from the sum of correct trials in 'size' trials
  gg_final <- gg %>% 
    group_by(id, group) %>% 
    summarise(resp_prop = response/size,
              n_trials = size) %>%
    ungroup %>% 
    mutate(group = ifelse(group == "g1", 1, 2)) %>%
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
        n_pp = n,
        size = size,
        prob = prob,
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
    # summary(aov_ss)[["Pr(>F)"]] returns the same output
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
    glm
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
#   ,control = glmerControl(maxit = 1000)
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
      n_pp = n,
      size = size,
      prob = prob,
      df = list(gg_final),
      glm_obj = glm$glm_obj, # save glm object for the convergence problem
      glmm_obj = glmm$glmm_obj # save glmm object as na
    )
  }
  return(data)
}

########################################################################

## rerun_logistic_bs_id(n, size, prob, sd = 0.5, reruns)
## run for SEVERAL times 'reruns'
## simulate data with logistic regression considering individual differences
## with one BETWEEN-SUBJECT factor with two levels / two groups
## and return p-values from two analysis for this simulation

##     n: number of observations, AKA participants
##  size: number of trials per participant
##  prob: probability of success on each trial
##    sd: standard deviation of individual variances around prob, default is 0.5
## rerun: times of rerunning

rerun_logistic_bs_id <- function(n, size, prob, sd, reruns) {
  result <- rerun(reruns, 
                  logistic_bs_id_one_run(n, size, prob, sd))
  return(result)
}








