########################################################################
##                      Ratcliff Diffusion Model                      ##
##                         Parameter Estimate                         ##
########################################################################

library("rtdists")
library("dplyr")
library("purrr")
library("tidyr")
library("ggplot2")
library("parallel")

## Principle

## Maximum likelihood estimation (MLE)
## searches for those parameter values that maximize the likelihood function 
## given a fixed data set; use probability density function (PDF) for fitting;
## assume data are fixed (not parameters - then the PDF integrates to 1); MLE 
## means to search for the parameters that maximize the PDF given fixed data

## In this case, PDF is called the likelihood function as it looses the property 
## to integrate to 1 across the search space

## Description

## We fit the Ratcliff diffusion model to all trials per participant
## using maximum likelihood estimation (MLE)

## Basic Parameters
##   a: an (ntp) & at (tp)
##   v: vng (ntp + gain) & vnl (ntp + loss) & 
##      vtg (tp + gain) & vtl (tp + loss)
##  t0: t0n (ntp) & t0t (tp)
##   z: z for all
##  sv: sv for all
## st0: st0 for all

## parameters (an, at, vng, vnl, vtg, vtl, t0n, t0t, z, sv, st0)

## ---------------------------------------------------------------------
## Load the data
## ---------------------------------------------------------------------

load("./data_final.rda")

length(unique(data_final$id))

## variable: response and id
data_final <- data_final %>% 
  mutate(response_gamble = factor(response_gamble,
                                  levels = c(0, 1),
                                  labels = c("lower", "upper")),
         id = factor(id, labels = 1:124))
data_final <- rename(data_final, "response" = "response_gamble")
data_final <- data_final %>% 
  mutate(id = as.double(id))

glimpse(data_final)
unique(data_final$id)
length(unique(data_final$id)) # 124

## see the distributions of rt in different conditions
data_final %>% 
  ggplot(aes(x = rt)) +
  geom_histogram(binwidth = 0.1, boundary = 0) +
  facet_grid(rows = vars(frame), cols = vars(time_pressure))

## ---------------------------------------------------------------------
## Multiple starting values - avoid Local Optima
## ---------------------------------------------------------------------

## To avoid the local optimum instead of the global optimum
## we usually run our optimization algorithm several times
## but the optimization algorithms here are deterministic
## meaning that same start values lead to identical results
## so run the algorithm with DIFFERENT starting values

## get_start_values()
## returns random start values of 
## (an, at, vng, vnl, vtg, vtl, t0n, t0t, z, sv, st0)

# runif(n, min, max) generates n random deviates in a uniform distribution
# rnorm(n, mean, sd) generates n random deviates in a normal distribution

ntp <- data_final %>% filter(time_pressure == "ntp")
min(ntp$rt) # [1] 0.204063
tp <- data_final %>% filter(time_pressure == "tp")
min(tp$rt) # [1] 0.200433

get_start_values <- function() {
  c(
    an = runif(1, 1, 4),
    at = runif(1, 1, 4),
    vng = rnorm(1, 0, 2), # only v is from norm dis
    vnl = rnorm(1, 0, 2),
    vtg = rnorm(1, 0, 2),
    vtl = rnorm(1, 0, 2),
    t0n = runif(1, 0, 0.204), # always smaller than the smallest rt
    t0t = runif(1, 0, 0.2), 
    z = runif(1, 0.4, 0.6), # use the relative scale; around 0.5
    sv = runif(1, 0, 2), # always positive
    st0 = runif(1, 0, 0.5) # always positive
  )
}

get_start_values()

## ---------------------------------------------------------------------
## Log-Likelihood Function
## ---------------------------------------------------------------------

## ll_diffusion(pars, df)
## return a positive value which is the negative sum of log-likelihoods 
## of certain data set given a set of parameters

## Input
##  pars: a set of starting values of parameters
## df_pp: data of one participant

## Output
## densities in four conditions

ll_diffusion <- function(pars, df_pp) {
  ## densities is a scalar value - the negative sum of log-likelihoods
  
  ## split the data frame into four conditions for each participant
  df_ng <- df_pp %>% filter(time_pressure == "ntp" & frame == "gain")
  df_nl <- df_pp %>% filter(time_pressure == "ntp" & frame == "loss")
  df_tg <- df_pp %>% filter(time_pressure == "tp" & frame == "gain")
  df_tl <- df_pp %>% filter(time_pressure == "tp" & frame == "loss")
  
  ## calculate densities for each condition
  ## densities_ng = no time pressure & gain
  densities_ng <- ddiffusion(df_ng$rt, response = df_ng$response, 
                             a = pars["an"],
                             v = pars["vng"], 
                             t0 = pars["t0n"],
                             z = pars["z"]*pars["an"],
                             sv = pars["sv"],
                             st0 = pars["st0"])
  ## densities_nl = no time pressure & loss
  densities_nl <- ddiffusion(df_nl$rt, response = df_nl$response, 
                             a = pars["an"],
                             v = pars["vnl"], 
                             t0 = pars["t0n"],
                             z = pars["z"]*pars["an"],
                             sv = pars["sv"],
                             st0 = pars["st0"])
  ## densities_tg = time pressure & gain
  densities_tg <- ddiffusion(df_tg$rt, response = df_tg$response, 
                             a = pars["at"],
                             v = pars["vtg"], 
                             t0 = pars["t0t"],
                             z = pars["z"]*pars["at"],
                             sv = pars["sv"],
                             st0 = pars["st0"])
  ## densities_tl = time pressure & loss
  densities_tl <- ddiffusion(df_tg$rt, response = df_tg$response, 
                             a = pars["at"],
                             v = pars["vtl"], 
                             t0 = pars["t0t"],
                             z = pars["z"]*pars["at"],
                             sv = pars["sv"],
                             st0 = pars["st0"])
  
  ## combine all densities
  densities <- c(densities_ng, densities_nl, 
                 densities_tg, densities_tl)
  
  ## if an individual densities is 0, returns 1e6
  ## or will be error from -log(0)
  ## e.g., 't0' is larger than an observed RT
  if (any(densities == 0)) return(1e6)
  
  return(-sum(log(densities)))
}

## ---------------------------------------------------------------------
## Maximum likelihood estimate
## ---------------------------------------------------------------------

## Use an optimization algorithm using numerical methods (gradient descent) 
## to find the maximum likelihood estimate for the current data

## nlminb(start, objective, ..., lower = -Inf, upper = Inf)
## do unconstrained and box-constrained optimization using PORT routines to
## minimize the negative log-likelihood by changing initial parameter values
## iteratively until a pre-specified convergence criterion is reached (i.e., 
## when changing parameters doesn't decrease negative sum of log-likelihoods 
## any more); it allows specifying parameters bounds (e.g., a and t0 > 0)

##     start: a set of start values for the search (e.g., par1)
## objective: the "objective" function that should be minimized (ll_diffusion)
##       ...: further arguments to be supplied to objective.
##     lower: vectors of lower bounds
##     upper: vectors of upper bounds

## Generally
## split the data by each participant
## loop optimization over each participant on four conditions

## prepare #############################################################
n_participants <- length(unique(data_final$id))
df_all <- vector("list", n_participants)
res_all <- vector("list", n_participants)

## rerun function
reruns_fitting <- function(reruns, df_pp){
  result <- map(1:reruns, 
                ~ nlminb(start = get_start_values(), 
                         objective = ll_diffusion, 
                         df_pp = df_pp,
                         # c(an, at, vng, vnl, vtg, vtl, t0n, t0t, z, sv, st0)
                         lower = c(0, 0, -Inf, -Inf, -Inf, -Inf, 0, 0, 0, 0, 0),
                         upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 1, Inf, Inf)))
  return(result)
}

## abstract the data for each participant
for (i in 1:n_participants) {
  df_all[[i]] <- data_final %>% filter(id == i)
}

## run 20 fitting for each pp ##########################################
res_all <- mclapply(df_all, reruns_fitting, reruns = 20, mc.cores = 12)

save(res_all, file = "res_all.rda")

## ---------------------------------------------------------------------
## Get parameter values, logLik, and convergence for each fit
## ---------------------------------------------------------------------

## create an empty tibble for parameters ###############################
pars_all <- tibble(id = NA, fit = NA,
                   an = NA, at = NA, 
                   vng = NA, vnl = NA, vtg = NA, vtl = NA,
                   t0n = NA, t0t = NA,
                   z = NA, sv = NA, st0 = NA,
                   logLik = NA, convergence = NA)
pars_all

## get estimated parameter values from each fit for each participant 
for (j in 1:n_participants) {
  for (k in 1:20) {
    ## get estimated parameter values, logLik, convergence for each fit
    id = j
    fit = k
    pars <- res_all[[j]][[k]][["par"]]
    logLik <- -res_all[[j]][[k]][["objective"]]
    convergence <- res_all[[j]][[k]][["convergence"]]
    new_fit <- c(id, fit, pars, logLik, convergence)
    pars_all <- rbind(pars_all, new_fit)
  }
}

## remove the first row
pars_all <- pars_all[-1, ]
pars_all <- as.data.frame(pars_all)
save(pars_all, file = "pars_all.rda")

## ---------------------------------------------------------------------
## Pick "correct" fit for each participant
## ---------------------------------------------------------------------

glimpse(pars_all)

mle <- tibble(id = NA, fit = NA,
              an = NA, at = NA, 
              vng = NA, vnl = NA, vtg = NA, vtl = NA,
              t0n = NA, t0t = NA,
              z = NA, sv = NA, st0 = NA,
              logLik = NA, convergence = NA)

for (i in 1:n_participants) {
  tmp <- pars_all %>% filter(id == i)
  new_mle <- tmp %>% slice(which.max(logLik))
  mle <- rbind(mle, new_mle)
}

## remove the first row
mle <- mle[-1, ]
mle

## check the identifiability
mle_identifiable <- tibble(id = NA, fit = NA,
                           an = NA, at = NA, 
                           vng = NA, vnl = NA, vtg = NA, vtl = NA,
                           t0n = NA, t0t = NA,
                           z = NA, sv = NA, st0 = NA,
                           logLik = NA, convergence = NA)

for (i in 1:n_participants) {
  ## tmp: 20 fits for each participant
  fit20_tmp <- pars_all %>% filter(id == i)
  
  ## get the mle for each participant
  mle_tmp <- mle %>% filter(id == i)
  mle_tmp <- mle_tmp %>% as.data.frame()
  
  ## get 2nd mle for each participant
  fit19_tmp <- fit20_tmp %>% filter(fit != mle_tmp$fit)
  mle2_tmp <- fit19_tmp %>% slice(which.max(logLik))
  
  if (abs(mle_tmp$logLik - mle2_tmp$logLik) < 0.0001) {
    identifiability <- as.vector(abs(mle_tmp[, 3:13] - mle2_tmp[, 3:13]) > 0.01)
    mle_tmp[, 3:13][identifiability] <- NA
  }
  
  ## compare mle and 2nd mle
  mle_identifiable <- rbind(mle_identifiable, mle_tmp)
}

## remove the first row
mle_identifiable <- mle_identifiable[-1, ]
mle_identifiable

## ---------------------------------------------------------------------
## Calculate the mean and covariance matrix
## ---------------------------------------------------------------------

## parameters (an, at, vng, vnl, vtg, vtl, t0n, t0t, z, sv, st0)

mle_identifiable %>% 
  select(c(vnl, an, z, t0n, sv, st0)) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(vars(name), scales = "free") + 
  labs(title = "all")

## We only include parameters in ntp condition bc we don't need time pressure.
## But we use data from all the trials to estimate parameters since we believe
## those parameters that are consistent across all conditions will be more 
## representative in this way (not very sure).
## And we also use ntp trials when fitting but the problematic outliers are not
## significantly reduced. So we stick to doing the fit with all trials - more
## trials, more accurate model fitting and MLE
## And we only used vnl instead of collasping vnl & vng bc they are hugely
## different from each other (positive vs. negative), so we used vnl whose 
## distribution (histogram) is more on the postive side.

## inclusion range: according to observations of the univariate distributions
## vnl < 5
##  an < 3
##   z [0.35, 0.75]
## t0n [0.05, 0.6]
##  sv < 2.5
## st0 < 0.6

mle_identifiable_normal <- mle_identifiable %>% 
  select(c(vnl, an, z, t0n, sv, st0)) %>% 
  mutate(
    vnl = ifelse(vnl < 5, vnl, NA),
    an = ifelse(an < 3, an, NA),
    z = ifelse(z > 0.35 & z < 0.75, z, NA),
    t0n = ifelse(t0n > 0.05 & t0n < 0.6, t0n, NA),
    sv = ifelse(sv < 2.5, sv, NA),
    st0 = ifelse(st0 < 0.6, st0, NA)
  )

mu <- mle_identifiable_normal %>% 
  summarise(across(everything(), ~mean(., na.rm = TRUE)))

mu
# A tibble: 1 × 6
#     vnl    an     z   t0n    sv   st0
#   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 0.975  2.01 0.542 0.280 0.639 0.279

cov <- mle_identifiable_normal %>% 
  cov(use = "pairwise.complete.obs")

cor <- mle_identifiable_normal %>% 
  cor(use = "pairwise.complete.obs")

save(mu, file = "mu.rda")
save(cov, file = "cov.rda")
save(cor, file = "cor.rda")
