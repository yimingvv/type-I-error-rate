########################################################################
##                 Type I Error Rate: ANOVA, GLM & GLMM               ##
##                  Data Simulated from Logistic Model                ##
##                     with Individual Differences                    ##
##                           Between Subject                          ##
########################################################################

## A Thorough Analysis of Type I Error Rate

## ---------------------------------------------------------------------
## Prepare
## ---------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(binom)
library(purrr)
library(lme4)
library(graphics)
library(ggthemes)
theme_set(theme_bw(base_size = 12))


## -----------------------------------------------------------------------
## BS100
## -----------------------------------------------------------------------

## load data
load(file = "./results/log_bs_30pp100tr_1000runs.rda")
log100_bs <- simulated_data

########################################################################
## Set check labels for convergence problems ###########################
########################################################################

## CASE 1: mean_prop_real == 1 #########################################
log100_bs$manual_ch1 <- ifelse(log100_bs$mean_prop_real == 1, 
                                TRUE, FALSE)

## CASE 2: Inconvergence in GLM ########################################
log100_bs$glm_ch1 <- map_lgl(log100_bs$glm_obj, ~ .$converged == FALSE)

## CASE 3: Inconvergence in GLMM #######################################
## LABEL 1
log100_bs$glmm_ch1 <- map_lgl(log100_bs$glmm_obj,
                               ~ isTRUE(.@optinfo$conv$lme4$code == -1))
## LABEL 2
log100_bs$glmm_ch2 <- map_lgl(log100_bs$glmm_obj, isSingular)

## Release the workspace
log100_bs$glm_obj = NULL
log100_bs$glmm_obj = NULL

########################################################################
## Create type I error rate tibble #####################################
########################################################################

# ----------------- Overall ----------------- 

aov_error <- log100_bs %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- log100_bs %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- log100_bs %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glmm_error_Conver <- log100_bs %>% 
  filter(glmm_ch1 == FALSE & glmm_ch2 == FALSE) %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_log100_bs_all <- 
  bind_rows(aov_error, glm_error,
            glmm_error, glmm_error_Conver)

# save the file
save(error_log100_bs_all, file = "./results/error_log100_bs_all.rda")

# plot
p_log100_bs_all <- 
  ggplot(error_log100_bs_all,
         aes(x = analysis, y = mean)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.03)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "probability", y = "Type I Error Rate",
       title = "Logistic Simulation bs100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed") 
p_log100_bs_all

# ----------------- Group by prob ----------------- 

aov_error <- log100_bs %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- log100_bs %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- bs100_1000 %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glmm_error_Conver <- log100_bs %>% 
  filter(glmm_ch1 == FALSE & glmm_ch2 == FALSE) %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_log100_bs_prob <- 
  bind_rows(aov_error, glm_error,
            glmm_error, glmm_error_Conver)

# save the file
save(error_log100_bs_prob, file = "./results/error_log100_bs_prob.rda")

# plot
p_log100_bs_prob <- 
  ggplot(error_log100_bs_prob,
         aes(x = prob, y = mean, shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.03)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "probability", y = "Type I Error Rate",
       title = "Logistic Simulation bs100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed") 
p_log100_bs_prob

## -----------------------------------------------------------------------
## BS1
## -----------------------------------------------------------------------

load(file = "./results/log_bs_30pp1tr_1000runs.rda")
log1_bs <- bs1_1000

########################################################################
## Set check labels for convergence problems ###########################
########################################################################

## CASE 1: mean_prop_real == 1 #########################################
log1_bs$manual_ch1 <- ifelse(log1_bs$mean_prop_real == 1, TRUE, FALSE)

## CASE 2: Inconvergence in GLM ########################################
log1_bs$glm_ch1 <- TRUE
log1_bs$glm_ch1[!log1_bs$manual_ch1] <- 
  map_lgl(log1_bs$glm_obj[!log1_bs$manual_ch1], ~ !.$converged)

## CASE 3: Inconvergence in GLMM #######################################
## Recalculate glmm_obj
# new_fit = TRUE needs to recalculate glmm
log1_bs <- log1_bs %>% 
  mutate(new_fit = if_else(
    # mean_prop_real != 1 but glmm_obj = 1
    map_lgl(glmm_obj, is.numeric) & (mean_prop_real != 1), 
    TRUE, FALSE)
  )

log1_bs$glmm_obj[log1_bs$new_fit] <- 
  map(.x = log1_bs$df[log1_bs$new_fit], 
      .f = ~ glmer(
        resp_prop ~ group + (1|id),
        data = .,
        weights = n_trials,
        family = binomial
      )
  )

# sometimes this calculation may fail again
# we can temporarily save the file to explore more
# we succeeded in this case

## LABEL 1
log1_bs$glmm_ch1 <- TRUE
log1_bs$glmm_ch1[!log1_bs$manual_ch1] <- 
  map_lgl(log1_bs$glmm_obj[!log1_bs$manual_ch1], 
          ~ isTRUE(.@optinfo$conv$lme4$code == -1))

## LABEL 2
log1_bs$glmm_ch2 <- TRUE
log1_bs$glmm_ch2[!log1_bs$manual_ch1] <- 
  map_lgl(log1_bs$glmm_obj[!log1_bs$manual_ch1], 
          isSingular)

## Release the workspace
log1_bs$glm_obj <- NULL
log1_bs$glmm_obj <- NULL

########################################################################
## Create type I error rate tibble #####################################
########################################################################

# ----------------- Overall ----------------- 

aov_error <- log1_bs %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- log1_bs %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- log1_bs %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glm_error_Conver <- log1_bs %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver")

glmm_error_Conver <- log1_bs %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_log1_bs_all <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# save the file
save(error_log1_bs_all, file = "./results/error_log1_bs_all.rda")

# plot
p_log1_bs_all <-
  ggplot(error_log1_bs_all,
         aes(x = analysis, y = mean)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.03)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "probability", y = "Type I Error Rate",
       title = "Logistic Simulation bs1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_log1_bs_all

# ----------------- Group by prob ----------------- 

aov_error <- log1_bs %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- log1_bs %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- log1_bs %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glm_error_Conver <- log1_bs %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver")

glmm_error_Conver <- log1_bs %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_log1_bs_prob <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# save the file
save(error_log1_bs_prob, file = "./results/error_log1_bs_prob.rda")

# prob
p_log1_bs_prob <-
  ggplot(error_log1_bs_prob,
         aes(x = prob, y = mean,
             shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.03)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "probability", y = "Type I Error Rate",
       title = "Logistic Simulation bs1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_log1_bs_prob




