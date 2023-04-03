########################################################################
##                 Type I Error Rate: ANOVA, GLM & GLMM               ##
##                  Data Simulated from Logistic Model                ##
##                     with Individual Differences                    ##
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

## BS1 ######
load(file = "~/Logistic Simulation/logistic_simulation_id/log_bs_id_30pp1tr_1000runs.rda")
bs1_1000

## -----------------------------------------------------------------------
## BS100
## -----------------------------------------------------------------------

## load data
load(file = "./log_bs_id_30pp100tr_1000runs.rda")
bs100_1000 <- simulated_data

########################################################################
## Set check labels for convergence problems ###########################
########################################################################

## CASE 1: mean_prop_real == 1 #########################################
bs100_1000$manual_ch1 <- ifelse(bs100_1000$mean_prop_real == 1, TRUE, FALSE)

## CASE 2: Inconvergence in GLM ########################################
bs100_1000$glm_ch1 <- map_lgl(bs100_1000$glm_obj, ~ .$converged == FALSE)

## CASE 3: Inconvergence in GLMM #######################################
## LABEL 1
bs100_1000$glmm_ch1 <- map_lgl(bs100_1000$glmm_obj,
                               ~ isTRUE(.@optinfo$conv$lme4$code == -1))
## LABEL 2
bs100_1000$glmm_ch2 <- map_lgl(bs100_1000$glmm_obj, isSingular)

## Release the workspace
mvt100_bs_100$glm_obj = NULL
mvt100_bs_100$glmm_obj = NULL

########################################################################
## Create type I error rate tibble #####################################
########################################################################

# ----------------- Overall ----------------- 

aov_error <- bs100_1000 %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- bs100_1000 %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- bs100_1000 %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glmm_error_Conver <- bs100_1000 %>% 
  filter(glmm_ch1 == FALSE & glmm_ch2 == FALSE) %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_bs100_1000_all <- 
  bind_rows(aov_error, glm_error,
            glmm_error, glmm_error_Conver)

# save the file
save(error_bs100_1000_all, file = "./error_bs100_1000_all.rda")

# plot
p_bs100_1000_all <- 
  ggplot(error_bs100_1000_all,
         aes(x = analysis, y = mean)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.03)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "probability", y = "Type I Error Rate",
       title = "Logistic Simulation bs100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed") 
p_bs100_1000_all

# ----------------- Group by prob ----------------- 

aov_error <- bs100_1000 %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- bs100_1000 %>% 
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

glmm_error_Conver <- bs100_1000 %>% 
  filter(glmm_ch1 == FALSE & glmm_ch2 == FALSE) %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_bs100_1000_prob <- 
  bind_rows(aov_error, glm_error,
            glmm_error, glmm_error_Conver)

# plot
p_bs100_1000_prob <- 
  ggplot(error_bs100_1000_prob,
         aes(x = prob, y = mean, shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.03)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "probability", y = "Type I Error Rate",
       title = "Logistic Simulation bs100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed") 
p_bs100_1000_prob

## -----------------------------------------------------------------------
## BS1
## -----------------------------------------------------------------------

########################################################################
## Set check labels for convergence problems ###########################
########################################################################

## CASE 1: mean_prop_real == 1 #########################################
bs1_1000$manual_ch1 <- ifelse(bs1_1000$mean_prop_real == 1, TRUE, FALSE)

## CASE 2: Inconvergence in GLM ########################################
bs1_1000$glm_ch1 <- TRUE
bs1_1000$glm_ch1[!bs1_1000$manual_ch1] <- 
  map_lgl(bs1_1000$glm_obj[!bs1_1000$manual_ch1], ~ !.$converged)

## CASE 3: Inconvergence in GLMM #######################################
## Recalculate glmm_obj
# new_fit = TRUE needs to recalculate glmm
bs1_1000 <- bs1_1000 %>% 
  mutate(new_fit = if_else(
    # mean_prop_real != 1 but glmm_obj = 1
    map_lgl(glmm_obj, is.numeric) & (mean_prop_real != 1), 
    TRUE, FALSE)
  )

bs1_1000$glmm_obj[bs1_1000$new_fit] <- 
  map(.x = bs1_1000$df[bs1_1000$new_fit], 
      .f = ~ glmer(
        resp_prop ~ group + (1|id),
        data = .,
        weights = n_trials,
        family = binomial
      )
  )

# sometimes this calculation may fail again
# we can temporarily save the file to explore more
tmp <- bs1_1000 %>% 
  filter(new_fit)
save(tmp, file = "./bs1_1000_abnormal_cases.rda",
     compress = "xz")
# We succeeded in this case

## LABEL 1
bs1_1000$glmm_ch1 <- TRUE
bs1_1000$glmm_ch1[!bs1_1000$manual_ch1] <- 
  map_lgl(bs1_1000$glmm_obj[!bs1_1000$manual_ch1], 
          ~ isTRUE(.@optinfo$conv$lme4$code == -1))

## LABEL 2
bs1_1000$glmm_ch2 <- TRUE
bs1_1000$glmm_ch2[!bs1_1000$manual_ch1] <- 
  map_lgl(bs1_1000$glmm_obj[!bs1_1000$manual_ch1], 
          isSingular)

## Release the workspace
bs1_1000$glm_obj <- NULL
bs1_1000$glmm_obj <- NULL

########################################################################
## Create type I error rate tibble #####################################
########################################################################

# ----------------- Overall ----------------- 

aov_error <- bs1_1000 %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- bs1_1000 %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- bs1_1000 %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glm_error_Conver <- bs1_1000 %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver")

glmm_error_Conver <- bs1_1000 %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  group_by(size) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_bs1_1000_all <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# save the file
save(error_bs1_1000_all, file = "./error_bs1_1000_all.rda")

# plot
p_bs1_1000_all <-
  ggplot(error_bs1_1000_all,
         aes(x = analysis, y = mean)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.03)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "probability", y = "Type I Error Rate",
       title = "Logistic Simulation bs1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_bs1_1000_all

# ----------------- Group by prob ----------------- 

aov_error <- bs1_1000 %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- bs1_1000 %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- bs1_1000 %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glm_error_Conver <- bs1_1000 %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver")

glmm_error_Conver <- bs1_1000 %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  group_by(prob) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_bs1_1000_prob <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# prob
p_bs1_1000_prob <-
  ggplot(error_bs1_1000_prob,
         aes(x = prob, y = mean,
             shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.03)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "probability", y = "Type I Error Rate",
       title = "Logistic Simulation bs1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_bs1_1000_prob




