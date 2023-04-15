########################################################################
##                 Type I Error Rate: ANOVA, GLM & GLMM               ##
##               Data Simulated from Drift Diffusion Model            ##
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
theme_set(theme_bw(base_size = 12))

## ---------------------------------------------------------------------
## MVT100: Between-Subject
## ---------------------------------------------------------------------

## load data
load(file = "./multi_ddm_bs_30pp100tr_1000runs.rda")
mvt100_bs_1000 <- sim_mvt

########################################################################
## Set check labels for convergence problems ###########################
########################################################################

## CASE 1: mean_prop_real == 1 #########################################
mvt100_bs_100$manual_ch1 <- 
  ifelse(mvt100_bs_100$mean_prop_real == 1, TRUE, FALSE)

## CASE 2: Inconvergence in GLM ########################################
mvt100_bs_100$glm_ch1 <- TRUE
mvt100_bs_100$glm_ch1[!mvt100_bs_100$manual_ch1] <- 
  map_lgl(mvt100_bs_100$glm_obj[!mvt100_bs_100$manual_ch1],
          ~ .$converged == FALSE)

## CASE 3: Inconvergence in GLMM #######################################
# LABEL 1 
mvt100_bs_100$glmm_ch1 <- TRUE
mvt100_bs_100$glmm_ch1[!mvt100_bs_100$manual_ch1] <- 
  map_lgl(mvt100_bs_100$glmm_obj[!mvt100_bs_100$manual_ch1],
          ~ isTRUE(.@optinfo$conv$lme4$code == -1))
# LABEL 2
mvt100_bs_100$glmm_ch2 <- TRUE
mvt100_bs_100$glmm_ch2[!mvt100_bs_100$manual_ch1] <- 
  map_lgl(mvt100_bs_100$glmm_obj[!mvt100_bs_100$manual_ch1],
          isSingular)

## Release the workspace
mvt100_bs_100$glm_obj = NULL
mvt100_bs_100$glmm_obj = NULL

########################################################################
## Create type I error rate tibble #####################################
########################################################################

# ----------------- Overall ----------------- 
aov_error <- mvt100_bs_1000 %>%
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- mvt100_bs_1000 %>%
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- mvt100_bs_1000 %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

# glm_error_Conver: no such case

glmm_error_Conver <- mvt100_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_mvt100_bs_all <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glmm_error_Conver)

# save the file
save(error_mvt100_bs_all, file = "./error_mvt100_bs_all.rda")

# Plot
p_mvt100_bs_all <- 
  ggplot(error_mvt100_bs_all, aes(x = analysis, y = mean)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper)
                  #                  , position = position_dodge(0.5)
  ) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "Analysis", y = "Type I Error Rate",
       title = "Diffusion Simulation mvt100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt100_bs_all

# ----------------- Group by v ----------------- 

aov_error <- mvt100_bs_1000 %>%
  group_by(v_mean) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- mvt100_bs_1000 %>%
  group_by(v_mean) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- mvt100_bs_1000 %>% 
  group_by(v_mean) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glmm_error_Conver <- mvt100_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  group_by(v_mean) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_mvt100_bs_v <- 
  bind_rows(aov_error, glm_error, glmm_error, glmm_error_Conver)
error_mvt100_bs_v

# save the file
save(error_mvt1_bs_all, file = "./error_mvt1_bs_all.rda")

# Plot
p_mvt100_bs_v <- 
  ggplot(error_mvt100_bs_v, aes(x = v_mean, y = mean,
                              shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.5)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = 0:8, minor_breaks = NULL) + 
  labs(x = "Location of v", y = "Type I Error Rate",
       title = "Diffusion Simulation mvt100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt100_bs_v

# ----------------- Group by sv ----------------- 

aov_error <- mvt100_bs_1000 %>%
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- mvt100_bs_1000 %>%
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- mvt100_bs_1000 %>% 
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glmm_error_Conver <- mvt100_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_mvt100_bs_sv <- 
  bind_rows(aov_error, glm_error, glmm_error, 
            glmm_error_Conver)

error_mvt100_bs_sv

# Plot
p_mvt100_bs_sv <- 
  ggplot(error_mvt100_bs_sv, aes(x = sv_mean, y = mean,
                               shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.3)) +
  coord_cartesian(ylim = c(0, 1)) + 
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2,
                                2.5, 3, 3.5, 4),
                     minor_breaks = NULL) +
  labs(x = "Location of sv", y = "Type I Error Rate",
       title = "Diffusion Simulation mvt100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt100_bs_sv

# ----------------- Group by mean_prop_real ----------------- 

aov_error <- mvt100_bs_1000 %>%
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

glm_error <- mvt100_bs_1000 %>%
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

glmm_error <- mvt100_bs_1000 %>% 
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

glmm_error_Conver <- mvt100_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

error_mvt100_bs_90prop <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glmm_error_Conver)

error_mvt100_bs_90prop

# Plot
p_mvt100_bs_90prop <- 
  ggplot(error_mvt100_bs_90prop, aes(x = more_90, y = mean,
                                   shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.3)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "Response proportion", y = "Type I Error Rate",
       title = "Diffusion Simulation mvt100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt100_bs_90prop

# ----------------- Group by differences in resp_prop ----------------- 

# multiple ifelse to transfer diff_props into several integers
range_diff <- function(diff_props){
  round_diff <- ifelse(
    abs(diff_props) >= 0 & abs(diff_props) < 0.02, 0.02,
    ifelse(abs(diff_props) >= 0.02 & abs(diff_props) < 0.04, 0.04,
           ifelse(abs(diff_props) >= 0.04 & abs(diff_props) < 0.06, 0.06,
                  ifelse(abs(diff_props) >= 0.06 & abs(diff_props) < 0.08, 0.08,
                         ifelse(abs(diff_props) >= 0.08 & abs(diff_props) < 0.1, 0.1,
                                ifelse(abs(diff_props) >= 0.1 & abs(diff_props) < 0.12, 0.12,
                                       ifelse(abs(diff_props) >= 0.12 & abs(diff_props) < 0.14, 0.14,
                                              ifelse(abs(diff_props) >= 0.14 & abs(diff_props) < 0.16, 0.16,
                                                     ifelse(abs(diff_props) >= 0.16 & abs(diff_props) < 0.18, 0.18,
                                                            ifelse(abs(diff_props) >= 0.18, 0.2, NA))))))))))
}

aov_error <- mvt100_bs_1000 %>%
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  group_by(range_diff_props) %>% ###### new group_by() here or error (?)
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- mvt100_bs_1000 %>%
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  # group_by(range_diff_props) %>% ##### here no need ?
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- mvt100_bs_1000 %>% 
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  # group_by(range_diff_props) %>% ###### here no need ?
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glmm_error_Conver <- mvt100_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  group_by(range_diff_props) %>% ###### new group_by() here or error (?)
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_mvt100_bs_range_diff <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glmm_error_Conver)

error_mvt100_bs_range_diff

# Plot
p_mvt100_bs_diff <- 
  ggplot(error_mvt100_bs_range_diff,
         aes(x = range_diff_props, y = mean,
             shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.01)) +
  coord_cartesian(ylim = c(0, 1)) + 
  scale_x_continuous(breaks = c(0.02, 0.04, 0.06, 0.08, 0.1, 
                                0.12, 0.14, 0.16, 0.18, 0.2),
                     minor_breaks = NULL) + 
  labs(x = "Range of differences in response proportion",
       y = "Type I Error Rate",
       title = "Diffusion Simulation mvt100 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt100_bs_diff

## ---------------------------------------------------------------------
## MVT1: Between-Subject
## ---------------------------------------------------------------------

## load data
load(file = "./multi_ddm_bs_30pp1tr_1000runs.rda")
mvt1_bs_1000 <- sim_mvt

########################################################################
## Set check labels for convergence problems ###########################
########################################################################

## CASE 1: mean_prop_real == 1 #########################################
mvt1_bs_1000$manual_ch1 <- 
  ifelse(mvt1_bs_1000$mean_prop_real == 1, TRUE, FALSE)

## CASE 2: Inconvergence in GLM ########################################
mvt1_bs_1000$glm_ch1 <- TRUE
mvt1_bs_1000$glm_ch1[!mvt1_bs_1000$manual_ch1] <- 
  map_lgl(mvt1_bs_1000$glm_obj[!mvt1_bs_1000$manual_ch1], ~ !.$converged)

## CASE 3: Inconvergence in GLMM #######################################
## Recalculate glmm_obj
mvt1_bs_1000 <- mvt1_bs_1000 %>% 
  mutate(new_fit = if_else(
    # mean_prop_real != 1 but glmm_obj = 1
    map_lgl(glmm_obj, is.numeric) & (mean_prop_real != 1), 
    TRUE, FALSE)
  )

sum(mvt1_bs_1000$new_fit)
     
mvt1_bs_1000$glmm_obj[mvt1_bs_1000$new_fit] <- 
  map(.x = mvt1_bs_1000$df[mvt1_bs_1000$new_fit], 
      .f = ~ glmer(
        resp_prop ~ group + (1|id),
        data = .,
        weights = n_trials,
        family = binomial
      )
  )

# sometimes this calculation may fail again
# we can temporarily save the file to explore more
tmp <- mvt1_bs_1000 %>% 
  filter(new_fit)
save(tmp, file = "./mvt1_1000_abnormal_cases.rda",
     compress = "xz")

# We failed to recalculate, so we set a manual_ch2 to mark these cases
mvt1_bs_1000 <- mvt1_bs_1000 %>% 
  mutate(manual_ch2 = if_else(
    map_lgl(glmm_obj, is.numeric) & (mean_prop_real != 1), 
    TRUE, FALSE)
  )

# combine two manual_ch into one for glmm
mvt1_bs_1000 <- mvt1_bs_1000 %>% 
  mutate(manual_ch = if_else(
    (manual_ch1 == TRUE) | (manual_ch2 == TRUE), 
    TRUE, FALSE)
  )

## LABEL 1
mvt1_bs_1000$glmm_ch1 <- TRUE
mvt1_bs_1000$glmm_ch1[!mvt1_bs_1000$manual_ch] <- 
  map_lgl(mvt1_bs_1000$glmm_obj[!mvt1_bs_1000$manual_ch], 
          ~ isTRUE(.@optinfo$conv$lme4$code == -1))
## LABEL 2
mvt1_bs_1000$glmm_ch2 <- TRUE
mvt1_bs_1000$glmm_ch2[!mvt1_bs_1000$manual_ch] <- 
  map_lgl(mvt1_bs_1000$glmm_obj[!mvt1_bs_1000$manual_ch],
          isSingular)

## Release the workspace
mvt1_bs_1000$glm_obj <- NULL
mvt1_bs_1000$glmm_obj <- NULL

########################################################################
## Create type I error rate tibble #####################################
########################################################################

# ----------------- Overall ----------------- 
aov_error <- mvt1_bs_1000 %>%
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- mvt1_bs_1000 %>%
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- mvt1_bs_1000 %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver")

glmm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_mvt1_bs_all <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# Plot
p_mvt1_bs_all <- 
  ggplot(error_mvt1_bs_all, aes(x = analysis, y = mean)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "Analysis", y = "Type I Error Rate",
       title = "Diffusion Simulation mvt1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt1_bs_all

# ----------------- Group by v ----------------- 

aov_error <- mvt1_bs_1000 %>%
  group_by(v_mean) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- mvt1_bs_1000 %>%
  group_by(v_mean) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- mvt1_bs_1000 %>% 
  group_by(v_mean) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  group_by(v_mean) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver")

glmm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  group_by(v_mean) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_mvt1_bs_v <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# Plot
p_mvt1_bs_v <- 
  ggplot(error_mvt1_bs_v, aes(x = v_mean, y = mean,
                              shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.5)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = 0:8, minor_breaks = NULL) + 
  labs(x = "Location of v", y = "Type I Error Rate",
       title = "Diffusion Simulation mvt1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt1_bs_v

# ----------------- Group by sv ----------------- 

aov_error <- mvt1_bs_1000 %>%
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- mvt1_bs_1000 %>%
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- mvt1_bs_1000 %>% 
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver")

glmm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  group_by(sv_mean) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_mvt1_bs_sv <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# Plot
p_mvt1_bs_sv <- 
  ggplot(error_mvt1_bs_sv, aes(x = sv_mean, y = mean,
                              shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.3)) +
  coord_cartesian(ylim = c(0, 1)) + 
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2,
                                2.5, 3, 3.5, 4),
                     minor_breaks = NULL) +
  labs(x = "Location of sv", y = "Type I Error Rate",
       title = "Diffusion Simulation mvt1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt1_bs_sv

# ----------------- Group by mean_prop_real ----------------- 

aov_error <- mvt1_bs_1000 %>%
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

glm_error <- mvt1_bs_1000 %>%
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

glmm_error <- mvt1_bs_1000 %>% 
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

glm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

glmm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  mutate(more_90 = ifelse(mean_prop_real >= 0.9, 1, 0)) %>% 
  group_by(more_90) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver",
         more_90 = factor(more_90, levels = c(0, 1),
                          labels = c("Less than 90", "More than 90")))

error_mvt1_bs_90prop <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# Plot
p_mvt1_bs_90prop <- 
  ggplot(error_mvt1_bs_90prop, aes(x = more_90, y = mean,
                               shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.3)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "Response proportion", y = "Type I Error Rate",
       title = "Diffusion Simulation mvt1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt1_bs_90prop

# ----------------- Group by differences in resp_prop ----------------- 

# multiple ifelse to transfer diff_props into several integers
range_diff <- function(diff_props){
  round_diff <- ifelse(
    abs(diff_props) >= 0 & abs(diff_props) < 0.02, 0.02,
    ifelse(abs(diff_props) >= 0.02 & abs(diff_props) < 0.04, 0.04,
           ifelse(abs(diff_props) >= 0.04 & abs(diff_props) < 0.06, 0.06,
                  ifelse(abs(diff_props) >= 0.06 & abs(diff_props) < 0.08, 0.08,
                         ifelse(abs(diff_props) >= 0.08 & abs(diff_props) < 0.1, 0.1,
                                ifelse(abs(diff_props) >= 0.1 & abs(diff_props) < 0.12, 0.12,
                                       ifelse(abs(diff_props) >= 0.12 & abs(diff_props) < 0.14, 0.14,
                                              ifelse(abs(diff_props) >= 0.14 & abs(diff_props) < 0.16, 0.16,
                                                     ifelse(abs(diff_props) >= 0.16 & abs(diff_props) < 0.18, 0.18,
                                                            ifelse(abs(diff_props) >= 0.18, 0.2, NA))))))))))
}

aov_error <- mvt1_bs_1000 %>%
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(aov_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "aov")

glm_error <- mvt1_bs_1000 %>%
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm")

glmm_error <- mvt1_bs_1000 %>% 
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm")

glm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch1 != TRUE & glm_ch1 != TRUE) %>% 
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(glm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glm_Conver")

glmm_error_Conver <- mvt1_bs_1000 %>% 
  filter(manual_ch != TRUE & glmm_ch1 != TRUE & glmm_ch2 != TRUE) %>% 
  mutate(range_diff_props = range_diff(diff_props = diff_props)) %>% 
  group_by(range_diff_props) %>% 
  summarise(below_05 = sum(glmm_p < 0.05),
            n = n()) %>% 
  mutate(binom.logit(x = below_05, n = n)) %>% 
  mutate(analysis = "glmm_Conver")

error_mvt1_bs_range_diff <- 
  bind_rows(aov_error, glm_error, glmm_error,
            glm_error_Conver, glmm_error_Conver)

# Plot
p_mvt1_bs_diff <- 
  ggplot(error_mvt1_bs_range_diff,
         aes(x = range_diff_props, y = mean,
             shape = analysis)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.01)) +
  coord_cartesian(ylim = c(0, 1)) + 
  scale_x_continuous(breaks = c(0.02, 0.04, 0.06, 0.08, 0.1, 
                                0.12, 0.14, 0.16, 0.18, 0.2),
                     minor_breaks = NULL) + 
  labs(x = "Range of differences in response proportion",
       y = "Type I Error Rate",
       title = "Diffusion Simulation mvt1 (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
p_mvt1_bs_diff

## ---------------------------------------------------------------------
## MVT: Repeated-Measure
## ---------------------------------------------------------------------








