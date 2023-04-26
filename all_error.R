########################################################################
##                 Type I Error Rate: ANOVA, GLM & GLMM               ##
##                                   ALL                              ##
########################################################################

## ---------------------------------------------------------------------
## Prepare
## ---------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(binom)
library(purrr)
library(lme4)
theme_set(theme_bw(base_size = 14))

## ---------------------------------------------------------------------
## Between-Subject 1 trial pp
## ---------------------------------------------------------------------

load("./logistic simulation/error_log1_bs_all.rda")
load("./diffusion simulation/error_ddm1_bs_all.rda")

error_bs1_1000_all
error_bs1_1000_all$size <- NULL
error_bs1_1000_all$generation <- "logistic"

error_mvt1_bs_all
error_mvt1_bs_all$generation <- "diffusion"

error_1tr_bs <- rbind(error_bs1_1000_all, error_mvt1_bs_all)
error_1tr_bs

error_1tr_bs <- error_1tr_bs %>% 
  mutate(generation = factor(generation,
                             levels = c("logistic", "diffusion"),
                             labels = c("Binomial logistic model",
                                        "Ratcliff diffusion model")),
         analysis = factor(analysis,
                           levels = c("aov", "glm", "glm_noConver",
                                      "glmm", "glmm_noConver"),
                           labels = c("ANOVA", "GLM all", "GLM converged",
                                      "GLMM all", "GLMM converged")
         ))
error_1tr_bs

p_1tr_bs <- 
  ggplot(error_1tr_bs, 
         aes(x = analysis, y = mean, shape = generation)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.5)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Candidate models", y = "Type I Error Rate",
       title = "Between-subject with 1 trial per participant (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_shape_discrete(name = "Data generation model") +
  theme(axis.text = element_text(colour = "black"))

p_1tr_bs

## ---------------------------------------------------------------------
## Between-Subject 100 trials pp
## ---------------------------------------------------------------------

load("~/logistic simulation/error_log100_bs_all.rda")
load("./error_mvt100_bs_all.rda")

error_bs100_1000_all
error_bs100_1000_all$size <- NULL
error_bs100_1000_all$generation <- "logistic"

error_mvt100_bs_all
error_mvt100_bs_all$generation <- "diffusion"

error_100tr_bs <- rbind(error_bs100_1000_all, error_mvt100_bs_all)
error_100tr_bs


error_100tr_bs <- error_100tr_bs %>% 
  mutate(generation = factor(generation,
                             levels = c("logistic", "diffusion"),
                             labels = c("Binomial logistic model",
                                        "Ratcliff diffusion model")))
error_100tr_bs

p_100tr_bs <- 
  ggplot(error_100tr_bs, aes(x = analysis, y = mean, shape = generation)) +
  geom_pointrange(mapping = aes(ymin = lower, ymax = upper), 
                  position = position_dodge(0.5)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Candidate models", y = "Type I Error Rate",
       title = "Between-subject with 100 trials per participant (1000 samples)") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_shape_discrete(name = "Data generation model") +
  scale_x_discrete(breaks = c("aov", "glm", "glmm", "glmm_noConver"),
                   labels = c("ANOVA", "GLM all/converged", 
                              "GLMM all", "GLMM converged")) +
  theme(axis.text = element_text(colour = "black"))

p_100tr_bs




