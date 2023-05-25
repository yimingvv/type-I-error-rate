########################################################################
##                      Ratcliff Diffusion Model                      ##
##                         Preprocessing Data                         ##
########################################################################

## Data
## Experiment 1
## Guo, L., Trueblood, J. S., & Diederich, A. (2017). 
## Thinking Fast Increases Framing Effects in Risky Decision Making. 
## Psychological Science, 28(4), 530–543. https://doi.org/10.1177/0956797616689092

library(readxl)
library(dplyr)
library(ggplot2)

## ---------------------------------------------------------------------
## Load Data
## ---------------------------------------------------------------------

## 4 different variants manipulating the color scheme of presentation
## only 3 of those datasets were available to us.
## So, we used these 3 datasets to estimate the parameter

# v1: version1 - original red & green presentation - 44 pp
v1 <- read_excel("UpdatedFiles_Org_RG.xlsx", sheet = "file1.txt")
v1 <- rename(v1, "game" = "Game #")
v1 <- rename(v1, "id" = "Subj #")
v1 <- rename(v1, "response" = "Respose")
v1 <- rename(v1, "rt" = "RT")

# v2: version2 - grey-scale (black & white) presentation - 44 pp
v2 <- read_excel("UpdatedFiles_BW.xlsx", sheet = "file1.txt")
v2 <- rename(v2, "game" = "Game #")
v2 <- rename(v2, "id" = "Subj #")
v2 <- rename(v2, "response" = "Respose")
v2 <- rename(v2, "rt" = "RT")

# v3: version3 - random positioned red & green - 53 pp
v3 <- read_excel("UpdatedFiles_RandRG.xlsx", sheet = "file1.txt")
v3 <- rename(v3, "game" = "Game #")
v3 <- rename(v3, "id" = "Subj #")
v3 <- rename(v3, "response" = "Response")
v3 <- rename(v3, "rt" = "RT")

## ---------------------------------------------------------------------
## Code Data
## ---------------------------------------------------------------------

## Coding Key
#  Game 1-72	 Gain + TP
#  Game 73-144	Loss + TP
#  Game 145-152	Gain + TP + catch
#  Game 153-160	Loss + TP + catch
#  Game 161-232	Gain + noTP
#  Game 233-304	Loss + noTP
#  Game 305-312	Gain + noTP + catch
#  Game 313-320	Loss + noTP + catch
## Response: 1 - Sure (A); 2 - Gamble (B)
## All rt in (sec)	

## create labels for v1
# create labels for frame: gain/loss
v1 <- v1 %>% 
  mutate(frame = ifelse(game %in% c(1:72) |
                          game %in% c(145:152) | 
                          game %in% c(161:232) | 
                          game %in% c(305:312) ,
                        "gain", "loss"))
# create labels for time pressure: tp/ntp
v1 <- v1 %>% 
  mutate(time_pressure = ifelse(game %in% c(1:160),
                                "tp", "ntp"))
# create labels for catch (1) or (0)
v1 <- v1 %>% 
  mutate(catch = ifelse(game %in% c(145:160) |
                          game %in% c(305:320),
                        1, 0))
# create labels for version
v1 <- v1 %>% mutate(version = 1)

## Create labels for v2
# create labels for frame: gain/loss
v2 <- v2 %>% 
  mutate(frame = ifelse(game %in% c(1:72) |
                          game %in% c(145:152) | 
                          game %in% c(161:232) | 
                          game %in% c(305:312) ,
                        "gain", "loss"))
# create labels for time pressure: tp/ntp
v2 <- v2 %>% 
  mutate(time_pressure = ifelse(game %in% c(1:160),
                                "tp", "ntp"))
# create labels for catch (1) or (0)
v2 <- v2 %>% 
  mutate(catch = ifelse(game %in% c(145:160) |
                          game %in% c(305:320),
                        1, 0))
# create labels for version
v2 <- v2 %>% mutate(version = 2)

## Create labels for v3
# create labels for frame: gain/loss
v3 <- v3 %>% 
  mutate(frame = ifelse(game %in% c(1:72) |
                          game %in% c(145:152) | 
                          game %in% c(161:232) | 
                          game %in% c(305:312) ,
                        "gain", "loss"))
# create labels for time pressure: tp/ntp
v3 <- v3 %>% 
  mutate(time_pressure = ifelse(game %in% c(1:160),
                                "tp", "ntp"))
# create labels for catch (1) or (0)
v3 <- v3 %>% 
  mutate(catch = ifelse(game %in% c(145:160) |
                          game %in% c(305:320),
                        1, 0))
# create labels for version
v3 <- v3 %>% mutate(version = 3)

## create unique id for each participant
n_distinct(v1$id) # 49
n_distinct(v2$id) # 49
n_distinct(v3$id) # 53

v2 <- v2 %>% mutate(id = id + 49)
v3 <- v3 %>% mutate(id = id - 300 + 49 + 49 + 1)

## combine all the data sets
allv <- bind_rows(v1, v2, v3)

# create labels for specific catch trial type
# The catch trials had nonequivalent "sure" and "gamble" options
# one of which had a significantly larger expected value.
allv <- allv %>% 
  mutate(catch_type = ifelse(game %in% c(145:152) | # Gain + TP + catch
                                game %in% c(305:312), # Gain + nTP + catch
                              'sure', # gain catch - sure
                              ifelse(game %in% c(153:160) | # Loss + TP + catch
                                       game %in% c(313:320), # Loss + nTP + catch
                                     'gamble', # loss catch - gamble
                                     0)),
         response_gamble = ifelse(response == 2, 1, 0))
allv$response = NULL

## Label
allv <- allv %>% 
  mutate(id = factor(id),
         game = factor(game),
         frame = factor(frame),
         time_pressure = factor(time_pressure),
         catch = factor(catch,
                        levels = c(1, 0),
                        labels = c("catch", "nocatch")),
         version = factor(version),
         catch_type = factor(catch_type))
glimpse(allv)

## ---------------------------------------------------------------------
## Exclude catch trials
## ---------------------------------------------------------------------

allv %>% 
  filter(catch == "catch") %>% 
  group_by(id) %>% 
  summarise(n = n())
# 32 catch trials for each pp

allv_nocatch <- allv %>% 
  filter(catch == "nocatch")

## ---------------------------------------------------------------------
## Remove participants whose props_gamble is too small (< 0.05)
## ---------------------------------------------------------------------

## Participants with very low probability to choose the gamble option

filter_gamble05 <- allv_nocatch %>% 
  group_by(id) %>% 
  summarise(props_gamble = mean(response_gamble)) %>% 
  mutate(props_gamble_05 = ifelse(props_gamble <= 0.05, 1, 0)) %>% 
  filter(props_gamble_05 == 1)
unique(filter_gamble05$id)
# [1] 7   59  61  104 120

allv_nocatch_gamble <- allv_nocatch %>% 
  filter(!id %in% filter_gamble05$id)

## ---------------------------------------------------------------------
## Calculate trials with abnormal rt
## ---------------------------------------------------------------------

## See the rt distribution in two time_pressure conditions
allv_nocatch %>% 
  group_by(time_pressure) %>% 
  summarise(c02 = mean(rt < 0.2),
            c1 = mean(rt > 1),
            c25 = mean(rt > 2.5),
            c3 = mean(rt > 3),
            c35 = mean(rt > 3.5),
            c4 = mean(rt > 4), 
            c6 = mean(rt > 6), 
            c10 = mean(rt > 10))

## Normal range
## tp [0.2, 1]; ntp [0.2, 3.5] (90% included)
## 1000ms or 1s is the required range of the task

## Exclude participants 
## with >= 15% trials outside the chosen time window
allv_rt_proportion <- allv_nocatch_gamble %>% 
  mutate(rt_normal = ifelse(time_pressure == "tp" & 
                              rt >= 0.2 & rt <= 1,
                            1, 
                            ifelse(time_pressure == "ntp" &
                                     rt >= 0.2 & rt <= 3.5, 
                                   1, 0))) %>% 
  group_by(id) %>% 
  summarise(rt_abnormal_prop = 1 - mean(rt_normal)) %>% 
  mutate(rt15ab = ifelse(rt_abnormal_prop >= 0.15, 1, 0))

allv_rt_proportion %>% filter(rt15ab == 1) # 22

filter_rt15ab <- allv_rt_proportion %>% 
  filter(rt15ab == 1)

allv_nocatch_gamble_rt1 <- allv_nocatch_gamble %>% 
  filter(!id %in% filter_rt15ab$id)

allv_nocatch_gamble_rt1 %>% 
  group_by(time_pressure) %>% 
  summarise(c02 = mean(rt < 0.2),
            c1 = mean(rt > 1),
            c35 = mean(rt > 3.5))

## Exclude trials 
## whose rt is outside the chosen time window
data_final <- allv_nocatch_gamble_rt1 %>% 
  mutate(rt_normal = ifelse(time_pressure == "tp" & 
                              rt >= 0.2 & rt <= 1,
                            1, 
                            ifelse(time_pressure == "ntp" &
                                     rt >= 0.2 & rt <= 3.5, 
                                   1, 0))) %>% 
  filter(rt_normal == 1)

## see the final data
glimpse(data_final)

data_final$catch <- NULL
data_final$catch_type <- NULL
data_final$rt_normal <- NULL

glimpse(data_final)

length(unique(data_final$id)) # 124

data_final %>%
  group_by(time_pressure, frame) %>% 
  summarise(n = length(unique(id))) # within-subject design

save(data_final, file = "data_final.rda")

