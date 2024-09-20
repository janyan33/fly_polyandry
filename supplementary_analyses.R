## setwd("C:/Users/jy33/OneDrive/Desktop/R/fly_polyandry")

library(tidyverse) # data wrangling
library(ggplot2); theme_set(theme_classic()) # plotting
library(lme4) # mixed-modelling
library(glmmTMB) # mixed-modelling
library(DHARMa) # diagnostic plots
library(car) # model outputs
library(emmeans) # posthoc analyses
library(ggeffects)

My_Theme = theme(
           axis.title.x = element_text(size = 16),
           axis.text.x = element_text(size = 16),
           axis.title.y = element_text(size = 16), 
           axis.text.y = element_text(size = 16))

############# PILOT GLUED VS. UNGLUED MALES' COURTSHIP BEHAVIOURS ###############
# Results described in methods of the main ms but fig not included bc its minor
glue_data <- read.csv("data/glued_pilot.csv", stringsAsFactors = TRUE)

ggplot(data = glue_data, aes(x = treatment, y = prop_courtship, fill = treatment)) + 
       geom_boxplot() + geom_jitter(width = 0.18, height = 0, alpha = 0.5, size = 2) + 
       My_Theme + theme(legend.position = "none") + ylim(0, 1) +
       ylab("Proportion of trial spent courting") + xlab("Male type") +
       scale_fill_manual(values=c("grey80", "grey80")) 

glue_mod <- glmmTMB(data = glue_data, asin(sqrt(prop_courtship)) ~ treatment + 
                   (1|observer) + (1|day))

plot(simulateResiduals(glue_mod)) # looks good, I arcsin sqrt transformed to satisfy the QQ plot

Anova(glue_mod) 

##################### OFFSPRING PRODUCTION OVER TIME ########################
offspring_batches <- read.csv("data/offspring_over_time.csv", stringsAsFactors = TRUE) %>% 
                     unite("rep_ID", c(replicate, treatment, vial), sep = "_", remove = FALSE) %>% # create ID column
                     filter(rep_ID != "1_low_3" & rep_ID != "1_low_11") %>% # remove rep 1 flies with issues
                     filter(rep_ID != "2_low_4" & rep_ID != "2_medium_16") %>%  # remove rep 2 flies with issues
                     mutate(daily_offspring = total_offspring/days)  # create offspring production rate column

offspring_batches$treatment <- factor(offspring_batches$treatment, 
                                      levels = c("low", "medium", "high")) # re-order factor levels

offspring_batches$female_age <- factor(offspring_batches$female_age, 
                                       levels = c("4-5", "6-7", "8-9", "10-11", "12-13", # re-order factor levels
                                                  "14-15", "16-19", "20-23", "24-27", "28-31",
                                                  "32-35", "36-39", "40-46", "47-53"))

ggplot(data = offspring_batches, aes(x = female_age, y = daily_offspring, fill = treatment)) + geom_boxplot() + 
       scale_fill_manual(values=c("#caf0f8", "#468faf", "#01497c")) + 
       My_Theme + ylab("Daily offspring produced per female") + xlab("Female age (days)")

offspring_batches_block1 <- offspring_batches %>% 
                            filter(female_age == "4-5" | female_age == "6-7" | female_age == "8-9" | 
                                   female_age == "10-11" | female_age == "12-13")

# Zoomed in plot for talk
ggplot(data = offspring_batches_block1, aes(x = female_age, y = daily_offspring, fill = treatment)) + 
       geom_boxplot() + 
       scale_fill_manual(values=c("#caf0f8", "#468faf", "#01497c")) + 
       My_Theme + ylab("Daily offspring produced per female") + xlab("Female age (days)")

## Offspring produced as a function of Low females' last mating opportunity analyses
time_dat <- read.csv("data/offspring_by_batch.csv", stringsAsFactors = TRUE) %>% 
            unite("rep_ID", c(replicate, treatment, vial), sep = "_", remove = FALSE) %>% # create ID column
            filter(rep_ID != "1_low_3" & rep_ID != "1_low_11") %>% # remove rep 1 flies with issues
            filter(rep_ID != "2_low_4" & rep_ID != "2_medium_16") %>%  # remove rep 2 flies with issues
            mutate(daily_offspring = total_offspring/days) %>%  # create offspring production rate column
            filter(treatment != "medium") %>% # Remove medium, keeping high as the ref level
            filter(days_since_mating_opp != "NA")

time_dat$treatment <- factor(time_dat$treatment, 
                             levels = c("low", "high")) # re-order factor levels

time_mod <- glmmTMB(total_offspring ~ treatment*days_since_mating_opp + (1|replicate/rep_ID), 
                    data = time_dat)

Anova(time_mod)

test(emtrends(time_mod, pairwise ~ treatment, var = "days_since_mating_opp")) # get simple effects

# Figure showing sperm/seminal fluid depletion
time_dat$days_since_mating_opp <- as.factor(time_dat$days_since_mating_opp)

ggplot(data = time_dat, aes(x = days_since_mating_opp, y = daily_offspring, fill = treatment)) + geom_boxplot() + 
       scale_fill_manual(values=c("#caf0f8", "#01497c")) + 
       My_Theme + ylab("Daily offspring produced per female") + xlab("Days since last mating opportunity") +
       facet_grid(~ replicate)

########################## COURTSHIP REPLICATE 1 ############################
courtship_1_data <- read.csv("data/r1_courtship.csv") %>% 
                    unite("fem_ID", c(treatment, vial), sep = "_", remove = FALSE) %>% # create ID column
                    filter(fem_ID != "low_3" & fem_ID != "low_11") %>% 
                    filter(mating_latency > 900 | is.na(mating_latency)) %>% # exclude trials where mating latency <15 minutes (89)
                    mutate(bouts_per_hour = courtship_bouts) %>% 
                    mutate(bout_max = ifelse(courtship_bouts == 5 , 1, 0)) %>% 
                    filter(trial_start != "dead")

courtship_1_data$treatment <- factor(courtship_1_data$treatment, # re-order levels for treatment
                                   levels = c("low", "medium", "high")) 

courtship_1_data$fem_ID <- as.factor(courtship_1_data$fem_ID) # turn female ID into factor

######### Female attractiveness analyses
## FIGURE 
ggplot(data = courtship_1_data, aes(x = female_age, y = bout_max)) + 
       geom_smooth(method = "glm", method.args = list(family = binomial), size = 1, 
                   aes(color = treatment, linetype = treatment), fill = "grey80") +
       My_Theme + ylab("Probability of females receiving \n five courtship bouts") + 
       xlab("Female age (days)") + ylim(0, 1) +
       scale_x_continuous(breaks = c(14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34)) + 
       scale_color_manual(values=c("#90D5FF", "#0F70A9", "#01497c"))

## MODEL 
rep_1_courtship_age_mod <- glmmTMB(data = courtship_1_data, 
                              bout_max ~ treatment*female_age + (1|fem_ID),
                              family = binomial())

plot(simulateResiduals(rep_1_courtship_age_mod)) # Mild departure, but overall good
Anova(rep_1_courtship_age_mod)

######### Courtship rate of glued vs. unglued analyses
rep_1_glue_mod <- glmmTMB(data = courtship_1_data, 
                                 bout_max ~ male_type + (1|fem_ID), 
                                 family = binomial())

plot(simulateResiduals(rep_1_glue_mod)) # Looks good
Anova(rep_1_glue_mod)

########################## COURTSHIP REPLICATE 2 ############################
courtship_2_data <- read.csv("data/r2_courtship_trials.csv", stringsAsFactors = TRUE) %>% 
                    unite("ID", c(treatment, fem_ID), sep = "_", remove = FALSE) %>% 
                    filter(female_age != 12)

courtship_2_data$treatment <- factor(courtship_2_data$treatment, # re-order levels for treatment
                                   levels = c("low", "medium", "high")) 

courtship_2_data$ID <- as.factor(courtship_2_data$ID)

### Courtship duration (proportion of trial where females were courted by males)
## FIGURE 
courtship_2_data$female_age <- as.factor(courtship_2_data$female_age)

ggplot(data = courtship_2_data, aes(x = female_age, y = prop_courtship, fill = treatment)) + 
       geom_boxplot() + scale_fill_manual(values=c("#B2D7ED", "#078AD6", "#054469")) + 
       geom_jitter(height = 0, width = 0, alpha = 0.5, size = 3, aes()) + facet_grid(~treatment) +
       My_Theme + ylim(0, 1) + ylab("Proportion of trial courted by male") + xlab("Female age (days)")

## MODEL
r2_courtship_mod <- glmmTMB(data = courtship_2_data, prop_courtship ~ treatment*female_age + 
                            (1|ID/treatment))

plot(simulateResiduals(r2_courtship_mod)) # looks good

Anova(r2_courtship_mod)

# Investigate what's driving significant interaction
em_mod <- emmeans(r2_courtship_mod, ~ female_age*treatment)
test(pairs(em_mod, by = "treatment"))

courtship_low_data <- courtship_2_data %>% 
                      filter(treatment == "low")

r2_courtship_low_mod <- glmmTMB(data = courtship_low_data, prop_courtship ~ female_age + 
                              (1|ID/treatment))

Anova(r2_courtship_low_mod)

