## setwd("C:/Users/jy33/OneDrive/Desktop/R/fly_polyandry")

library(tidyverse)
library(ggplot2); theme_set(theme_classic())
library(lme4)
library(DHARMa)
library(car)
library(ggsci)
library(glmmTMB)
library(survival)
library(survminer)
library(coxme)
library(emmeans)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.title.y = element_text(size = 16), 
  axis.text.y = element_text(size = 16))

# Import main fitness data
fitness_data <- read.csv("data/total_offspring_per_fly.csv", stringsAsFactors = TRUE)

fitness_data$treatment <- factor(fitness_data$treatment, # re-order levels for treatment
                                 levels = c("low", "medium", "high")) 

################ MATING RATE ##################
### PLOT CODE ###
ggplot(data = fitness_data, aes(x = treatment, y = mating_rate, fill = treatment)) +
       geom_hline(yintercept = 0.125, linetype = 2, color = "#B2D7ED", linewidth = 1) + 
       geom_hline(yintercept = 0.25, linetype = 2, color = "#078AD6", linewidth = 1) +
       geom_hline(yintercept = 0.5, linetype = 2, color = "#01497c", linewidth = 1) +
       geom_boxplot(alpha = 0.9, outlier.colour = NA) + 
       scale_fill_manual(values=c("#B2D7ED", "#078AD6", "#01497c")) +
       geom_jitter(width = 0.15, alpha = 0.4, height = 0) + 
       ylab("Mating rate") + xlab("Treatment") + ylim(0, 0.5) + 
       theme(legend.position="none") + My_Theme

### MODEL CODE ###
mating_rate_mod <- glmmTMB(data = fitness_data, mating_rate ~ treatment + (1|replicate),
                        dispformula = ~ treatment)

plot(simulateResiduals(mating_rate_mod)) # diagnostic plots; not ideal, explanation below
# this model violates uniformity but it's expected here given the low treatment data 
# it's also better than the alternative models I attempted to fit: 
# removing the dispformula violates homogeneity of var (worse than violating uniformity)
# a binomial model with cbind(matings, days_not_mated) as the response variable leads to a singular fit and violates QQ plot/homogeneity of var/uniformity

Anova(mating_rate_mod)

################## LONGEVITY ###################
curve <- survfit(Surv(longevity, status) ~ treatment, data = fitness_data)

### PLOT CODE ###
ggsurvplot(curve, data = fitness_data, legend = "right", linetype = 1, size = 2, 
           xlim = c(0, 80), color = "strata", 
           palette = (c("#B2D7ED", "#078AD6", "#01497c")),
           font.x = 22, font.y = 22, font.tickslab = 16) + 
           xlab("Female age (days)") + ylab("Proportion of females alive")

### MODEL CODE ###
long_model <- coxme(Surv(longevity, status) ~ treatment + (1|replicate), 
                    data = fitness_data)

Anova(long_model)

################## OFFSPRING ###################
### PLOT CODE ###
# Cumulative offspring plot
ggplot(data = fitness_data, aes(x = treatment, y = total_offspring, fill = treatment)) + 
       geom_boxplot(alpha = 0.9, outlier.colour = NA) + 
       scale_fill_manual(values=c("#B2D7ED", "#078AD6", "#054469")) + 
       scale_color_manual(values=c("#B2D7ED", "#078AD6", "#054469"))+ My_Theme + 
       ylab("Offspring per female") + xlab("Treatment") +
       scale_y_continuous(limits = c(0, 1600), 
                          breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600)) + 
       geom_jitter(width = 0.1, alpha = 0.4) + theme(legend.position="none")

### MODEL CODE ###
offspring_mod <- glmmTMB(data = fitness_data, total_offspring ~ treatment + (1|replicate))

plot(simulateResiduals(offspring_mod)) # # diagnostic plots; looks good

Anova(offspring_mod)




