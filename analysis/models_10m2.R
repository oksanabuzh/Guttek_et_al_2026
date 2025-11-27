# libraries
library(tidyverse)
library(nlme)
library(car) # for vif()
library(lme4)
library(lmerTest)
#library(cowplot)
# library(patchwork)
library(ggeffects) # to get predictions
library(sjPlot) # 'plot_model' function
library(performance) 

library(emmeans) # for posthoc test - pairvise comparisons
library(multcomp) # for posthoc tests wit the letters


# Read data

Dat1_1m2 <- read_csv("data/Diversity_10m2.csv") %>% 
  mutate(MowFreq=fct_relevel(MowFreq,"regular", "reduced", "reduced_sown")) %>% 
  mutate(MonthLabel=fct_relevel(MonthLabel,"March", "May", "July")) %>% 
  mutate(neophyte_prop = neophyte / SR)
  



str(Dat1_1m2)
names(Dat1_1m2)


# Correlation -----------------------------------------------------------------

corl1 <- round(cor(Dat1_1m2 %>% 
                     dplyr::select(cover, biomass, SR, evenness, shannon) %>% 
                     rename("species richness" =SR,
                            "Shannon diversity" =shannon ),
                   method = c("pearson"), use = "pairwise.complete.obs"), 2)

corl1

library(ggcorrplot)
ggcorrplot(corl1,
           hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3,
           colors = c("red", "white", "blue"))



# SR ------------------------------------------------------------------------
m1 <-glm (SR ~ MonthLabel * MowFreq , data = Dat1_1m2,  
            family = "poisson")


# check the model
summary(m1)


# check overdispersion
performance::check_overdispersion(m1)
#OR
sum(residuals(m1, type = "pearson")^2) / df.residual(m1)
# No overdispersion detected.

# check predictor effects
Anova(m1)
#or
drop1(m1,  test = "Chi")

# interaction is not significant ---- REMOVE it

m2 <-glm(SR ~ MonthLabel + MowFreq, data = Dat1_1m2,  
            family = "poisson")

performance::check_overdispersion(m2)

Anova(m2)
drop1(m2,  test = "Chi")

