# libraries
library(tidyverse)
library(nlme)
library(car) # for vif()
library(lme4)
library(lmerTest)
library(cowplot)
library(patchwork)
library(ggeffects) # to get predictions
library(sjPlot) # 'plot_model' function
library(performance) 
library(emmeans) # for posthoc test - pairvise comparisons
library(multcomp) # for posthoc tests wit the letters
library(multcompView)
library(MuMIn)
library(r2glmm)
library(DHARMa)
library(ggcorrplot)

# Read data --------------------------------------------------------------------



Phenophase_compos<- read_csv("data/Phenological_composition.csv")


Dat1_1m2 <- read_csv("data/Diversity_1m2.csv") %>% 
  mutate(neophyte_prop = neophyte / SR) %>% 
  left_join(Phenophase_compos,
            by=c("PlotNo", "MowFreq", "MonthLabel", "Subplot" )
            )%>% 
  mutate(MowFreq=fct_relevel(MowFreq,"regular", "reduced", "reduced_sown")) %>% 
  mutate(MonthLabel=fct_relevel(MonthLabel,"March", "May", "July", "September")) 

str(Dat1_1m2)
names(Dat1_1m2)

# Correlation ------------------------------------------------------------------

corl1 <- round(cor(Dat1_1m2 %>% 
                     dplyr::select(cover, biomass, SR, 
                                   evenness, shannon,
                                   phen_Richness, phen_Shannon, phen_evenness) %>% 
                     rename("species richness" =SR,
                            "Shannon diversity" =shannon,
                            "Phenological richness" =phen_Richness,
                            "Phenological Shannon diversity" =phen_Shannon,
                   "Phenological evennessy" =phen_evenness),
                   method = c("pearson"), use = "pairwise.complete.obs"), 2)

corl1

ggcorrplot(corl1,
           hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3,
           colors = c("red", "white", "blue"))



# SR ---------------------------------------------------------------------------
m1_SR <-glmer (SR ~ MonthLabel * MowFreq + (1|PlotNo), data = Dat1_1m2,  
            family = "poisson")

# check multicolinearity
vif(m1_SR)

# check the model
summary(m1_SR)

# check overdispersion
performance::check_overdispersion(m1_SR)
#OR
sum(residuals(m1_SR, type = "pearson")^2) / df.residual(m1_SR)
# No overdispersion detected.

# check predictor effects
Anova(m1_SR)
#or
drop1(m1_SR,  test = "Chi")
# interaction is not significant ---- REMOVE it

m2_SR <-glmer(SR ~ MonthLabel + MowFreq + (1|PlotNo), data = Dat1_1m2,  
            family = "poisson")

check_collinearity(m2_SR)

Anova(m2_SR)
drop1(m2_SR,  test = "Chi")


## R2 (for paper)---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_SR)
# R2m is marginal (for fixed predictors) coefficients of determination.
# R2c is conditional (for fixed and random predictors).
#           R2m       R2c
# delta     0.3577211 0.5536261
#           (this one) (not this)

# Partial R2 for fixed effects
r2glmm::r2beta(m2_SR,  partial = T)


## plots ------------------------------------------------------------------------

# get letters for MowFreq boxplot
emmeans_m2_SR <- cld(emmeans(m2_SR, list(pairwise ~ MowFreq)), 
                  Letters = letters) %>% arrange(MowFreq)
emmeans_m2_SR


# add to the ggplot

ggplot(Dat1_1m2,aes(x=MowFreq,y=SR,col=MowFreq)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MowFreq), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  geom_text(data=emmeans_m2_SR,
            aes(x=MowFreq, y=c(30, 40, 38),
                label=emmeans_m2_SR$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Mowing frequency" , y="Species richness",
       color="Mowing frequency", fill="Mowing frequency")
# MowingFrequency (Treatment is not just mowing but sowing as well so careful) is a significant predictor of Species richness
# Sowing seems to be the main factor, a sown AND regular category would have been nice



#Plot for month
# get letters for MonthLabel boxplot
emmeans_m2_SR2 <- cld(emmeans(m2_SR, list(pairwise ~ MonthLabel)), 
                     Letters = letters) %>% arrange(MonthLabel)
emmeans_m2_SR2

ggplot(Dat1_1m2,aes(x=MonthLabel,y=SR,col=MonthLabel)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MonthLabel), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_SR2,
            aes(x=MonthLabel, y=c(38, 38, 35, 35),
                label=emmeans_m2_SR2$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Month of sampling" , y="Species richness") 


# evenness  ------------------------------------------------------------------------

m1_ev <-lmer (log(evenness) ~ MonthLabel * MowFreq + (1|PlotNo), data = Dat1_1m2)

# assumptions
plot(m1_ev)
qqnorm(resid(m1_ev))
qqline(resid(m1_ev))


# check model
summary(m1_ev)
car::Anova(m1_ev)
# or
drop1(m1_ev, test="F")
# no interaction

m2_ev <-lmer (log(evenness) ~ MonthLabel + MowFreq + (1|PlotNo), data = Dat1_1m2)
Anova(m2_ev)


## Task 2 ----------------------------------------------------------------------
#Plot for MowFreq
emmeans_m2_ev <- cld(emmeans(m2_ev, list(pairwise ~ MowFreq)), 
                     Letters = letters) %>% arrange(MowFreq)
emmeans_m2_ev


# add to the ggplot

ggplot(Dat1_1m2,aes(x=MowFreq,y=evenness,col=MowFreq)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MowFreq), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_ev,
            aes(x=MowFreq, y=c(13, 13, 15),
                label=emmeans_m2_ev$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Mowing frequency" , y="Species evenness",
       color="Mowing frequency", fill="Mowing frequency")

#Plot for month
# get letters for MonthLabel boxplot
emmeans_m2_ev2 <- cld(emmeans(m2_ev, list(pairwise ~ MonthLabel)), 
                      Letters = letters) %>% arrange(MonthLabel)
emmeans_m2_ev2

ggplot(Dat1_1m2,aes(x=MonthLabel,y=evenness,col=MonthLabel)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MonthLabel), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_ev2,
            aes(x=MonthLabel, y=c(15, 15, 15, 15),
                label=emmeans_m2_ev2$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Month of sampling" , y="Species evenness") 


# biomass  ------------------------------------------------------------------------


m1_mass <-lmer (log(biomass)   ~ MonthLabel * MowFreq + (1|PlotNo), data = Dat1_1m2 )

# assumptions
plot(m1_mass)
qqnorm(resid(m1_mass))
qqline(resid(m1_mass))


# check model
summary(m1_mass)
car::Anova(m1_mass)
# or
drop1(m1_mass, test="F")

#Plot for interaction ---------------------
emmeans_m1_mass <- cld(emmeans(m1_mass, list(pairwise ~ MowFreq:MonthLabel)), 
                     Letters = letters) %>% arrange(MowFreq, MonthLabel)
emmeans_m1_mass

ggplot(Dat1_1m2,aes(x=MowFreq,y=biomass,col=MonthLabel)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MonthLabel), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.2, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m1_mass,
            aes(x=MowFreq, y=c(4000, 4000, 4000, 4000, 4000, 
                               4000, 4000, 4000, 4000, 4000, 
                               4000, 4000),
                label=emmeans_m1_mass$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Mowing frequency" , y="total biomass (cover*height)",
       color="Mowing frequency", fill="Mowing frequency")

# neophyte_prop ----------------------------------------------------------------

m1_neo <-glmer (neophyte_prop ~ MonthLabel + MowFreq + (1|PlotNo), data = Dat1_1m2,  
            family = "binomial", weights = SR)


performance::check_convergence(m1_neo)

# check the model
summary(m1_neo)


# check overdispersion
performance::check_overdispersion(m1_neo)
#'* not the same as the function below, maybe because the DHARMa package is required? ***********
#OR
sum(residuals(m1_neo, type = "pearson")^2) / df.residual(m1_neo)
# No overdispersion detected.

# check predictor effects
Anova(m1_neo)
#or
drop1(m1_neo,  test = "Chi")

## Task 4 ----------------------------------------------------------------------
#Plot for MowFreq
emmeans_m1_neo <- cld(emmeans(m1_neo, list(pairwise ~ MowFreq)), 
                       Letters = letters) %>% arrange(MowFreq)
emmeans_m1_neo

ggplot(Dat1_1m2,aes(x=MowFreq,y=neophyte_prop,col=MowFreq)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MowFreq), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m1_neo,
            aes(x=MowFreq, y=c(0.2, 0.2, 0.2),
                label=emmeans_m1_neo$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Mowing frequency" , y="Neophyte propability in %",
       color="Mowing frequency", fill="Mowing frequency")

#Plot for month
# get letters for MonthLabel boxplot
emmeans_m1_neo2 <- cld(emmeans(m1_neo, list(pairwise ~ MonthLabel)), 
                        Letters = letters) %>% arrange(MonthLabel)
emmeans_m1_neo2

ggplot(Dat1_1m2,aes(x=MonthLabel,y=neophyte_prop,col=MonthLabel)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MonthLabel), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m1_neo2,
            aes(x=MonthLabel, y=c(0.2, 0.2, 0.2, 0.2),
                label=emmeans_m1_neo2$.group),vjust=0.5, hjust=0.6, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Month of sampling" , y="Neophyte probability in %") 




# phen_Richness ---------------------------------------------------------------------------
m1_phen_Richness <-glmer (phen_Richness ~ MonthLabel * MowFreq + (1|PlotNo), data = Dat1_1m2,  
                          family = "poisson")

# check multicolinearity
vif(m1_phen_Richness)

# check the model
summary(m1_phen_Richness)

# check overdispersion
performance::check_overdispersion(m1_phen_Richness)
#OR
sum(residuals(m1_phen_Richness, type = "pearson")^2) / df.residual(m1_phen_Richness)
# No overdispersion detected.

# check predictor effects
Anova(m1_phen_Richness)
#or
drop1(m1_phen_Richness,  test = "Chi")
# interaction is not significant ---- REMOVE it

m2_phen_Richness <-glmer(phen_Richness ~ MonthLabel + MowFreq + (1|PlotNo), data = Dat1_1m2,  
                         family = "poisson")

check_collinearity(m2_phen_Richness)

Anova(m2_phen_Richness)
drop1(m2_phen_Richness,  test = "Chi")


## R2 (for paper)---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_phen_Richness)
# R2m is marginal (for fixed predictors) coefficients of determination.
# R2c is conditional (for fixed and random predictors).
#           R2m       R2c
# delta     0.3577211 0.5536261
#           (this one) (not this)

# Partial R2 for fixed effects
r2glmm::r2beta(m2_phen_Richness,  partial = T)


## plots ------------------------------------------------------------------------

# get letters for MowFreq boxplot
emmeans_m2_phen_Richness <- cld(emmeans(m2_phen_Richness, list(pairwise ~ MowFreq)), 
                                Letters = letters) %>% arrange(MowFreq)
emmeans_m2_phen_Richness


# add to the ggplot

ggplot(Dat1_1m2,aes(x=MowFreq,y=phen_Richness,col=MowFreq)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MowFreq), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  geom_text(data=emmeans_m2_phen_Richness,
            aes(x=MowFreq, y=c(7, 7, 7),
                label=emmeans_m2_phen_Richness$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Mowing frequency" , y="Phenological richness",
       color="Mowing frequency", fill="Mowing frequency")
# MowingFrequency (Treatment is not just mowing but sowing as well so careful) is a significant predictor of Species richness
# Sowing seems to be the main factor, a sown AND regular category would have been nice



#Plot for month
# get letters for MonthLabel boxplot
emmeans_m2_phen_Richness2 <- cld(emmeans(m2_phen_Richness, list(pairwise ~ MonthLabel)), 
                                 Letters = letters) %>% arrange(MonthLabel)
emmeans_m2_phen_Richness2

ggplot(Dat1_1m2,aes(x=MonthLabel,y=phen_Richness,col=MonthLabel)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MonthLabel), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_phen_Richness2,
            aes(x=MonthLabel, y=c(5, 5, 5, 5),
                label=emmeans_m2_phen_Richness2$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Month of sampling" , y="Phenological richness") 


# phen_evenness  ------------------------------------------------------------------------

m1_ev <-lmer (log(phen_evenness) ~ MonthLabel * MowFreq + (1|PlotNo), data = Dat1_1m2)

# assumptions
plot(m1_ev)
qqnorm(resid(m1_ev))
qqline(resid(m1_ev))


# check model
summary(m1_ev)
car::Anova(m1_ev)
# or
drop1(m1_ev, test="F")



#Plot for interaction
emmeans_m1_ev <- cld(emmeans(m1_ev, list(pairwise ~ MowFreq:MonthLabel)), 
                     Letters = letters) %>% arrange(MowFreq, MonthLabel)
emmeans_m1_ev


# add to the ggplot

ggplot(Dat1_1m2,aes(x=MowFreq,y=phen_evenness,col=MonthLabel)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MonthLabel), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.2, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "top") +
  geom_text(data=emmeans_m1_ev,
            aes(x=MowFreq, y=c(4, 4, 4,
                               4, 4, 4,
                               4, 4, 4,
                               4, 4, 4),
                label=emmeans_m1_ev$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Mowing frequency" , y="Phenological evenness",
       color="Mowing frequency", fill="Mowing frequency")

