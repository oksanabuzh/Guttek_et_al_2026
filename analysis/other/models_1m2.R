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

Mowing_data <- read_csv("data/raw_data/mowing_events_2025_DB.csv") %>% 
  dplyr::select(-mowing_events_2025) %>% 
  pivot_longer(cols = c(September,	July,	May,	March),
               names_to = "Month",
               values_to = "n_mow_events_befre_sampling") %>% 
  relocate(Month, .after=Subplot)


Funct_div <- read_csv("data/processed_data/Functional_Diversity_1m2.csv") 


Dat1_1m2 <- read_csv("data/processed_data/Diversity_phenology_1m2.csv") %>% 
  left_join(Funct_div, by=c("PlotNo", "Subplot", "Month")) %>% 
  left_join(Mowing_data, by=c("PlotNo", "Subplot", "Month")) %>%
  relocate(c(MowFreq, n_mow_events_befre_sampling), .after=Subplot) %>%
  mutate(MowFreq=ifelse(MowFreq == "reduced_sown", "reduced & sowing", MowFreq)) %>% 
  mutate(MowFreq=fct_relevel(MowFreq,"regular", 
                             "reduced", 
                             "reduced & sowing")) %>% 
  mutate(Month=fct_relevel(Month,"March", "May", "July", "September")) %>% 
  mutate(PlotNo=factor(PlotNo))

str(Dat1_1m2)
names(Dat1_1m2)

# assigning colors for the plots
MowFreq_col <- c("#F8766D", "#00B0F6","#00BA38")
Month_col <-  c( "orange", "#287271", "#6D326D","brown") 
  



# Exploratory plots ------------------------------------------------------------------

# Exploratory statistics -------------------------------------------------------

# With Carolin's data:

Dat1_1m2 %>%
  dplyr::select(SR, shrub_cov_buff20, tree_cov_buff20, tree_cover_pct, green_cover_pct,
                blue_cover_pct, protected_biotopes_polygons_cover_pct,
                protected_biotopes_polygons_count, protected_biotopes_lines_count,
                protected_biotopes_points_count, protected_biotopes_count_sum,
                road_density_km_per_ha, parking_area_pct, impervious_pct,
                ulc_pct, poi_count, building_density, building_footprint_cover) %>%
  pivot_longer(-SpecRichness, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = SpecRichness)) +
  geom_point( ) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  labs(x = NULL, y = "Species Richness")




# SR ---------------------------------------------------------------------------
m1_SR <-glmer (SR ~ Month * MowFreq +  n_mow_events_befre_sampling + 
                 (1|PlotNo), data = Dat1_1m2,  
            family = "poisson")

# check overdispersion
check_overdispersion(m1_SR)
# No overdispersion detected.

# check interaction
drop1(m1_SR,  test = "Chi")
# interaction is not significant, remove the interaction term

m2_SR <-glmer(SR ~ Month + MowFreq + n_mow_events_befre_sampling + 
                (1|PlotNo), data = Dat1_1m2,  
              family = "poisson")


# check multicolinearity
check_collinearity(m2_SR)

# check the model
summary(m2_SR)

# check overdispersion
check_overdispersion(m2_SR)
# No overdispersion detected.

# check predictor effects
Anova(m2_SR)
#or
drop1(m2_SR,  test = "Chi")

## R2 (for paper)---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_SR)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_SR,  partial = T)


## plots ------------------------------------------------------------------------

### MowFreq ------

max.SR_MowFreq <-  Dat1_1m2 %>% 
  summarise(max_SR = max(SR), .by = MowFreq)

emmeans_m2_SR <- cld(emmeans(m2_SR, list(pairwise ~ MowFreq)), 
                  Letters = letters) %>% arrange(MowFreq)
emmeans_m2_SR


ggplot(Dat1_1m2,aes(x=MowFreq,y=SR,col=MowFreq)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=MowFreq), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  geom_text(data=emmeans_m2_SR,
            aes(x=MowFreq, y=max.SR_MowFreq$max_SR+4,
                label=emmeans_m2_SR$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  scale_color_manual(values = MowFreq_col) +
  scale_fill_manual(values = MowFreq_col) +
  labs(x="Management" , y="Species richness",
       color="Management", fill="Management")



# Month ------

max.SR_Month <-  Dat1_1m2 %>% 
  summarise(max_SR = max(SR), .by = Month)

emmeans_m2_SR2 <- cld(emmeans(m2_SR, list(pairwise ~ Month)), 
                     Letters = letters) %>% arrange(Month)
emmeans_m2_SR2

ggplot(Dat1_1m2,aes(x=Month,y=SR,col=Month)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=Month), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  geom_text(data=emmeans_m2_SR2,
            aes(x=Month, y=max.SR_Month$max_SR+4,
                label=emmeans_m2_SR2$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Month of sampling" , y="Species richness") 


# evenness  ------------------------------------------------------------------------

m1_ev <-lmer (log(evenness) ~ Month * MowFreq + (1|PlotNo), data = Dat1_1m2)

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

m2_ev <-lmer (log(evenness) ~ Month + MowFreq + (1|PlotNo), data = Dat1_1m2)
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
# get letters for Month boxplot
emmeans_m2_ev2 <- cld(emmeans(m2_ev, list(pairwise ~ Month)), 
                      Letters = letters) %>% arrange(Month)
emmeans_m2_ev2

ggplot(Dat1_1m2,aes(x=Month,y=evenness,col=Month)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=Month), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_ev2,
            aes(x=Month, y=c(15, 15, 15, 15),
                label=emmeans_m2_ev2$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Month of sampling" , y="Species evenness") 


# biomass  ------------------------------------------------------------------------


m1_mass <-lmer (log(biomass)   ~ Month * MowFreq + (1|PlotNo), data = Dat1_1m2 )

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
emmeans_m1_mass <- cld(emmeans(m1_mass, list(pairwise ~ MowFreq:Month)), 
                     Letters = letters) %>% arrange(MowFreq, Month)
emmeans_m1_mass

ggplot(Dat1_1m2,aes(x=MowFreq,y=biomass,col=Month)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=Month), pch=21,
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

m1_neo <-glmer (neophyte_prop ~ Month + MowFreq + (1|PlotNo), data = Dat1_1m2,  
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
# get letters for Month boxplot
emmeans_m1_neo2 <- cld(emmeans(m1_neo, list(pairwise ~ Month)), 
                        Letters = letters) %>% arrange(Month)
emmeans_m1_neo2

ggplot(Dat1_1m2,aes(x=Month,y=neophyte_prop,col=Month)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=Month), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m1_neo2,
            aes(x=Month, y=c(0.2, 0.2, 0.2, 0.2),
                label=emmeans_m1_neo2$.group),vjust=0.5, hjust=0.6, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Month of sampling" , y="Neophyte probability in %") 




# phen_Richness ---------------------------------------------------------------------------
m1_phen_Richness <-glmer (phen_Richness ~ Month * MowFreq + (1|PlotNo), data = Dat1_1m2,  
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

m2_phen_Richness <-glmer(phen_Richness ~ Month + MowFreq + (1|PlotNo), data = Dat1_1m2,  
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
# get letters for Month boxplot
emmeans_m2_phen_Richness2 <- cld(emmeans(m2_phen_Richness, list(pairwise ~ Month)), 
                                 Letters = letters) %>% arrange(Month)
emmeans_m2_phen_Richness2

ggplot(Dat1_1m2,aes(x=Month,y=phen_Richness,col=Month)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=Month), pch=21,
              size=2, alpha=0.4, stroke = 0.8, 
              position=position_jitterdodge(jitter.width = 0.4, 
                                            jitter.height = 0)) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_phen_Richness2,
            aes(x=Month, y=c(5, 5, 5, 5),
                label=emmeans_m2_phen_Richness2$.group),vjust=0.5, hjust=0.5, 
            size=4, col="black" , position=position_dodge(0)) +
  labs(x="Month of sampling" , y="Phenological richness") 


# phen_evenness  ------------------------------------------------------------------------

m1_ev <-lmer (log(phen_evenness) ~ Month * MowFreq + (1|PlotNo), data = Dat1_1m2)

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
emmeans_m1_ev <- cld(emmeans(m1_ev, list(pairwise ~ MowFreq:Month)), 
                     Letters = letters) %>% arrange(MowFreq, Month)
emmeans_m1_ev


# add to the ggplot

ggplot(Dat1_1m2,aes(x=MowFreq,y=phen_evenness,col=Month)) + 
  geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  geom_point( aes(fill=Month), pch=21,
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

