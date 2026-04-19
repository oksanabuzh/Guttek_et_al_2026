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

library(conflicted)
# Prefer dplyr's select whenever there is a conflict
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Data -------------------------------------------------------------------------

## Plot data ------------------------------------------------------------------

Cover_data <- read_csv("data/raw_data/BC_2025_Cover_Data.csv") %>% 
  filter(Scale_m2 == 1) %>%
  dplyr::select(-c("Date", "Scale_m2", "Remarks", "10m_Max_Cryptogam_Height",
                   "Litter_Cover_litte_from_2025_mowing",
                   "Stones_Rocks_Cover"))

names(Cover_data)

summary(Cover_data)



Mowing_data <- read_csv("data/raw_data/mowing_events_2025_DB.csv") %>% 
  dplyr::select(-mowing_events_2025) %>% 
  pivot_longer(cols = c(September,	July,	May,	March),
               names_to = "Month",
               values_to = "n_mow_events_befre_sampling") %>% 
  relocate(Month, .after=Subplot)



## Diversity data -------------------------------------------------------------
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

## GIS data: ----------------------------------------------------------------

# Data measured in the field:
Fild_data <- read_csv("data/processed_data/GIS_data/GIS_field_measurements.csv") 
# GIS main data (patch type, area, distances):
GIS_main <- read_csv("data/processed_data/GIS_data/GISmain_PatchTypeArea_Distances.csv") 
# Sky-view data:
Sky_view <- read_csv("data/processed_data/GIS_data/Sky_view_factor.csv") 

# GIS data in buffers:
GIS_250 <- read_csv("data/processed_data/GIS_data/GIS_buffer250m.csv")%>% 
  dplyr::select(-protected_biotopes_lines_count,
                -protected_biotopes_points_count,
                -poi_count) %>% 
  rename("Urban_Land_Cover_pct"=ulc_pct)

names(GIS_250)


GIS_500 <- read_csv("data/processed_data/GIS_data/GIS_buffer500m.csv")%>% 
  dplyr::select(-protected_biotopes_lines_count,
                -protected_biotopes_points_count,
                -poi_count) %>% 
  rename("Urban_Land_Cover_pct"=ulc_pct)

GIS_1000 <- read_csv("data/processed_data/GIS_data/GIS_buffer1000m.csv")%>% 
  dplyr::select(-protected_biotopes_lines_count,
                -protected_biotopes_points_count,
                -poi_count) %>% 
  rename("Urban_Land_Cover_pct"=ulc_pct)

Biotop_diversity <- read_csv("data/processed_data/GIS_data/Biotop_Diversity.csv") 




## Merged data ------------------------------------------------------------------

Data_1m2 <- Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, n_mow_events_befre_sampling, 
                SR, biomass, evenness, 
                phen_Richness, phen_evenness, FRic, FEve, FDis, 
                Neophytes_SR_propr, Neophytes_mass_propr) %>% 
  left_join(Cover_data %>% 
            select(PlotNo, Subplot, Month, Bare_Ground_Cover, Litter_Cover_Total,
                   Litter_Cover_last_year_litter) %>% 
              rename(Litter_Cover=Litter_Cover_last_year_litter),
            by = c("PlotNo", "Subplot", "Month")) %>% 
  left_join(Fild_data %>% # mapping field data (by Carolin)
              select(PlotNo, Subplot, slope_degr, dist_tree), 
            by = c("PlotNo", "Subplot")) %>% 
  left_join(GIS_main %>%
              select(PlotNo, patch_biotope_area_sqm) %>% 
              rename(patch_size_m2=patch_biotope_area_sqm), 
            by = c("PlotNo")) %>%
  left_join(Biotop_diversity %>% select(PlotNo, Biotop_richness_specific),
            by = c("PlotNo")) %>% 
  left_join(Sky_view %>% 
              rename(sky_view_factor=skyview_factor_rayman),
            by = c("PlotNo"))  %>% 
  left_join(GIS_500 %>% 
            select(PlotNo, protected_biotopes_polygons_cover_pct, 
                   green_cover_pct, # (correlation with impervious serfices)
                   road_density_km_per_ha) %>% 
              rename(protected_cover_pct=protected_biotopes_polygons_cover_pct),
            by = c("PlotNo")) %>% 
  mutate(protected_cover_pct_log = log1p(protected_cover_pct),
         green_cover_pct_log=log1p(green_cover_pct),
         patch_size_m2_scaled = scale(patch_size_m2)) %>% 
  mutate(MowFreq=fct_relevel(MowFreq,"regular", 
                             "reduced", 
                             "reduced & sowing")) %>% 
  mutate(Month=fct_relevel(Month,"March", "May", "July", "September")) %>% 
  mutate(PlotNo=factor(PlotNo))

names(Dat1_1m2)

# Predictors (selected):

` # plot data
Month * MowFreq +  n_mow_events_befre_sampling +
  # field data from vegetation survays (cover data)
  Bare_Ground_Cover + Litter_Cover + 
  # mapping in field
  slope_degr + dist_tree + sky_view_factor +
  # GIS data
  # patch_size_m2  # correlates strongly with MowFreq
  Biotop_richness_specific + 
  # 500 m buffer
  protected_cover_pct + green_cover_pct + road_density_km_per_ha`

# assigning colors for the plots
MowFreq_col <- c("#F8766D", "#00B0F6","#00BA38")
Month_col <-  c( "orange", "#287271", "#6D326D","brown") 


# 4) Neophyte proportion ----------------------------------------------------------
# Neophytes_mass_propr is the proportion of neophytes in the total biomass, so it is a value between 0 and 1.

## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, Neophytes_mass_propr, SR,
         n_mow_events_befre_sampling,
           Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
           sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, Neophytes_mass_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = log1p(Neophytes_mass_propr))) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes biomass proportion")


## Models: ----------------------------
Data_1m2 %>% 
  select(biomass, Neophytes_mass_propr)
# Neophytes_mass_propr is the proportion of neophytes in the total biomass, so it is a value between 0 and 1.
# beta regression is appropriate for modeling proportions, but it cannot handle values exactly equal to 0 or 1.
# we use zero-inflated beta regression to account for the excess zeros in the data (plots with no neophytes).

library(glmmTMB)
m1_neo <- glmmTMB(Neophytes_mass_propr ~ 
                    # SR + 
                    Month*MowFreq +  
                    n_mow_events_befre_sampling +
                    Litter_Cover + 
                    road_density_km_per_ha + 
                    patch_size_m2_scaled +
                    # can be removed:
                    Bare_Ground_Cover + 
                    slope_degr + dist_tree + 
                    sky_view_factor +
                    Biotop_richness_specific + 
                    log1p(green_cover_pct) + 
                    (1|PlotNo),
                  family = beta_family(),
                  ziformula = ~1,
                  # weights = biomass, # In glmmTMB with beta family, weights are interpreted as precision weights (inverse variance), not as sample sizes like in binomial models. If we would want to weight by total biomass as a measure of reliability, this would work. But this is not what we want.
                  data = Data_1m2)

# check model assumptions
check_convergence(m1_neo)
check_model(m1_neo)
check_overdispersion(m1_neo)

# check interactions
drop1(m1_neo, test = "Chisq")

# interaction is not significant, remove the interaction term
m2_neo <- glmmTMB(Neophytes_mass_propr ~ 
                   # SR + 
                    Month + MowFreq +  
                    n_mow_events_befre_sampling +
                    Litter_Cover + 
                    road_density_km_per_ha + 
                    patch_size_m2_scaled +
                   # can be removed:
                    Bare_Ground_Cover + 
                    slope_degr + dist_tree + 
                    sky_view_factor +
                    # patch_size_m2_scaled +  # correlates with MowFreq
                    Biotop_richness_specific + 
                    log1p(green_cover_pct) + 
                  (1|PlotNo),
                  family = beta_family(),
                  ziformula = ~1,
                  data = Data_1m2)

summary(m2_neo)

# check model assumptions
check_convergence(m2_neo)
check_model(m2_neo)
check_collinearity(m2_neo)
check_overdispersion(m2_neo)

# check predictor effects
drop1(m2_neo, test = "Chisq")
# Anova(m2_neo, type = "II")


## Plots ------------------------------------------------------------------------

##  road_density_km_per_ha  ---------
Data_1m2$road_density_km_per_ha %>% 
  summary()

m2_neo_Road_pred <- ggpredict(m2_neo, 
                              terms = c("road_density_km_per_ha[0.06:0.151, by=0.001]")) %>%
  as.data.frame()

ggplot(m2_neo_Road_pred, aes(x = x, y = predicted)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(road_density_km_per_ha, Neophytes_mass_propr, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Road density (km/ha)", 
       y = "% biomass of neophytes") 

##  Litter_Cover  ---------
Data_1m2$Litter_Cover %>% 
  summary()

m2_neo_Litter_pred <- ggpredict(m2_neo, terms = c("Litter_Cover[2.00:65.00, by=0.01]")) %>%
  as.data.frame()

ggplot(m2_neo_Litter_pred, aes(x = x, y = predicted)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, Neophytes_mass_propr, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Litter cover, %", 
       y = "% biomass of neophytes") 

##  mowing events before sampling  ---------

m2_neo_Mow_pred <- ggpredict(m2_neo, terms = c("n_mow_events_befre_sampling[-0.2:3.2, by=0.0001]")) %>%
  as.data.frame()

ggplot(m2_neo_Mow_pred, aes(x = x, y = predicted)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
               aes(n_mow_events_befre_sampling, Neophytes_mass_propr, color=Month), 
              pch=19, size=1.5, alpha=0.6, stroke = 0.8,
              position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Mowing events before sampling", 
       y = "% biomass of neophytes") 


##  Patch size  ---------
Data_1m2$patch_size_m2_scaled %>% 
  summary()

m2_neo_Patch_pred <- ggpredict(m2_neo, 
                              terms = c("patch_size_m2_scaled[-1.25:2.28, by=0.001]")) %>%
  as.data.frame() %>% 
  mutate(x_tr = x * sd(Data_1m2$patch_size_m2) + mean(Data_1m2$patch_size_m2))

ggplot(m2_neo_Patch_pred, aes(x = x_tr, y = predicted)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(patch_size_m2, Neophytes_mass_propr, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = expression("Patch size, m"^2),
       y = "% biomass of neophytes") 

## MowFreq -------

emmeans(m2_neo, list(pairwise ~ MowFreq))

emmeans_m2_neo_MowFreq1 <- cld(emmeans(m2_neo, list(pairwise ~ MowFreq)), 
                              Letters = letters) %>% 
  arrange(MowFreq)

Neo_max1 <-  Data_1m2 %>% 
  summarise(max=max(Neophytes_mass_propr), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = Neophytes_mass_propr)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "% biomass of neophytes") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_neo_MowFreq1 %>% 
              left_join(Neo_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+0.2,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_neo, list(pairwise ~ Month))

emmeans_m2_neo_Month <- cld(emmeans(m2_neo, list(pairwise ~ Month)), 
                               Letters = letters) %>% 
  arrange(Month)

Neo_max2 <-  Data_1m2 %>% 
  summarise(max=max(Neophytes_mass_propr), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = Neophytes_mass_propr)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "% biomass of neophytes") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_neo_Month %>% 
              left_join(Neo_max2, by=c("Month")),
            aes(x=Month, y=max+0.1,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_neo, list(pairwise ~ MowFreq | Month))

emmeans_m2_neo_MowFreq <- cld(emmeans(m2_neo, list(pairwise ~ MowFreq | Month)), 
                          Letters = letters) %>% 
  arrange(MowFreq)

Neo_max <-  Data_1m2 %>% 
  summarise(max=max(Neophytes_mass_propr), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = Neophytes_mass_propr)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "% biomass of neophytes") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_neo_MowFreq %>% 
              left_join(Neo_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+0.2,
                label=.group),
            size=3.5, col="black") 




# -------------------------------------------------------------------------------
m1_neo <-lmerTest::lmer((Neophytes_mass_propr) ~ 
                          SR + 
                          Month + MowFreq +  
                          n_mow_events_befre_sampling +          
                          Litter_Cover + 
                          Bare_Ground_Cover + 
                          slope_degr + dist_tree + 
                          sky_view_factor +
                          patch_size_m2 + Biotop_richness_specific + 
                          log1p(green_cover_pct) + 
                          road_density_km_per_ha  + 
                          (1|PlotNo), weights = biomass, data = Data_1m2)
## check interactions ------
# from lmerTest package
drop1(m1_neo) # <- use F and  P values to determine which predictors are significant 


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

