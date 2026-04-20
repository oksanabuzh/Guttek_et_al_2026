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
         patch_size_m2_scaled = scale(patch_size_m2),
         road_density_km_per_ha_scaled = scale(road_density_km_per_ha)) %>% 
  mutate(MowFreq=fct_relevel(MowFreq,"regular", 
                             "reduced", 
                             "reduced & sowing")) %>% 
  mutate(Month=fct_relevel(Month,"March", "May", "July", "September")) %>% 
  mutate(PlotNo=factor(PlotNo))

names(Data_1m2)

Data_1m2 %>% 
  write_csv("data/processed_data/Data_1m2_analysis.csv")

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

# 1) Biomass -----------------------------------------------------------------------


## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, biomass,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, biomass), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = biomass)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL)


## Models: ----------------------------
Data_1m2 %>% 
  select(biomass)

m1_biomass <- lmerTest::lmer(biomass ~ 
                               MowFreq*Month +  
                               scale(n_mow_events_befre_sampling) +
                               scale(Litter_Cover) + 
                               scale(road_density_km_per_ha) + 
                               scale(patch_size_m2) +
                               Bare_Ground_Cover + 
                               scale(slope_degr) + scale(dist_tree) + 
                               scale(sky_view_factor) +
                               scale(Biotop_richness_specific) + 
                               scale(green_cover_pct) + 
                               (1|PlotNo),
                             data = Data_1m2)


# check model assumptions
check_convergence(m1_biomass)
check_model(m1_biomass)
check_collinearity(m1_biomass)
# check interactions
drop1(m1_biomass)

# interaction is  significant

# remove patch_size_m2
m2_biomass <- lmerTest::lmer(biomass ~ 
                               MowFreq*Month +  
                               scale(n_mow_events_befre_sampling) +
                               scale(Litter_Cover) + 
                               scale(road_density_km_per_ha) + 
                               # log(patch_size_m2) +
                               Bare_Ground_Cover + 
                               scale(slope_degr) + scale(dist_tree) + 
                               scale(sky_view_factor) +
                               scale(Biotop_richness_specific) + 
                               log1p(green_cover_pct) + 
                               (1|PlotNo),
                             data = Data_1m2)



summary(m2_biomass)
# check model assumptions
check_convergence(m2_biomass)
check_model(m2_biomass)

anova(m1_biomass, m2_biomass)
# m2_biomass

# check predictor effects
drop1(m2_biomass)
# Anova(m2_biomass, type = "II")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_biomass)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_biomass,  partial = T)

Mod_results_biomass <- drop1(m2_biomass) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
 # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m2_biomass,  partial = T) %>% as.data.frame() %>% 
              rename(Driver="Effect") %>% 
              select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "biomass",.before= Driver)

Mod_results_biomass %>% 
  write_csv("results/LMM_biomass.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m2_biomass))


##  road_density_km_per_ha  ---------
Data_1m2$biomass %>% 
  summary()

# m2_biomass_Road_pred <- ggpredict(m2_biomass, terms = c("road_density_km_per_ha")) %>% as.data.frame()
# Very large CI
# "ggeffects" includes both fixed and random effect uncertainty
# "effects" Includes only fixed effect uncertainty
library(effects)
m2_biomass_Road_pred <- Effect("road_density_km_per_ha", m2_biomass) %>% 
  as.data.frame()

ggplot(m2_biomass_Road_pred, aes(x = road_density_km_per_ha, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(road_density_km_per_ha, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Road density (km/ha)", 
       y = "Biomass") 

##  Litter_Cover  ---------
Data_1m2$Litter_Cover %>% 
  summary()

# m2_biomass_Litter_pred <- ggpredict(m2_biomass, terms = c("Litter_Cover[2.00:65.00, by=0.01]")) %>% as.data.frame()
m2_biomass_Litter_pred <- Effect("Litter_Cover", m2_biomass) %>% 
  as.data.frame()

ggplot(m2_biomass_Litter_pred, aes(Litter_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Litter cover, %", 
       y = "Biomass") 

##  mowing events before sampling  ---------

m2_biomass_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_biomass) %>% 
  as.data.frame()

ggplot(m2_biomass_Mow_pred, aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Biomass") 

##  Biotop_richness_specific  ---------
Data_1m2$Biotop_richness_specific %>% 
  summary()

# m2_biomass_Litter_pred <- ggpredict(m2_biomass, terms = c("Biotop_richness_specific")) %>% as.data.frame()
m2_biomass_Biotop_pred <- Effect("Biotop_richness_specific", m2_biomass) %>% 
  as.data.frame()

ggplot(m2_biomass_Biotop_pred, aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Biomass") 


##  dist_tree  ---------
Data_1m2$dist_tree %>% 
  summary()

m2_biomass_Patch_pred <- Effect("dist_tree", m2_biomass) %>% 
  as.data.frame() 

ggplot(m2_biomass_Patch_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Biomass") 



##  green_cover_pct  ---------
Data_1m2$green_cover_pct %>% 
  summary()

m2_biomass_green_pred <- Effect("green_cover_pct", m2_biomass) %>% 
  as.data.frame() 

ggplot(m2_biomass_green_pred, aes(x = green_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(green_cover_pct, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Green cover (%)",
       y = "Biomass") 

## MowFreq -------

emmeans(m2_biomass, list(pairwise ~ MowFreq))

emmeans_m2_biomass_MowFreq1 <- cld(emmeans(m2_biomass, list(pairwise ~ MowFreq)), 
                                   Letters = letters) %>% 
  arrange(MowFreq)

biomass_max1 <-  Data_1m2 %>% 
  summarise(max=max(biomass), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = biomass)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Biomass") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_biomass_MowFreq1 %>% 
              left_join(biomass_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+3,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_biomass, list(pairwise ~ Month))

emmeans_m2_biomass_Month <- cld(emmeans(m2_biomass, list(pairwise ~ Month)), 
                                Letters = letters) %>% 
  arrange(Month)

biomass_max2 <-  Data_1m2 %>% 
  summarise(max=max(biomass), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = biomass)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Biomass") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_biomass_Month %>% 
              left_join(biomass_max2, by=c("Month")),
            aes(x=Month, y=max+5,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_biomass, list(pairwise ~ MowFreq | Month))

emmeans_m2_biomass_MowFreq <- cld(emmeans(m2_biomass, list(pairwise ~ MowFreq | Month)), 
                                  Letters = letters) %>% 
  arrange(MowFreq)

biomass_max <-  Data_1m2 %>% 
  summarise(max=max(biomass), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = biomass)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Biomass") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_biomass_MowFreq %>% 
              left_join(biomass_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+1120,
                label=.group),
            size=3.5, col="black") 





# 2) Species diversity  ----------------------------------------------------------

# 2.1) Species richness -----------------------

## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, SR,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, SR), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = SR)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL)


## Models: ----------------------------
Data_1m2 %>% 
  select(SR)

m1_SR <- glmer(SR ~ 
                 Month * MowFreq +  
                 scale(n_mow_events_befre_sampling) +
                 scale(Litter_Cover) + 
                 scale(road_density_km_per_ha) + 
                 scale(patch_size_m2) +
                 Bare_Ground_Cover + 
                 scale(slope_degr) + scale(dist_tree) + 
                 scale(sky_view_factor) +
                 scale(Biotop_richness_specific) + 
                 scale(green_cover_pct) + 
                 (1|PlotNo),
               family = poisson,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)),  
               data = Data_1m2)


# check model assumptions
check_convergence(m1_SR)
check_model(m1_SR)
check_overdispersion(m1_SR)
check_collinearity(m1_SR)
# check interactions
drop1(m1_SR, test = "Chisq")

# interaction is not significant, remove the interaction term
m2_SR_a <- glmer(SR ~ 
                   Month + MowFreq +  
                   scale(n_mow_events_befre_sampling) +
                   scale(Litter_Cover) + 
                   scale(road_density_km_per_ha) + 
                   scale(patch_size_m2) +
                   scale(Bare_Ground_Cover) + 
                   scale(slope_degr) + scale(dist_tree) + 
                   scale(sky_view_factor) +
                   scale(Biotop_richness_specific) + 
                   log1p(green_cover_pct) + 
                   (1|PlotNo),
                 family = poisson,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)),  
                 data = Data_1m2)

# remove  patch_size_m2
m2_SR <- glmer(SR ~ 
                 Month + MowFreq +  
                 scale(n_mow_events_befre_sampling) +
                 scale(Litter_Cover) + 
                 scale(road_density_km_per_ha) + 
                 #  scale(patch_size_m2) +
                 scale(Bare_Ground_Cover) + 
                 scale(slope_degr) + scale(dist_tree) + 
                 scale(sky_view_factor) +
                 scale(Biotop_richness_specific) + 
                 log1p(green_cover_pct) + 
                 (1|PlotNo),
               family = poisson,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)),  
               data = Data_1m2)


# compare model without patch_size_m2
anova(m2_SR_a, m2_SR)
# keep simpler model

summary(m2_SR)

# check model assumptions
check_convergence(m2_SR)
check_model(m2_SR)
check_collinearity(m2_SR)
check_overdispersion(m2_SR)

# check predictor effects
drop1(m2_SR, test = "Chisq")
# Anova(m2_SR, type = "II")


## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_SR)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_SR,  partial = T)

# r2beta has problem with scale directly in teh model, rerun model
Data_1m2_dummy <-Data_1m2 %>% 
  mutate(
    n_mow_events_scaled = scale(n_mow_events_befre_sampling),
    Litter_Cover_scaled = scale(Litter_Cover),
    road_density_scaled = scale(road_density_km_per_ha),
    Bare_Ground_scaled = scale(Bare_Ground_Cover),
    slope_degr_scaled = scale(slope_degr),
    dist_tree_scaled = scale(dist_tree),
    sky_view_factor_scaled = scale(sky_view_factor),
    Biotop_richness_scaled = scale(Biotop_richness_specific),
    green_cover_log = log1p(green_cover_pct))

m2_SR_dummy <- glmer(SR ~ 
                 Month + MowFreq +  
                 n_mow_events_scaled +
                 Litter_Cover_scaled + 
                 road_density_scaled + 
                 Bare_Ground_scaled + 
                 slope_degr_scaled + dist_tree_scaled + 
                 sky_view_factor_scaled +
                 Biotop_richness_scaled + 
                 green_cover_log + 
                 (1|PlotNo),
               family = poisson,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)),  
               data = Data_1m2_dummy)

r2glmm::r2beta(m2_SR_dummy,  partial = T)

Mod_results_SR <- drop1(m2_SR_dummy, test = "Chisq") %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-AIC) %>% 
  filter(Driver != "<none>") %>% 
  rename("DF"="npar",
         "Chi"= LRT) %>%
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  left_join(
    r2glmm::r2beta(m2_SR_dummy,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "SR", .before= Driver)


Mod_results_SR %>% 
  write_csv("results/GLMM_SpecRich.csv")

## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m2_SR))




##  road_density_km_per_ha  ---------
Data_1m2$SR %>% 
  summary()

# m2_SR_Road_pred <- ggpredict(m2_SR, terms = c("road_density_km_per_ha")) %>% as.data.frame()
# Very large CI
# "ggeffects" includes both fixed and random effect uncertainty
# "effects" Includes only fixed effect uncertainty
library(effects)
m2_SR_Road_pred <- Effect("road_density_km_per_ha", m2_SR) %>% 
  as.data.frame()

ggplot(m2_SR_Road_pred, aes(x = road_density_km_per_ha, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(road_density_km_per_ha, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Road density (km/ha)", 
       y = "Species richness") 


##  Biotop_richness_specific  ---------
Data_1m2$Biotop_richness_specific %>% 
  summary()

# m2_SR_Litter_pred <- ggpredict(m2_SR, terms = c("Biotop_richness_specific")) %>% as.data.frame()
m2_SR_Biotop_pred <- Effect("Biotop_richness_specific", m2_SR) %>% 
  as.data.frame()

ggplot(m2_SR_Biotop_pred, aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Species richness") 


##  slope_degr  ---------
Data_1m2$slope_degr %>% 
  summary()

# m2_SR_Litter_pred <- ggpredict(m2_SR, terms = c("slope_degr[2.00:65.00, by=0.01]")) %>% as.data.frame()
m2_SR_slope_pred <- Effect("slope_degr", m2_SR) %>% 
  as.data.frame()

ggplot(m2_SR_slope_pred, aes(slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope", 
       y = "Species richness") 

##  mowing events before sampling  ---------

m2_SR_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_SR) %>% 
  as.data.frame()

ggplot(m2_SR_Mow_pred, aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Species richness") 



## MowFreq -------

emmeans(m2_SR, list(pairwise ~ MowFreq))

emmeans_m2_SR_MowFreq1 <- cld(emmeans(m2_SR, list(pairwise ~ MowFreq)), 
                              Letters = letters) %>% 
  arrange(MowFreq)

SR_max1 <-  Data_1m2 %>% 
  summarise(max=max(SR), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = SR)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Species richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_SR_MowFreq1 %>% 
              left_join(SR_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+3,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_SR, list(pairwise ~ Month))

emmeans_m2_SR_Month <- cld(emmeans(m2_SR, list(pairwise ~ Month)), 
                           Letters = letters) %>% 
  arrange(Month)

SR_max2 <-  Data_1m2 %>% 
  summarise(max=max(SR), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = SR)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Species richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_SR_Month %>% 
              left_join(SR_max2, by=c("Month")),
            aes(x=Month, y=max+5,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_SR, list(pairwise ~ MowFreq | Month))

emmeans_m2_SR_MowFreq <- cld(emmeans(m2_SR, list(pairwise ~ MowFreq | Month)), 
                             Letters = letters) %>% 
  arrange(MowFreq)

SR_max <-  Data_1m2 %>% 
  summarise(max=max(SR), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = SR)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Species richness") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_SR_MowFreq %>% 
              left_join(SR_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+5,
                label=.group),
            size=3.5, col="black") 




# 2.2) evenness -----------------------------------------------------------------------


## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, evenness,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, evenness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = evenness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Biomass")


## Models: ----------------------------
Data_1m2 %>% 
  select(evenness)

m1_evenness <- lmerTest::lmer(evenness ~ 
                                MowFreq*Month +  
                                scale(n_mow_events_befre_sampling) +
                                scale(Litter_Cover) + 
                                scale(road_density_km_per_ha) + 
                                scale(patch_size_m2) +
                                Bare_Ground_Cover + 
                                scale(slope_degr) + scale(dist_tree) + 
                                scale(sky_view_factor) +
                                scale(Biotop_richness_specific) + 
                                scale(green_cover_pct) + 
                                (1|PlotNo),
                              data = Data_1m2)


# check model assumptions
check_convergence(m1_evenness)
check_model(m1_evenness)
check_collinearity(m1_evenness)
# check interactions
drop1(m1_evenness)

# remove interaction 
# remove patch_size_m2

m2_evenness_a <- lmerTest::lmer(evenness ~ 
                                  MowFreq+Month +  
                                  scale(n_mow_events_befre_sampling) +
                                  scale(Litter_Cover) + 
                                  scale(road_density_km_per_ha) + 
                                  # scale(patch_size_m2) +
                                  Bare_Ground_Cover + 
                                  scale(slope_degr) + scale(dist_tree) + 
                                  scale(sky_view_factor) +
                                  scale(Biotop_richness_specific) + 
                                  log1p(green_cover_pct) + 
                                  (1|PlotNo),
                                data = Data_1m2)

# keep patch_size_m2
m2_evenness <- lmerTest::lmer(evenness ~ 
                                MowFreq+Month +  
                                scale(n_mow_events_befre_sampling) +
                                scale(Litter_Cover) + 
                                scale(road_density_km_per_ha) + 
                                scale(patch_size_m2) +
                                Bare_Ground_Cover + 
                                scale(slope_degr) + scale(dist_tree) + 
                                scale(sky_view_factor) +
                                scale(Biotop_richness_specific) + 
                                log1p(green_cover_pct) + 
                                (1|PlotNo),
                              data = Data_1m2)



summary(m2_evenness)
# check model assumptions
check_convergence(m2_evenness)
check_model(m2_evenness)

anova(m2_evenness_a, m2_evenness)
# m2_evenness # with patch_size_m2

# check predictor effects
drop1(m2_evenness)
# Anova(m1_evenness, type = "II")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_evenness)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_evenness,  partial = T)

Mod_results_evenness <- drop1(m2_evenness) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m2_evenness,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "evenness",.before= Driver)

Mod_results_evenness %>% 
  write_csv("results/LMM_evenness.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m2_evenness))


##  road_density_km_per_ha  ---------
Data_1m2$evenness %>% 
  summary()

# m2_evenness_Road_pred <- ggpredict(m2_evenness, terms = c("road_density_km_per_ha")) %>% as.data.frame()
# Very large CI
# "ggeffects" includes both fixed and random effect uncertainty
# "effects" Includes only fixed effect uncertainty
library(effects)
m2_evenness_Road_pred <- Effect("road_density_km_per_ha", m2_evenness) %>% 
  as.data.frame()

ggplot(m2_evenness_Road_pred, aes(x = road_density_km_per_ha, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(road_density_km_per_ha, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Road density (km/ha)", 
       y = "Hill-Simpson diversity") 

##  patch_size_m2  ---------
Data_1m2$patch_size_m2 %>% 
  summary()

# m2_evenness_Litter_pred <- ggpredict(m2_evenness, terms = c("patch_size_m2[2.00:65.00, by=0.01]")) %>% as.data.frame()
m2_evenness_Litter_pred <- Effect("patch_size_m2", m2_evenness) %>% 
  as.data.frame()

ggplot(m2_evenness_Litter_pred, aes(patch_size_m2, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(patch_size_m2, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = expression("Patch size, m"^2),
       y = "Hill-Simpson diversity") 

##  mowing events before sampling  ---------

m2_evenness_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_evenness) %>% 
  as.data.frame()

ggplot(m2_evenness_Mow_pred, aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Hill-Simpson diversity") 

##  Biotop_richness_specific  ---------
Data_1m2$Biotop_richness_specific %>% 
  summary()

# m2_evenness_Litter_pred <- ggpredict(m2_evenness, terms = c("Biotop_richness_specific")) %>% as.data.frame()
m2_evenness_Biotop_pred <- Effect("Biotop_richness_specific", m2_evenness) %>% 
  as.data.frame()

ggplot(m2_evenness_Biotop_pred, aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Hill-Simpson diversity") 


##  slope_degr  ---------
Data_1m2$slope_degr %>% 
  summary()

m2_evenness_Patch_pred <- Effect("slope_degr", m2_evenness) %>% 
  as.data.frame() 

ggplot(m2_evenness_Patch_pred, aes(x = slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope",
       y = "Hill-Simpson diversity") 



##  sky_view_factor  ---------
Data_1m2$sky_view_factor %>% 
  summary()

m2_evenness_green_pred <- Effect("sky_view_factor", m2_evenness) %>% 
  as.data.frame() 

ggplot(m2_evenness_green_pred, aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Hill-Simpson diversity") 

## MowFreq -------

emmeans(m2_evenness, list(pairwise ~ MowFreq))

emmeans_m2_evenness_MowFreq1 <- cld(emmeans(m2_evenness, list(pairwise ~ MowFreq)), 
                                    Letters = letters) %>% 
  arrange(MowFreq)

evenness_max1 <-  Data_1m2 %>% 
  summarise(max=max(evenness), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = evenness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Hill-Simpson diversity") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_evenness_MowFreq1 %>% 
              left_join(evenness_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+3,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_evenness, list(pairwise ~ Month))

emmeans_m2_evenness_Month <- cld(emmeans(m2_evenness, list(pairwise ~ Month)), 
                                 Letters = letters) %>% 
  arrange(Month)

evenness_max2 <-  Data_1m2 %>% 
  summarise(max=max(evenness), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = evenness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Hill-Simpson diversity") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_evenness_Month %>% 
              left_join(evenness_max2, by=c("Month")),
            aes(x=Month, y=max+5,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_evenness, list(pairwise ~ MowFreq | Month))

emmeans_m2_evenness_MowFreq <- cld(emmeans(m2_evenness, list(pairwise ~ MowFreq | Month)), 
                                   Letters = letters) %>% 
  arrange(MowFreq)

evenness_max <-  Data_1m2 %>% 
  summarise(max=max(evenness), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = evenness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Hill-Simpson diversity") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_evenness_MowFreq %>% 
              left_join(evenness_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+5,
                label=.group),
            size=3.5, col="black") 

# 3) Phenological DIversity -----------------------------------------------------------------------


# 3.1) Phenological richness -----------------------------------------------------------------------

## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, phen_Richness,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, phen_Richness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Richness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Biomass")


## Models: ----------------------------
Data_1m2 %>% 
  select(phen_Richness)

m1_phen_Richness <- lmerTest::lmer(phen_Richness ~ 
                                     MowFreq*Month +  
                                     scale(n_mow_events_befre_sampling) +
                                     scale(Litter_Cover) + 
                                     scale(road_density_km_per_ha) + 
                                     scale(patch_size_m2) +
                                     Bare_Ground_Cover + 
                                     scale(slope_degr) + scale(dist_tree) + 
                                     scale(sky_view_factor) +
                                     scale(Biotop_richness_specific) + 
                                     scale(green_cover_pct) + 
                                     (1|PlotNo),
                                   data = Data_1m2)


# check model assumptions
check_convergence(m1_phen_Richness)
check_model(m1_phen_Richness)
check_collinearity(m1_phen_Richness)
# check interactions
drop1(m1_phen_Richness)

# remove patch_size_m2

m2_phen_Richness <- lmerTest::lmer(phen_Richness ~ 
                                     MowFreq*Month +  
                                     scale(n_mow_events_befre_sampling) +
                                     scale(Litter_Cover) + 
                                     scale(road_density_km_per_ha) + 
                                     # scale(patch_size_m2) +
                                     Bare_Ground_Cover + 
                                     scale(slope_degr) + scale(dist_tree) + 
                                     scale(sky_view_factor) +
                                     scale(Biotop_richness_specific) + 
                                     log1p(green_cover_pct) + 
                                     (1|PlotNo),
                                   data = Data_1m2)


summary(m2_phen_Richness)
# check model assumptions
check_convergence(m2_phen_Richness)
check_model(m2_phen_Richness)

anova(m1_phen_Richness, m2_phen_Richness)
# m2_phen_Richness # with patch_size_m2

# check predictor effects
drop1(m2_phen_Richness)
# Anova(m1_phen_Richness, type = "II")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_phen_Richness)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_phen_Richness,  partial = T)

Mod_results_Phen_Richness <- drop1(m2_phen_Richness) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m2_phen_Richness,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "phen_Richness",.before= Driver)

Mod_results_Phen_Richness %>% 
  write_csv("results/LMM_Phenol_Richness.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m2_phen_Richness))

##  mowing events before sampling  ---------

m2_phen_Richness_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_phen_Richness) %>% 
  as.data.frame()

ggplot(m2_phen_Richness_Mow_pred, aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, phen_Richness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.2, height = 0.15)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Phenological richness") 

## MowFreq -------

emmeans(m2_phen_Richness, list(pairwise ~ MowFreq))

emmeans_m2_phen_Richness_MowFreq1 <- cld(emmeans(m2_phen_Richness, list(pairwise ~ MowFreq)), 
                                         Letters = letters) %>% 
  arrange(MowFreq)

phen_Richness_max1 <-  Data_1m2 %>% 
  summarise(max=max(phen_Richness), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = phen_Richness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Phenological richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_phen_Richness_MowFreq1 %>% 
              left_join(phen_Richness_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+1,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_phen_Richness, list(pairwise ~ Month))

emmeans_m2_phen_Richness_Month <- cld(emmeans(m2_phen_Richness, list(pairwise ~ Month)), 
                                      Letters = letters) %>% 
  arrange(Month)

phen_Richness_max2 <-  Data_1m2 %>% 
  summarise(max=max(phen_Richness), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = phen_Richness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Phenological richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_phen_Richness_Month %>% 
              left_join(phen_Richness_max2, by=c("Month")),
            aes(x=Month, y=max+1,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_phen_Richness, list(pairwise ~ MowFreq | Month))

emmeans_m2_phen_Richness_MowFreq <- cld(emmeans(m2_phen_Richness, list(pairwise ~ MowFreq | Month)), 
                                        Letters = letters) %>% 
  arrange(MowFreq) %>% 
  mutate(.group = ifelse(Month == "May" & MowFreq== "reduced & sowing", 
                         "b", .group)) # very close to significant difference

phen_Richness_max <-  Data_1m2 %>% 
  summarise(max=max(phen_Richness), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = phen_Richness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Phenological richness") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_phen_Richness_MowFreq %>% 
              left_join(phen_Richness_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+1,
                label=.group),
            size=3.5, col="black") 



# 3.2) Phenological evenness  -----------------------------------------------------------------------


## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, phen_evenness,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, phen_evenness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_evenness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL)


## Models: ----------------------------
Data_1m2 %>% 
  select(phen_evenness)

m1_phen_evenness <- lmerTest::lmer(phen_evenness ~ 
                                     MowFreq*Month +  
                                     scale(n_mow_events_befre_sampling) +
                                     scale(Litter_Cover) + 
                                     scale(road_density_km_per_ha) + 
                                     scale(patch_size_m2) +
                                     Bare_Ground_Cover + 
                                     scale(slope_degr) + scale(dist_tree) + 
                                     scale(sky_view_factor) +
                                     scale(Biotop_richness_specific) + 
                                     scale(green_cover_pct) + 
                                     (1|PlotNo),
                                   data = Data_1m2)


# check model assumptions
check_convergence(m1_phen_evenness)
check_model(m1_phen_evenness)
check_collinearity(m1_phen_evenness)
# check interactions
drop1(m1_phen_evenness)

# interaction is  significant

# remove patch_size_m2
m2_phen_evenness <- lmerTest::lmer(phen_evenness ~ 
                                     MowFreq*Month +  
                                     scale(n_mow_events_befre_sampling) +
                                     scale(Litter_Cover) + 
                                     scale(road_density_km_per_ha) + 
                                     # log(patch_size_m2) +
                                     Bare_Ground_Cover + 
                                     scale(slope_degr) + scale(dist_tree) + 
                                     scale(sky_view_factor) +
                                     scale(Biotop_richness_specific) + 
                                     log1p(green_cover_pct) + 
                                     (1|PlotNo),
                                   data = Data_1m2)



summary(m2_phen_evenness)
# check model assumptions
check_convergence(m2_phen_evenness)
check_model(m2_phen_evenness)

anova(m1_phen_evenness, m2_phen_evenness)
# m2_phen_evenness

# check predictor effects
drop1(m2_phen_evenness)
# Anova(m2_phen_evenness, type = "II")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_phen_evenness)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_phen_evenness,  partial = T)

Mod_results_Phen_evenness <- drop1(m2_phen_evenness) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m2_phen_evenness,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "phen_evenness",.before= Driver)

Mod_results_Phen_evenness %>% 
  write_csv("results/LMM_Phen_evenness.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m2_phen_evenness))


##  road_density_km_per_ha  ---------
Data_1m2$phen_evenness %>% 
  summary()

# m2_phen_evenness_Road_pred <- ggpredict(m2_phen_evenness, terms = c("road_density_km_per_ha")) %>% as.data.frame()
# Very large CI
# "ggeffects" includes both fixed and random effect uncertainty
# "effects" Includes only fixed effect uncertainty
library(effects)
m2_phen_evenness_Road_pred <- Effect("road_density_km_per_ha", m2_phen_evenness) %>% 
  as.data.frame()

ggplot(m2_phen_evenness_Road_pred, aes(x = road_density_km_per_ha, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(road_density_km_per_ha, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Road density (km/ha)", 
       y = "Phenological evenness") 



##  mowing events before sampling  ---------

m2_phen_evenness_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_phen_evenness) %>% 
  as.data.frame()

ggplot(m2_phen_evenness_Mow_pred, aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Phenological evenness") 

##  Biotop_richness_specific  ---------
Data_1m2$Biotop_richness_specific %>% 
  summary()

# m2_phen_evenness_Litter_pred <- ggpredict(m2_phen_evenness, terms = c("Biotop_richness_specific")) %>% as.data.frame()
m2_phen_evenness_Biotop_pred <- Effect("Biotop_richness_specific", m2_phen_evenness) %>% 
  as.data.frame()

ggplot(m2_phen_evenness_Biotop_pred, aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Phenological evenness") 


##  dist_tree  ---------
Data_1m2$dist_tree %>% 
  summary()

m2_phen_evenness_Patch_pred <- Effect("dist_tree", m2_phen_evenness) %>% 
  as.data.frame() 

ggplot(m2_phen_evenness_Patch_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Phenological evenness") 


## MowFreq -------

emmeans(m2_phen_evenness, list(pairwise ~ MowFreq))

emmeans_m2_phen_evenness_MowFreq1 <- cld(emmeans(m2_phen_evenness, list(pairwise ~ MowFreq)), 
                                         Letters = letters) %>% 
  arrange(MowFreq)

phen_evenness_max1 <-  Data_1m2 %>% 
  summarise(max=max(phen_evenness), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = phen_evenness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Phenological evenness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_phen_evenness_MowFreq1 %>% 
              left_join(phen_evenness_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+3,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_phen_evenness, list(pairwise ~ Month))

emmeans_m2_phen_evenness_Month <- cld(emmeans(m2_phen_evenness, list(pairwise ~ Month)), 
                                      Letters = letters) %>% 
  arrange(Month)

phen_evenness_max2 <-  Data_1m2 %>% 
  summarise(max=max(phen_evenness), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = phen_evenness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Phenological evenness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_phen_evenness_Month %>% 
              left_join(phen_evenness_max2, by=c("Month")),
            aes(x=Month, y=max+5,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_phen_evenness, list(pairwise ~ MowFreq | Month))

emmeans_m2_phen_evenness_MowFreq <- cld(emmeans(m2_phen_evenness, list(pairwise ~ MowFreq | Month)), 
                                        Letters = letters) %>% 
  arrange(MowFreq)

phen_evenness_max <-  Data_1m2 %>% 
  summarise(max=max(phen_evenness), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = phen_evenness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Phenological evenness") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_phen_evenness_MowFreq %>% 
              left_join(phen_evenness_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+0.5,
                label=.group),
            size=3.5, col="black") 


# 4) Functional Diversity -----------------------------------------------------------------------

# 4.1) Functional Richness -----------------------------------------------------------------------

## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, FRic,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, FRic), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FRic)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL)


## Models: ----------------------------
Data_1m2 %>% 
  select(FRic)

m1_FRic <- lmerTest::lmer(FRic ~ 
                            MowFreq*Month +  
                            scale(n_mow_events_befre_sampling) +
                            scale(Litter_Cover) + 
                            scale(road_density_km_per_ha) + 
                            scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(Biotop_richness_specific) + 
                            scale(green_cover_pct) + 
                            (1|PlotNo),
                          data = Data_1m2)


# check model assumptions
check_convergence(m1_FRic)
check_model(m1_FRic)
check_collinearity(m1_FRic)
# check interactions
drop1(m1_FRic)

# interaction is  significant

# remove patch_size_m2
m2_FRic <- lmerTest::lmer(FRic ~ 
                            MowFreq*Month +  
                            scale(n_mow_events_befre_sampling) +
                            scale(Litter_Cover) + 
                            scale(road_density_km_per_ha) + 
                            # log(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(Biotop_richness_specific) + 
                            log1p(green_cover_pct) + 
                            (1|PlotNo),
                          data = Data_1m2)



summary(m2_FRic)
# check model assumptions
check_convergence(m2_FRic)
check_model(m2_FRic)

anova(m1_FRic, m2_FRic)
# m2_FRic

# check predictor effects
drop1(m2_FRic)
# Anova(m2_FRic, type = "II")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_FRic)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_FRic,  partial = T)

Mod_results_FRic <- drop1(m2_FRic) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m2_FRic,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "FRic",.before= Driver)

Mod_results_FRic %>% 
  write_csv("results/LMM_Functional_richness.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m2_FRic))


##  mowing events before sampling  ---------

m2_FRic_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_FRic) %>% 
  as.data.frame()

ggplot(m2_FRic_Mow_pred, aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Functional richness") 




##  Bare_Ground_Cover  ---------
Data_1m2$Bare_Ground_Cover %>% 
  summary()

m2_FRic_green_pred <- Effect("Bare_Ground_Cover", m2_FRic) %>% 
  as.data.frame() 

ggplot(m2_FRic_green_pred, aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground (%)",
       y = "Functional richness") 

## MowFreq -------

emmeans(m2_FRic, list(pairwise ~ MowFreq))

emmeans_m2_FRic_MowFreq1 <- cld(emmeans(m2_FRic, list(pairwise ~ MowFreq)), 
                                Letters = letters) %>% 
  arrange(MowFreq)

FRic_max1 <-  Data_1m2 %>% 
  summarise(max=max(FRic), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = FRic)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Functional richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FRic_MowFreq1 %>% 
              left_join(FRic_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+3,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_FRic, list(pairwise ~ Month))

emmeans_m2_FRic_Month <- cld(emmeans(m2_FRic, list(pairwise ~ Month)), 
                             Letters = letters) %>% 
  arrange(Month)

FRic_max2 <-  Data_1m2 %>% 
  summarise(max=max(FRic), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = FRic)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Functional richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FRic_Month %>% 
              left_join(FRic_max2, by=c("Month")),
            aes(x=Month, y=max+5,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_FRic, list(pairwise ~ MowFreq | Month))

emmeans_m2_FRic_MowFreq <- cld(emmeans(m2_FRic, list(pairwise ~ MowFreq | Month)), 
                               Letters = letters) %>% 
  arrange(MowFreq)

FRic_max <-  Data_1m2 %>% 
  summarise(max=max(FRic), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = FRic)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Functional richness") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FRic_MowFreq %>% 
              left_join(FRic_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+4,
                label=.group),
            size=3.5, col="black") 




# 4.2) Functional evenness -----------------------------------------------------------------------


## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, FEve,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, FEve), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FEve)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL)


## Models: ----------------------------
Data_1m2 %>% 
  select(FEve)

m1_FEve <- lmerTest::lmer(FEve ~ 
                            MowFreq*Month +  
                            scale(n_mow_events_befre_sampling) +
                            scale(Litter_Cover) + 
                            scale(road_density_km_per_ha) + 
                            scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(Biotop_richness_specific) + 
                            scale(green_cover_pct) + 
                            (1|PlotNo),
                          data = Data_1m2)


# check model assumptions
check_convergence(m1_FEve)
check_model(m1_FEve)
check_collinearity(m1_FEve)
# check interactions
drop1(m1_FEve)

# interaction is not significant
m2_FEve_a <- lmerTest::lmer(FEve ~ 
                              MowFreq+Month +  
                              scale(n_mow_events_befre_sampling) +
                              scale(Litter_Cover) + 
                              scale(road_density_km_per_ha) + 
                              scale(patch_size_m2) +
                              Bare_Ground_Cover + 
                              scale(slope_degr) + scale(dist_tree) + 
                              scale(sky_view_factor) +
                              scale(Biotop_richness_specific) + 
                              log1p(green_cover_pct) + 
                              (1|PlotNo),
                            data = Data_1m2)
# remove patch_size_m2
m2_FEve <- lmerTest::lmer(FEve ~ 
                            MowFreq+Month +  
                            scale(n_mow_events_befre_sampling) +
                            scale(Litter_Cover) + 
                            scale(road_density_km_per_ha) + 
                            # scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(Biotop_richness_specific) + 
                            log1p(green_cover_pct) + 
                            (1|PlotNo),
                          data = Data_1m2)



anova(m2_FEve, m2_FEve_a)
# check model assumptions
check_convergence(m2_FEve)
check_model(m2_FEve)


# m2_FEve

# check predictor effects
drop1(m2_FEve)
# Anova(m2_FEve, type = "II")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_FEve)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_FEve,  partial = T)

Mod_results_FEve <- drop1(m2_FEve) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m2_FEve,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "FEve",.before= Driver)

Mod_results_FEve %>% 
  write_csv("results/LMM_Funct_evenness.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m2_FEve))

##  mowing events before sampling  ---------

m2_FEve_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_FEve) %>% 
  as.data.frame()

ggplot(m2_FEve_Mow_pred, aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Functional evenness") 


## MowFreq -------

emmeans(m2_FEve, list(pairwise ~ MowFreq))

emmeans_m2_FEve_MowFreq1 <- cld(emmeans(m2_FEve, list(pairwise ~ MowFreq)), 
                                Letters = letters) %>% 
  arrange(MowFreq) %>% 
  mutate(.group = ifelse(MowFreq == "reduced", "b", .group)) %>% 
  mutate(.group = ifelse(MowFreq == "reduced & sowing", "ab", .group))

FEve_max1 <-  Data_1m2 %>% 
  summarise(max=max(FEve), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = FEve)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Functional evenness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FEve_MowFreq1 %>% 
              left_join(FEve_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+0.1,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_FEve, list(pairwise ~ Month))

emmeans_m2_FEve_Month <- cld(emmeans(m2_FEve, list(pairwise ~ Month)), 
                             Letters = letters) %>% 
  arrange(Month)

FEve_max2 <-  Data_1m2 %>% 
  summarise(max=max(FEve), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = FEve)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Functional evenness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FEve_Month %>% 
              left_join(FEve_max2, by=c("Month")),
            aes(x=Month, y=max+0.1,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_FEve, list(pairwise ~ MowFreq | Month))

emmeans_m2_FEve_MowFreq <- cld(emmeans(m2_FEve, list(pairwise ~ MowFreq | Month)), 
                               Letters = letters) %>% 
  arrange(MowFreq)

FEve_max <-  Data_1m2 %>% 
  summarise(max=max(FEve), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = FEve)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Functional evenness") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FEve_MowFreq %>% 
              left_join(FEve_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+1120,
                label=.group),
            size=3.5, col="black") 



# 4.3) Functional dispersion -----------------------------------------------------------------------


## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, FDis,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, road_density_km_per_ha, protected_cover_pct) %>% 
  pivot_longer(-c(Month, FDis), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FDis)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL)


## Models: ----------------------------
Data_1m2 %>% 
  select(FDis)

m1_FDis <- lmerTest::lmer(FDis ~ 
                            MowFreq*Month +  
                            scale(n_mow_events_befre_sampling) +
                            scale(Litter_Cover) + 
                            scale(road_density_km_per_ha) + 
                            scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(Biotop_richness_specific) + 
                            scale(green_cover_pct) + 
                            (1|PlotNo),
                          data = Data_1m2)


# check model assumptions
check_convergence(m1_FDis)
check_model(m1_FDis)
check_collinearity(m1_FDis)
# check interactions
drop1(m1_FDis)

# interaction is  not significant
m2_FDis_a <- lmerTest::lmer(FDis ~ 
                              MowFreq + Month +  
                              scale(n_mow_events_befre_sampling) +
                              scale(Litter_Cover) + 
                              scale(road_density_km_per_ha) + 
                              scale(patch_size_m2) +
                              Bare_Ground_Cover + 
                              scale(slope_degr) + scale(dist_tree) + 
                              scale(sky_view_factor) +
                              scale(Biotop_richness_specific) + 
                              log1p(green_cover_pct) + 
                              (1|PlotNo),
                            data = Data_1m2)

# remove patch_size_m2
m2_FDis <- lmerTest::lmer(FDis ~ 
                            MowFreq + Month +  
                            scale(n_mow_events_befre_sampling) +
                            scale(Litter_Cover) + 
                            scale(road_density_km_per_ha) + 
                            # scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(Biotop_richness_specific) + 
                            log1p(green_cover_pct) + 
                            (1|PlotNo),
                          data = Data_1m2)

anova(m2_FDis_a, m2_FDis)

summary(m2_FDis)
# check model assumptions
check_convergence(m2_FDis)
check_model(m2_FDis)

# m2_FDis

# check predictor effects
drop1(m2_FDis)
# Anova(m2_FDis, type = "II")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_FDis)
# Partial R2 for fixed effects
r2glmm::r2beta(m2_FDis,  partial = T)

Mod_results_FDis <- drop1(m2_FDis) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m2_FDis,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "FDis",.before= Driver)

Mod_results_FDis %>% 
  write_csv("results/LMM_Functional_dispersion.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m2_FDis))


##  mowing events before sampling  ---------

m2_FDis_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_FDis) %>% 
  as.data.frame()

ggplot(m2_FDis_Mow_pred, aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Functional dispersion") 

##  Biotop_richness_specific  ---------
Data_1m2$Biotop_richness_specific %>% 
  summary()

# m2_FDis_Litter_pred <- ggpredict(m2_FDis, terms = c("Biotop_richness_specific")) %>% as.data.frame()
m2_FDis_Biotop_pred <- Effect("Biotop_richness_specific", m2_FDis) %>% 
  as.data.frame()

ggplot(m2_FDis_Biotop_pred, aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Functional dispersion") 


##  dist_tree  ---------
Data_1m2$dist_tree %>% 
  summary()

m2_FDis_Patch_pred <- Effect("dist_tree", m2_FDis) %>% 
  as.data.frame() 

ggplot(m2_FDis_Patch_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Functional dispersion") 



##  Bare_Ground_Cover  ---------
Data_1m2$Bare_Ground_Cover %>% 
  summary()

m2_FDis_green_pred <- Effect("Bare_Ground_Cover", m2_FDis) %>% 
  as.data.frame() 

ggplot(m2_FDis_green_pred, aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Functional dispersion") 

## MowFreq -------

emmeans(m2_FDis, list(pairwise ~ MowFreq))

emmeans_m2_FDis_MowFreq1 <- cld(emmeans(m2_FDis, list(pairwise ~ MowFreq)), 
                                Letters = letters) %>% 
  arrange(MowFreq)

FDis_max1 <-  Data_1m2 %>% 
  summarise(max=max(FDis), .by = c(MowFreq))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = FDis)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Functional dispersion") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FDis_MowFreq1 %>% 
              left_join(FDis_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+0.5,
                label=.group),
            size=3.5, col="black") 


## Month -------

emmeans(m2_FDis, list(pairwise ~ Month))

emmeans_m2_FDis_Month <- cld(emmeans(m2_FDis, list(pairwise ~ Month)), 
                             Letters = letters) %>% 
  arrange(Month)

FDis_max2 <-  Data_1m2 %>% 
  summarise(max=max(FDis), .by = c(Month))

Data_1m2 %>% 
  ggplot(aes(x = Month, y = FDis)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Functional dispersion") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FDis_Month %>% 
              left_join(FDis_max2, by=c("Month")),
            aes(x=Month, y=max+0.5,
                label=.group),
            size=3.5, col="black")

## MowFreq * Month  -------

emmeans(m2_FDis, list(pairwise ~ MowFreq | Month))

emmeans_m2_FDis_MowFreq <- cld(emmeans(m2_FDis, list(pairwise ~ MowFreq | Month)), 
                               Letters = letters) %>% 
  arrange(MowFreq)

FDis_max <-  Data_1m2 %>% 
  summarise(max=max(FDis), .by = c(MowFreq, Month))

Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = FDis)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Functional dispersion") +
  facet_wrap(~Month) +
  scale_color_manual(values = Month_col) +
  theme(legend.position = "none") +
  geom_text(data=emmeans_m2_FDis_MowFreq %>% 
              left_join(FDis_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+0.5,
                label=.group),
            size=3.5, col="black") 


# 5) Neophyte proportion ----------------------------------------------------------
# Neophytes_mass_propr is the proportion of neophytes in the total biomass, so it is a value between 0 and 1.

## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, Neophytes_mass_propr, 
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
                    Bare_Ground_Cover + 
                    slope_degr + dist_tree + 
                    sky_view_factor +
                    Biotop_richness_specific + 
                    log1p(green_cover_pct) +
                  (1|PlotNo),
                  family = beta_family(),
                  ziformula = ~1,
                  data = Data_1m2)


# remove patch_size_m2_scaled
m2_neo_b <- glmmTMB(Neophytes_mass_propr ~ 
                    # SR + 
                    Month + MowFreq +  
                    n_mow_events_befre_sampling +
                    Litter_Cover + 
                    road_density_km_per_ha + 
                  #  patch_size_m2_scaled +
                    Bare_Ground_Cover + 
                    slope_degr + dist_tree + 
                    sky_view_factor +
                    Biotop_richness_specific + 
                    log1p(green_cover_pct) +
                    (1|PlotNo),
                  family = beta_family(),
                  ziformula = ~1,
                  data = Data_1m2)


anova(m2_neo, m2_neo_b)

summary(m2_neo)

# check model assumptions
check_convergence(m2_neo)
check_model(m2_neo)
check_collinearity(m2_neo)
check_overdispersion(m2_neo)

# model diagnostics based on the simulated residuals
m2_neo_simres <- DHARMa::simulateResiduals(m2_neo)
plot(m2_neo_simres)


# check predictor effects
drop1(m2_neo, test = "Chisq")
# Anova(m2_neo, type = "II")


## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m2_neo)
performance::r2(m2_neo)
# Partial R2 do not work with glmmTMB
# r2glmm::r2beta(m2_neo,  partial = T)
# No good way to calculate for beta glmmTMB


Mod_results_neo <- drop1(m2_neo, test = "Chisq") %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-AIC) %>% 
  filter(Driver != "<none>") %>% 
  rename("Chi"= LRT) 
# NO Partial R2
  # relocate raw "MowFreq:Month" in column "Driver" to the top
#  left_join(
#    r2glmm::r2beta(m2_neo,  partial = T) %>% as.data.frame() %>% 
#      rename(Driver="Effect") %>% 
#      select(Driver,  Rsq), by = "Driver") %>% 
#  mutate(Responce = "Neoph", .before= Driver)


Mod_results_neo %>% 
  write_csv("results/GLMM_Neophite_perc.csv")



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
