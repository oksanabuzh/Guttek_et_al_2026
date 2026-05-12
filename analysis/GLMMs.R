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
conflict_prefer("lmer", "lmerTest")

# R version
R.version.string

# Citation for R
citation()

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
                   impervious_pct,
                   road_density_km_per_ha) %>% 
              rename(protected_cover_pct=protected_biotopes_polygons_cover_pct),
            by = c("PlotNo")) %>% 
  mutate(protected_cover_pct_log = log1p(protected_cover_pct),
         green_cover_pct_log=log1p(green_cover_pct),
         patch_size_m2_scaled = as.vector(scale(patch_size_m2)),
         road_density_km_per_ha_scaled = as.vector(scale(road_density_km_per_ha))) %>% 
  mutate(MowFreq=fct_relevel(MowFreq,"regular", 
                             "reduced", 
                             "reduced & sowing")) %>% 
  mutate(Month=fct_relevel(Month,"March", "May", "July", "September")) %>% 
  mutate(PlotNo=factor(PlotNo))

names(Data_1m2)

Data_1m2 %>% 
  write_csv("data/processed_data/Data_1m2_analysis.csv")


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
         green_cover_pct_log, impervious_pct, protected_cover_pct) %>% 
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
                               scale(log1p(n_mow_events_befre_sampling)) +
                               scale(Litter_Cover) + 
                         #      scale(road_density_km_per_ha) + 
                         #      scale(patch_size_m2) +
                               Bare_Ground_Cover + 
                               scale(slope_degr) + scale(dist_tree) + 
                               scale(sky_view_factor) +
                               scale(Biotop_richness_specific) + 
                          scale(protected_cover_pct) + 
                           #    scale(green_cover_pct) + 
                        # scale(impervious_pct) +
                               (1|PlotNo),
                             data = Data_1m2)


# check model assumptions
check_convergence(m1_biomass)
check_model(m1_biomass)
check_collinearity(m1_biomass)
# check interactions
drop1(m1_biomass)

# interaction is  significant


## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m1_biomass)
# Partial R2 for fixed effects
r2glmm::r2beta(m1_biomass,  partial = T)

Mod_results_biomass <- drop1(m1_biomass) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
 # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m1_biomass,  partial = T) %>% as.data.frame() %>% 
              rename(Driver="Effect") %>% 
              select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "biomass",.before= Driver)

Mod_results_biomass %>% 
  write_csv("results/LMM_biomass.csv")

## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m1_biomass))


##  Mowing events before sampling  ---------

mow_events_range = seq(min(Data_1m2$n_mow_events_befre_sampling),
                       max(Data_1m2$n_mow_events_befre_sampling), 
                       by = 0.001)

m1_biomass_Mow_pred <- Effect("n_mow_events_befre_sampling", m1_biomass,
                              xlevels = list(n_mow_events_befre_sampling = mow_events_range)) %>% 
  as.data.frame()

plot_Biomass_Mow_ev <- ggplot(m1_biomass_Mow_pred, 
                              aes(x = n_mow_events_befre_sampling, y = fit)) +
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

plot_Biomass_Mow_ev

##  Litter_Cover  ---------

Litter_range = seq(min(Data_1m2$Litter_Cover),
                   max(Data_1m2$Litter_Cover), 
                   by = 0.01)

m1_biomass_Litter_pred <- Effect("Litter_Cover", 
                                 m1_biomass,
                                 xlevels = list(Litter_Cover = Litter_range)) %>% 
  as.data.frame()

plot_Biomass_Litter <- ggplot(m1_biomass_Litter_pred, aes(Litter_Cover, y = fit)) +
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

plot_Biomass_Litter

##  Bare_Ground_Cover  ---------

BareSoil_range = seq(min(Data_1m2$Bare_Ground_Cover),
                     max(Data_1m2$Bare_Ground_Cover), 
                     by = 0.001)

m1_Biomass_BareSoil_pred <- Effect("Bare_Ground_Cover", 
                                   m1_biomass,
                                   xlevels = list(Bare_Ground_Cover = BareSoil_range)) %>% 
  as.data.frame() 

plot_Biomass_BareSoil <- ggplot(m1_Biomass_BareSoil_pred, 
                                aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Biomass") 

plot_Biomass_BareSoil


##  Slope effect   -------------------------------------------------------------

Slope_range = seq(min(Data_1m2$slope_degr),
                  max(Data_1m2$slope_degr), 
                  by = 0.001)

m1_Biomass_slope_pred <- Effect("slope_degr", 
                                m1_biomass,
                                xlevels = list(slope_degr = Slope_range)) %>% 
  as.data.frame() 


plot_Biomass_Slope <-ggplot(m1_Biomass_slope_pred, aes(slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope", 
       y = "Biomass") 


plot_Biomass_Slope



## Ttree distance  ----------------------------------------------------------------
Tree_dist_range = seq(min(Data_1m2$dist_tree),
                      max(Data_1m2$dist_tree), 
                      by = 0.001)

m1_Biomass_TreeDist_pred <- Effect("dist_tree", 
                                   m1_biomass,
                                   xlevels = list(dist_tree = Tree_dist_range)) %>% 
  as.data.frame() 


plot_Biomass_TreeDist <-
  ggplot(m1_Biomass_TreeDist_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Biomass") 

plot_Biomass_TreeDist


##  Sky view   -----------------------------------------------------------

Sky_range = seq(min(Data_1m2$sky_view_factor),
                max(Data_1m2$sky_view_factor), 
                by = 0.001)

m1_Biomass_Sky_pred <- Effect("sky_view_factor", 
                              m1_biomass,
                              xlevels = list(sky_view_factor = Sky_range)) %>% 
  as.data.frame() 


plot_Biomass_Sky <- ggplot(m1_Biomass_Sky_pred, 
                           aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Biomass") 


plot_Biomass_Sky

##  Biotop_richness_specific  --------------------------------------------------


Biotop_range = seq(min(Data_1m2$Biotop_richness_specific),
                   max(Data_1m2$Biotop_richness_specific), 
                   by = 0.001)

m1_biomass_Biotop_pred <- Effect("Biotop_richness_specific", 
                                 m1_biomass,
                                 xlevels = list(Biotop_richness_specific = Biotop_range)) %>% 
  as.data.frame()

plot_Biomass_Biotop <- ggplot(m1_biomass_Biotop_pred, 
                              aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Biomass") 

plot_Biomass_Biotop


##  Protected areas  ----------------------------------------------------------


ProtecArea_range = seq(min(Data_1m2$protected_cover_pct),
                       max(Data_1m2$protected_cover_pct), 
                       by = 0.001)

m1_biomass_Protect_pred <- Effect("protected_cover_pct", 
                                  m1_biomass,
                                  xlevels = list(protected_cover_pct = ProtecArea_range))%>% 
  as.data.frame() 



plot_Biomass_ProtcArea <- ggplot(m1_biomass_Protect_pred, 
                                 aes(x = protected_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(protected_cover_pct, biomass, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Protected areas (%)",
       y = "Biomass") 

plot_Biomass_ProtcArea


## Management ---------------------------------------------------------------------

m1_biomass_no_int <- update(m1_biomass, . ~ . - MowFreq:Month) # we remove interaction for lotting main effcets
drop1(m1_biomass_no_int)

emmeans(m1_biomass_no_int, list(pairwise ~ MowFreq))

emmeans_m1_biomass_MowFreq1 <- cld(emmeans(m1_biomass_no_int, list(pairwise ~ MowFreq)), 
                                   Letters = letters) %>% 
  arrange(MowFreq)

biomass_max1 <-  Data_1m2 %>% 
  summarise(max=max(biomass), .by = c(MowFreq))

plot_Biomass_Management <- Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = biomass)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Biomass") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m1_biomass_MowFreq1 %>% 
              left_join(biomass_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+400,
                # marginally signififcant difference, thus replace with custem letters
                label=  c("a", "ab", "b'" ) # .group 
            ),
            size=3.5, col="black") 

plot_Biomass_Management

## Month ----------------------------------------------------------------------

emmeans(m1_biomass_no_int, list(pairwise ~ Month))

emmeans_m1_biomass_Month <- cld(emmeans(m1_biomass_no_int, list(pairwise ~ Month)), 
                                Letters = letters) %>% 
  arrange(Month)

biomass_max2 <-  Data_1m2 %>% 
  summarise(max=max(biomass), .by = c(Month))

plot_Biomass_Month <- Data_1m2 %>% 
  ggplot(aes(x = Month, y = biomass)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Biomass") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m1_biomass_Month %>% 
              left_join(biomass_max2, by=c("Month")),
            aes(x=Month, y=max+500,
                # marginally signififcant difference, thus replace with custem letters
                label=  c("a", "ab", "b'" )
                label=.group),
            size=3.5, col="black")


plot_Biomass_Month

## MowFreq * Month  ------------------------------------------------------------

emmeans(m1_biomass, list(pairwise ~ MowFreq | Month))

emmeans_m1_biomass_MowFreq <- cld(emmeans(m1_biomass, list(pairwise ~ MowFreq | Month)), 
                                  Letters = letters) %>% 
  arrange(MowFreq) %>% 
  as.tibble()

biomass_max <-  Data_1m2 %>% 
  summarise(max=max(biomass), .by = c(MowFreq, Month))

plot_Biomass_Month_Manag_intr <- Data_1m2 %>% 
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
  geom_text(data=emmeans_m1_biomass_MowFreq %>% 
              left_join(biomass_max, by=c("MowFreq", "Month"))%>% 
              mutate(.group = case_when(
                Month == "July" & MowFreq == "regular" ~ "a",
                Month == "July" & MowFreq == "reduced" ~ "b'",
                Month == "July" & MowFreq == "reduced & sowing" ~ "ab",
                TRUE ~ .group)),
            aes(x=MowFreq, y=max+1120,
                label=.group),
            size=3.5, col="black") 

plot_Biomass_Month_Manag_intr


## COMBINE PLOTS ---------------------------------------------------------------

library(patchwork)

combined_Biomass <- 
  plot_Biomass_Mow_ev + plot_Biomass_Litter +
  plot_Biomass_BareSoil + plot_Biomass_Slope + 
  plot_Biomass_TreeDist + plot_Biomass_Sky + 
  plot_Biomass_Biotop + plot_Biomass_ProtcArea + 
  # plot_Biomass_Management + plot_Biomass_Month +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = 'bold', size = 15),
    plot.tag.position = c(0.1, 1.05),
    plot.margin = margin(t = 22, r = 10, b = 10, l = 10)
  )

print(combined_Biomass)


ggsave("results/plots/Biomass.png", combined_Biomass, width = 7, height = 10, dpi = 450)




# 2) Species diversity  ----------------------------------------------------------

# 2.1) Species richness -----------------------

## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, SR,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, impervious_pct, protected_cover_pct) %>% 
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
                 MowFreq*Month +  
                 scale(n_mow_events_befre_sampling) +
                 scale(Litter_Cover) + 
                 #      scale(road_density_km_per_ha) + 
                 #      scale(patch_size_m2) +
                 Bare_Ground_Cover + 
                 scale(slope_degr) + scale(dist_tree) + 
                 scale(sky_view_factor) +
                 scale(Biotop_richness_specific) + 
                 scale(protected_cover_pct) + #  scale(impervious_pct) +
                 #   scale(green_cover_pct) + 
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
m2_SR <- glmer(SR ~ MowFreq + Month +  
                   scale(n_mow_events_befre_sampling) +
                   scale(Litter_Cover) + 
                   #      scale(road_density_km_per_ha) + 
                   #      scale(patch_size_m2) +
                   Bare_Ground_Cover + 
                   scale(log1p(slope_degr)) + scale(dist_tree) + 
                   scale(sky_view_factor) +
                  scale(log1p(Biotop_richness_specific)) + 
                  scale(protected_cover_pct) +  # scale(impervious_pct) +
                  #   scale(green_cover_pct) + 
                   (1|PlotNo),
                 family = poisson,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)),  
                 data = Data_1m2)


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
    Bare_Ground_scaled = scale(Bare_Ground_Cover),
    slope_degr_scaled = scale(log1p(slope_degr)) ,
    dist_tree_scaled = scale(dist_tree),
    sky_view_factor_scaled = scale(sky_view_factor),
    Biotop_richness_scaled = scale(log1p(Biotop_richness_specific)),
    protected_cover_scaled = scale(protected_cover_pct))

m2_SR_dummy <- glmer(SR ~ 
                 Month + MowFreq +  
                 n_mow_events_scaled +
                 Litter_Cover_scaled + 
                 Bare_Ground_scaled + 
                 slope_degr_scaled + 
                   dist_tree_scaled + 
                 sky_view_factor_scaled +
                 Biotop_richness_scaled + 
                   protected_cover_scaled + 
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

##  Mowing events before sampling  ---------

mow_events_range = seq(min(Data_1m2$n_mow_events_befre_sampling),
                       max(Data_1m2$n_mow_events_befre_sampling), 
                       by = 0.001)

m2_SR_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_SR,
                         xlevels = list(n_mow_events_befre_sampling = mow_events_range)) %>% 
  as.data.frame()

plot_SR_Mow_ev <- ggplot(m2_SR_Mow_pred, 
                         aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Species richness") 

plot_SR_Mow_ev

##  Litter_Cover  ---------

Litter_range = seq(min(Data_1m2$Litter_Cover),
                   max(Data_1m2$Litter_Cover), 
                   by = 0.01)

m2_SR_Litter_pred <- Effect("Litter_Cover", 
                            m2_SR,
                            xlevels = list(Litter_Cover = Litter_range)) %>% 
  as.data.frame()

plot_SR_Litter <- ggplot(m2_SR_Litter_pred, aes(Litter_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Litter cover, %", 
       y = "Species richness") 

plot_SR_Litter

##  Bare_Ground_Cover  ---------

BareSoil_range = seq(min(Data_1m2$Bare_Ground_Cover),
                     max(Data_1m2$Bare_Ground_Cover), 
                     by = 0.001)

m2_SR_BareSoil_pred <- Effect("Bare_Ground_Cover", 
                              m2_SR,
                              xlevels = list(Bare_Ground_Cover = BareSoil_range)) %>% 
  as.data.frame() 

plot_SR_BareSoil <- ggplot(m2_SR_BareSoil_pred, 
                           aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Species richness") 

plot_SR_BareSoil


##  Slope effect   -------------------------------------------------------------

Slope_range = seq(min(Data_1m2$slope_degr),
                  max(Data_1m2$slope_degr), 
                  by = 0.001)

m2_SR_slope_pred <- Effect("slope_degr", 
                           m2_SR,
                           xlevels = list(slope_degr = Slope_range)) %>% 
  as.data.frame() 


plot_SR_Slope <-ggplot(m2_SR_slope_pred, aes(slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope", 
       y = "Species richness") 


plot_SR_Slope



## Ttree distance  ----------------------------------------------------------------
Tree_dist_range = seq(min(Data_1m2$dist_tree),
                      max(Data_1m2$dist_tree), 
                      by = 0.001)

m2_SR_TreeDist_pred <- Effect("dist_tree", 
                              m2_SR,
                              xlevels = list(dist_tree = Tree_dist_range)) %>% 
  as.data.frame() 


plot_SR_TreeDist <-
  ggplot(m2_SR_TreeDist_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Species richness") 

plot_SR_TreeDist


##  Sky view   -----------------------------------------------------------

Sky_range = seq(min(Data_1m2$sky_view_factor),
                max(Data_1m2$sky_view_factor), 
                by = 0.001)

m2_SR_Sky_pred <- Effect("sky_view_factor", 
                         m2_SR,
                         xlevels = list(sky_view_factor = Sky_range)) %>% 
  as.data.frame() 


plot_SR_Sky <- ggplot(m2_SR_Sky_pred, 
                      aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Species richness") 


plot_SR_Sky

##  Biotop_richness_specific  --------------------------------------------------


Biotop_range = seq(min(Data_1m2$Biotop_richness_specific),
                   max(Data_1m2$Biotop_richness_specific), 
                   by = 0.001)

m2_SR_Biotop_pred <- Effect("Biotop_richness_specific", 
                            m2_SR,
                            xlevels = list(Biotop_richness_specific = Biotop_range)) %>% 
  as.data.frame()

plot_SR_Biotop <- ggplot(m2_SR_Biotop_pred, 
                         aes(Biotop_richness_specific, y = fit)) +
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

plot_SR_Biotop


##  Protected areas  ----------------------------------------------------------


ProtecArea_range = seq(min(Data_1m2$protected_cover_pct),
                       max(Data_1m2$protected_cover_pct), 
                       by = 0.001)

m2_SR_Protect_pred <- Effect("protected_cover_pct", 
                             m2_SR,
                             xlevels = list(protected_cover_pct = ProtecArea_range))%>% 
  as.data.frame() 



plot_SR_ProtcArea <- ggplot(m2_SR_Protect_pred, 
                            aes(x = protected_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(protected_cover_pct, SR, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Protected areas (%)",
       y = "Species richness") 

plot_SR_ProtcArea


## Management ---------------------------------------------------------------------

emmeans(m2_SR, list(pairwise ~ MowFreq))

emmeans_m2_SR_MowFreq1 <- cld(emmeans(m2_SR, list(pairwise ~ MowFreq)), 
                              Letters = letters) %>% 
  arrange(MowFreq)

SR_max1 <-  Data_1m2 %>% 
  summarise(max=max(SR), .by = c(MowFreq))

plot_SR_Management <- Data_1m2 %>% 
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
            aes(x=MowFreq, y=max+2, label= .group 
            ),
            size=3.5, col="black") 

plot_SR_Management

## Month ----------------------------------------------------------------------

emmeans(m2_SR, list(pairwise ~ Month))

emmeans_m2_SR_Month <- cld(emmeans(m2_SR, list(pairwise ~ Month)), 
                           Letters = letters) %>% 
  arrange(Month)

SR_max2 <-  Data_1m2 %>% 
  summarise(max=max(SR), .by = c(Month))

plot_SR_Month <- Data_1m2 %>% 
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
            aes(x=Month, y=max+2,
                label=.group),
            size=3.5, col="black")


plot_SR_Month

## MowFreq * Month  ------------------------------------------------------------

# use model with interactionns

emmeans(m1_SR, list(pairwise ~ MowFreq | Month))

emmeans_m2_SR_MowFreq <- cld(emmeans(m1_SR, list(pairwise ~ MowFreq | Month)), 
                             Letters = letters) %>% 
  arrange(MowFreq) %>% 
  as.tibble()

SR_max <-  Data_1m2 %>% 
  summarise(max=max(SR), .by = c(MowFreq, Month))

plot_SR_Month_Manag_intr <- Data_1m2 %>% 
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
              left_join(SR_max, by=c("MowFreq", "Month"))%>% 
              mutate(.group = case_when(
                Month == "July" & MowFreq == "regular" ~ "a",
                Month == "July" & MowFreq == "reduced" ~ "a",
                Month == "July" & MowFreq == "reduced & sowing" ~ "b'",
                TRUE ~ .group)),
            aes(x=MowFreq, y=max+2,
                label=.group),
            size=3.5, col="black") 

plot_SR_Month_Manag_intr


## COMBINE PLOTS ---------------------------------------------------------------

library(patchwork)

combined_SR <- 
  plot_SR_Mow_ev + plot_SR_Litter +
  plot_SR_BareSoil + plot_SR_Slope + 
  plot_SR_TreeDist + plot_SR_Sky + 
  plot_SR_Biotop + plot_SR_ProtcArea + 
  # plot_SR_Management + plot_SR_Month +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = 'bold', size = 15),
    plot.tag.position = c(0.1, 1.05),
    plot.margin = margin(t = 22, r = 10, b = 10, l = 10)
  )

print(combined_SR)


ggsave("results/plots/SR.png", combined_SR, width = 7, height = 10, dpi = 450)





# 2.2) evenness -----------------------------------------------------------------------


## Exploration: ----------------------------

Data_1m2 %>%
  select(Month, evenness,
         n_mow_events_befre_sampling,
         Bare_Ground_Cover, Litter_Cover, slope_degr, dist_tree, 
         sky_view_factor, patch_size_m2, Biotop_richness_specific, 
         green_cover_pct_log, impervious_pct, protected_cover_pct) %>% 
  pivot_longer(-c(Month, evenness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = evenness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Evenness")


## Models: ----------------------------
Data_1m2 %>% 
  select(evenness)

m1_evenness <- lmerTest::lmer(evenness ~ 
                                MowFreq*Month +  
                                scale(n_mow_events_befre_sampling) +
                                scale(Litter_Cover) + 
                                #      scale(road_density_km_per_ha) + 
                                #      scale(patch_size_m2) +
                                Bare_Ground_Cover + 
                                scale(slope_degr) + scale(dist_tree) + 
                                scale(sky_view_factor) +
                                scale(Biotop_richness_specific) + 
                                scale(protected_cover_pct) + 
                                #  scale(green_cover_pct) + 
                                scale(impervious_pct) +
                                (1|PlotNo),
                              data = Data_1m2)


# check model assumptions
check_convergence(m1_evenness)
check_model(m1_evenness)
check_collinearity(m1_evenness)
# check interactions
drop1(m1_evenness)

# remove interaction 

m2_evenness <- lmerTest::lmer(evenness ~ 
                                  MowFreq + Month +  
                                  scale(n_mow_events_befre_sampling) +
                                  scale(Litter_Cover) + 
                                  #      scale(road_density_km_per_ha) + 
                                  #      scale(patch_size_m2) +
                                  Bare_Ground_Cover + 
                                  scale(slope_degr) + scale(dist_tree) + 
                                  scale(sky_view_factor) +
                                  scale(Biotop_richness_specific) + 
                                protected_cover_pct +  # scale(impervious_pct) +
                                #   scale(green_cover_pct) + 
                                  (1|PlotNo),
                                data = Data_1m2)


summary(m2_evenness)
# check model assumptions
check_convergence(m2_evenness)
check_model(m2_evenness)
check_collinearity(m2_evenness)

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


##  Mowing events before sampling  ---------

mow_events_range = seq(min(Data_1m2$n_mow_events_befre_sampling),
                       max(Data_1m2$n_mow_events_befre_sampling), 
                       by = 0.001)

m2_evenness_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_evenness,
                               xlevels = list(n_mow_events_befre_sampling = mow_events_range)) %>% 
  as.data.frame()

plot_evenness_Mow_ev <- ggplot(m2_evenness_Mow_pred, 
                               aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Hill-Simpson diversity") 

plot_evenness_Mow_ev

##  Litter_Cover  ---------

Litter_range = seq(min(Data_1m2$Litter_Cover),
                   max(Data_1m2$Litter_Cover), 
                   by = 0.01)

m2_evenness_Litter_pred <- Effect("Litter_Cover", 
                                  m2_evenness,
                                  xlevels = list(Litter_Cover = Litter_range)) %>% 
  as.data.frame()

plot_evenness_Litter <- ggplot(m2_evenness_Litter_pred, aes(Litter_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Litter cover, %", 
       y = "Hill-Simpson diversity") 

plot_evenness_Litter

##  Bare_Ground_Cover  ---------

BareSoil_range = seq(min(Data_1m2$Bare_Ground_Cover),
                     max(Data_1m2$Bare_Ground_Cover), 
                     by = 0.001)

m2_evenness_BareSoil_pred <- Effect("Bare_Ground_Cover", 
                                    m2_evenness,
                                    xlevels = list(Bare_Ground_Cover = BareSoil_range)) %>% 
  as.data.frame() 

plot_evenness_BareSoil <- ggplot(m2_evenness_BareSoil_pred, 
                                 aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Hill-Simpson diversity") 

plot_evenness_BareSoil


##  Slope effect   -------------------------------------------------------------

Slope_range = seq(min(Data_1m2$slope_degr),
                  max(Data_1m2$slope_degr), 
                  by = 0.001)

m2_evenness_slope_pred <- Effect("slope_degr", 
                                 m2_evenness,
                                 xlevels = list(slope_degr = Slope_range)) %>% 
  as.data.frame() 


plot_evenness_Slope <-ggplot(m2_evenness_slope_pred, aes(slope_degr, y = fit)) +
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


plot_evenness_Slope



## Ttree distance  ----------------------------------------------------------------
Tree_dist_range = seq(min(Data_1m2$dist_tree),
                      max(Data_1m2$dist_tree), 
                      by = 0.001)

m2_evenness_TreeDist_pred <- Effect("dist_tree", 
                                    m2_evenness,
                                    xlevels = list(dist_tree = Tree_dist_range)) %>% 
  as.data.frame() 


plot_evenness_TreeDist <-
  ggplot(m2_evenness_TreeDist_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Hill-Simpson diversity") 

plot_evenness_TreeDist


##  Sky view   -----------------------------------------------------------

Sky_range = seq(min(Data_1m2$sky_view_factor),
                max(Data_1m2$sky_view_factor), 
                by = 0.001)

m2_evenness_Sky_pred <- Effect("sky_view_factor", 
                               m2_evenness,
                               xlevels = list(sky_view_factor = Sky_range)) %>% 
  as.data.frame() 


plot_evenness_Sky <- ggplot(m2_evenness_Sky_pred, 
                            aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Hill-Simpson diversity") 


plot_evenness_Sky

##  Biotop_richness_specific  --------------------------------------------------


Biotop_range = seq(min(Data_1m2$Biotop_richness_specific),
                   max(Data_1m2$Biotop_richness_specific), 
                   by = 0.001)

m2_evenness_Biotop_pred <- Effect("Biotop_richness_specific", 
                                  m2_evenness,
                                  xlevels = list(Biotop_richness_specific = Biotop_range)) %>% 
  as.data.frame()

plot_evenness_Biotop <- ggplot(m2_evenness_Biotop_pred, 
                               aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Hill-Simpson diversity") 

plot_evenness_Biotop


##  Protected areas  ----------------------------------------------------------


ProtecArea_range = seq(min(Data_1m2$protected_cover_pct),
                       max(Data_1m2$protected_cover_pct), 
                       by = 0.001)

m2_evenness_Protect_pred <- Effect("protected_cover_pct", 
                                   m2_evenness,
                                   xlevels = list(protected_cover_pct = ProtecArea_range))%>% 
  as.data.frame() 



plot_evenness_ProtcArea <- ggplot(m2_evenness_Protect_pred, 
                                  aes(x = protected_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(protected_cover_pct, evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Protected areas (%)",
       y = "Hill-Simpson diversity") 

plot_evenness_ProtcArea


## Management ---------------------------------------------------------------------

emmeans(m2_evenness, list(pairwise ~ MowFreq))

emmeans_m2_evenness_MowFreq1 <- cld(emmeans(m2_evenness, list(pairwise ~ MowFreq)), 
                                    Letters = letters) %>% 
  arrange(MowFreq)

evenness_max1 <-  Data_1m2 %>% 
  summarise(max=max(evenness), .by = c(MowFreq))

plot_evenness_Management <- Data_1m2 %>% 
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
            aes(x=MowFreq, y=max+2, label= .group 
            ),
            size=3.5, col="black") 

plot_evenness_Management

## Month ----------------------------------------------------------------------

emmeans(m2_evenness, list(pairwise ~ Month))

emmeans_m2_evenness_Month <- cld(emmeans(m2_evenness, list(pairwise ~ Month)), 
                                 Letters = letters) %>% 
  arrange(Month)

evenness_max2 <-  Data_1m2 %>% 
  summarise(max=max(evenness), .by = c(Month))

plot_evenness_Month <- Data_1m2 %>% 
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
            aes(x=Month, y=max+2,
                label=.group),
            size=3.5, col="black")


plot_evenness_Month

## MowFreq * Month  ------------------------------------------------------------

# use model with interactionns

emmeans(m1_evenness, list(pairwise ~ MowFreq | Month))

emmeans_m2_evenness_MowFreq <- cld(emmeans(m1_evenness, list(pairwise ~ MowFreq | Month)), 
                                   Letters = letters) %>% 
  arrange(MowFreq) %>% 
  as.tibble()

evenness_max <-  Data_1m2 %>% 
  summarise(max=max(evenness), .by = c(MowFreq, Month))

plot_evenness_Month_Manag_intr <- Data_1m2 %>% 
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
            aes(x=MowFreq, y=max+2,
                label=.group),
            size=3.5, col="black") 

plot_evenness_Month_Manag_intr


## COMBINE PLOTS ---------------------------------------------------------------

library(patchwork)

combined_evenness <- 
  plot_evenness_Mow_ev + plot_evenness_Litter +
  plot_evenness_BareSoil + plot_evenness_Slope + 
  plot_evenness_TreeDist + plot_evenness_Sky + 
  plot_evenness_Biotop + plot_evenness_ProtcArea + 
  # plot_evenness_Management + plot_evenness_Month +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = 'bold', size = 15),
    plot.tag.position = c(0.1, 1.05),
    plot.margin = margin(t = 22, r = 10, b = 10, l = 10)
  )

print(combined_SR)


ggsave("results/plots/evenness.png", combined_evenness, width = 7, height = 10, dpi = 450)


# 3) Phenological Diversity -----------------------------------------------------------------------


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
  labs(x = NULL, y = "phen_Richness")


## Models: ----------------------------
Data_1m2 %>% 
  select(phen_Richness)

m1_phen_Richness <- lmerTest::lmer(phen_Richness ~ 
                                     MowFreq*Month +  
                                     scale(n_mow_events_befre_sampling) +
                                     scale(Litter_Cover) + 
                                     #      scale(road_density_km_per_ha) + 
                                     #      scale(patch_size_m2) +
                                     Bare_Ground_Cover + 
                                     scale(log1p(slope_degr)) + 
                                     scale(dist_tree) + 
                                     scale(sky_view_factor) +
                                     scale(Biotop_richness_specific) + 
                                     scale(protected_cover_pct) + 
                                    # scale(impervious_pct) +
                                     #   scale(green_cover_pct) + 
                                     (1|PlotNo),
                                   data = Data_1m2)


# check model assumptions
check_convergence(m1_phen_Richness)
check_model(m1_phen_Richness)
check_collinearity(m1_phen_Richness)
# check interactions
drop1(m1_phen_Richness)

# check predictor effects
drop1(m1_phen_Richness)
# Anova(m1_phen_Richness, type = "III", test.statistic = "F")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m1_phen_Richness)
# Partial R2 for fixed effects
r2glmm::r2beta(m1_phen_Richness,  partial = T)

Mod_results_Phen_Richness <- drop1(m1_phen_Richness) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m1_phen_Richness,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "phen_Richness",.before= Driver)

Mod_results_Phen_Richness %>% 
  write_csv("results/LMM_Phenol_Richness.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m1_phen_Richness))

##  Mowing events before sampling  ---------

mow_events_range = seq(min(Data_1m2$n_mow_events_befre_sampling),
                       max(Data_1m2$n_mow_events_befre_sampling), 
                       by = 0.001)

m1_phen_Richness_Mow_pred <- Effect("n_mow_events_befre_sampling", m1_phen_Richness,
                                    xlevels = list(n_mow_events_befre_sampling = mow_events_range)) %>% 
  as.data.frame()

plot_PhenRichness_Mow_ev <- ggplot(m1_phen_Richness_Mow_pred, 
                                   aes(x = n_mow_events_befre_sampling, y = fit)) +
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

plot_PhenRichness_Mow_ev

##  Litter_Cover  ---------

Litter_range = seq(min(Data_1m2$Litter_Cover),
                   max(Data_1m2$Litter_Cover), 
                   by = 0.01)

m1_phen_Richness_Litter_pred <- Effect("Litter_Cover", 
                                       m1_phen_Richness,
                                       xlevels = list(Litter_Cover = Litter_range)) %>% 
  as.data.frame()

plot_PhenRichness_Litter <- ggplot(m1_phen_Richness_Litter_pred, aes(Litter_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, phen_Richness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0.15))  +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Litter cover, %", 
       y = "Phenological richness") 

plot_PhenRichness_Litter

##  Bare_Ground_Cover  ---------

BareSoil_range = seq(min(Data_1m2$Bare_Ground_Cover),
                     max(Data_1m2$Bare_Ground_Cover), 
                     by = 0.001)

m1_phen_Richness_BareSoil_pred <- Effect("Bare_Ground_Cover", 
                                         m1_phen_Richness,
                                         xlevels = list(Bare_Ground_Cover = BareSoil_range)) %>% 
  as.data.frame() 

plot_PhenRichness_BareSoil <- ggplot(m1_phen_Richness_BareSoil_pred, 
                                     aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, phen_Richness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Phenological richness") 

plot_PhenRichness_BareSoil


##  Slope effect   -------------------------------------------------------------

Slope_range = seq(min(Data_1m2$slope_degr),
                  max(Data_1m2$slope_degr), 
                  by = 0.001)

m1_phen_Richness_slope_pred <- Effect("slope_degr", 
                                      m1_phen_Richness,
                                      xlevels = list(slope_degr = Slope_range)) %>% 
  as.data.frame() 


plot_PhenRichness_Slope <-ggplot(m1_phen_Richness_slope_pred, aes(slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, phen_Richness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope", 
       y = "Phenological richness") 


plot_PhenRichness_Slope



## Ttree distance  ----------------------------------------------------------------
Tree_dist_range = seq(min(Data_1m2$dist_tree),
                      max(Data_1m2$dist_tree), 
                      by = 0.001)

m1_phen_Richness_TreeDist_pred <- Effect("dist_tree", 
                                         m1_phen_Richness,
                                         xlevels = list(dist_tree = Tree_dist_range)) %>% 
  as.data.frame() 


plot_PhenRichness_TreeDist <-
  ggplot(m1_phen_Richness_TreeDist_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, phen_Richness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Phenological richness") 

plot_PhenRichness_TreeDist


##  Sky view   -----------------------------------------------------------

Sky_range = seq(min(Data_1m2$sky_view_factor),
                max(Data_1m2$sky_view_factor), 
                by = 0.001)

m1_phen_Richness_Sky_pred <- Effect("sky_view_factor", 
                                    m1_phen_Richness,
                                    xlevels = list(sky_view_factor = Sky_range)) %>% 
  as.data.frame() 


plot_PhenRichness_Sky <- ggplot(m1_phen_Richness_Sky_pred, 
                                aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, phen_Richness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Phenological richness") 


plot_PhenRichness_Sky

##  Biotop_richness_specific  --------------------------------------------------


Biotop_range = seq(min(Data_1m2$Biotop_richness_specific),
                   max(Data_1m2$Biotop_richness_specific), 
                   by = 0.001)

m1_phen_Richness_Biotop_pred <- Effect("Biotop_richness_specific", 
                                       m1_phen_Richness,
                                       xlevels = list(Biotop_richness_specific = Biotop_range)) %>% 
  as.data.frame()

plot_PhenRichness_Biotop <- ggplot(m1_phen_Richness_Biotop_pred, 
                                   aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, phen_Richness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Phenological richness") 

plot_PhenRichness_Biotop


##  Protected areas  ----------------------------------------------------------


ProtecArea_range = seq(min(Data_1m2$protected_cover_pct),
                       max(Data_1m2$protected_cover_pct), 
                       by = 0.001)

m1_phen_Richness_Protect_pred <- Effect("protected_cover_pct", 
                                        m1_phen_Richness,
                                        xlevels = list(protected_cover_pct = ProtecArea_range))%>% 
  as.data.frame() 



plot_PhenRichness_ProtcArea <- ggplot(m1_phen_Richness_Protect_pred, 
                                      aes(x = protected_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(protected_cover_pct, phen_Richness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Protected areas (%)",
       y = "Phenological richness") 

plot_PhenRichness_ProtcArea


## Management ---------------------------------------------------------------------

m1_phen_Richness_no_int <- update(m1_phen_Richness, . ~ . - MowFreq:Month) # we remove interaction for lotting main effcets
drop1(m1_phen_Richness_no_int)

emmeans(m1_phen_Richness_no_int, list(pairwise ~ MowFreq))

emmeans_m1_phen_Richness_MowFreq1 <- cld(emmeans(m1_phen_Richness_no_int, list(pairwise ~ MowFreq)), 
                                         Letters = letters) %>% 
  arrange(MowFreq)

phen_Richness_max1 <-  Data_1m2 %>% 
  summarise(max=max(phen_Richness), .by = c(MowFreq))

plot_PhenRichness_Management <- Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = phen_Richness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0.1)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Phenological richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m1_phen_Richness_MowFreq1 %>% 
              left_join(phen_Richness_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+0.5,
                label= .group 
            ),
            size=3.5, col="black") 

plot_PhenRichness_Management

## Month ----------------------------------------------------------------------

emmeans(m1_phen_Richness_no_int, list(pairwise ~ Month))

emmeans_m1_phen_Richness_Month <- cld(emmeans(m1_phen_Richness_no_int, list(pairwise ~ Month)), 
                                      Letters = letters) %>% 
  arrange(Month)

phen_Richness_max2 <-  Data_1m2 %>% 
  summarise(max=max(phen_Richness), .by = c(Month))

plot_PhenRichness_Month <- Data_1m2 %>% 
  ggplot(aes(x = Month, y = phen_Richness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0.1)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Phenological richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m1_phen_Richness_Month %>% 
              left_join(phen_Richness_max2, by=c("Month")),
            aes(x=Month, y=max+1,
                # marginally signififcant difference, thus replace with custem letters
                #  label=  c("a", "ab", "b'" ),
                label=.group),
            size=3.5, col="black")


plot_PhenRichness_Month

## MowFreq * Month  ------------------------------------------------------------

emmeans(m1_phen_Richness, list(pairwise ~ MowFreq | Month))

emmeans_m1_phen_Richness_MowFreq <- cld(emmeans(m1_phen_Richness, list(pairwise ~ MowFreq | Month)), 
                                        Letters = letters) %>% 
  arrange(MowFreq) %>% 
  as.tibble()

phen_Richness_max <-  Data_1m2 %>% 
  summarise(max=max(phen_Richness), .by = c(MowFreq, Month))

plot_PhenRichness_Month_Manag_intr <- Data_1m2 %>% 
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
  geom_text(data=emmeans_m1_phen_Richness_MowFreq %>% 
              left_join(phen_Richness_max, by=c("MowFreq", "Month"))%>% 
              mutate(.group = case_when(
                Month == "May" & MowFreq == "regular" ~ "a",
                Month == "May" & MowFreq == "reduced" ~ "a",
                Month == "May" & MowFreq == "reduced & sowing" ~ "b",
                TRUE ~ .group)),
            aes(x=MowFreq, y=max+1,
                label=.group),
            size=3.5, col="black") 

plot_PhenRichness_Month_Manag_intr


## COMBINE PLOTS ---------------------------------------------------------------

library(patchwork)

combined_phen_Richness <- 
  plot_PhenRichness_Mow_ev + plot_PhenRichness_Litter +
  plot_PhenRichness_BareSoil + plot_PhenRichness_Slope + 
  plot_PhenRichness_TreeDist + plot_PhenRichness_Sky + 
  plot_PhenRichness_Biotop + plot_PhenRichness_ProtcArea + 
  # plot_PhenRichness_Management + plot_PhenRichness_Month +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = 'bold', size = 15),
    plot.tag.position = c(0.1, 1.05),
    plot.margin = margin(t = 22, r = 10, b = 10, l = 10)
  )

print(combined_phen_Richness)


ggsave("results/plots/phen_Richness.png", combined_phen_Richness, width = 7, height = 10, dpi = 450)



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
                                     #      scale(road_density_km_per_ha) + 
                                     #      scale(patch_size_m2) +
                                     Bare_Ground_Cover + 
                                     scale(slope_degr) + scale(dist_tree) + 
                                     scale(sky_view_factor) +
                                     scale(log1p(Biotop_richness_specific)) + 
                                     scale(protected_cover_pct) + 
                                    # scale(impervious_pct) +
                                     #   scale(green_cover_pct) + 
                                     (1|PlotNo),
                                   data = Data_1m2)


# check model assumptions
check_convergence(m1_phen_evenness)
check_model(m1_phen_evenness)
check_collinearity(m1_phen_evenness)
# check interactions
drop1(m1_phen_evenness)

# interaction is  significant

# check predictor effects
drop1(m1_phen_evenness)
# Anova(m1_phen_evenness, type = "III", test.statistic = "F")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m1_phen_evenness)
# Partial R2 for fixed effects
r2glmm::r2beta(m1_phen_evenness,  partial = T)

Mod_results_Phen_evenness <- drop1(m1_phen_evenness) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m1_phen_evenness,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "phen_evenness",.before= Driver)

Mod_results_Phen_evenness %>% 
  write_csv("results/LMM_Phen_evenness.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m1_phen_evenness))

##  Mowing events before sampling  ---------

mow_events_range = seq(min(Data_1m2$n_mow_events_befre_sampling),
                       max(Data_1m2$n_mow_events_befre_sampling), 
                       by = 0.001)

m1_phen_evenness_Mow_pred <- Effect("n_mow_events_befre_sampling", m1_phen_evenness,
                                    xlevels = list(n_mow_events_befre_sampling = mow_events_range)) %>% 
  as.data.frame()

plot_PhenEvenness_Mow_ev <- ggplot(m1_phen_evenness_Mow_pred, 
                                   aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.2, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Phenological evenness") 

plot_PhenEvenness_Mow_ev

##  Litter_Cover  ---------

Litter_range = seq(min(Data_1m2$Litter_Cover),
                   max(Data_1m2$Litter_Cover), 
                   by = 0.01)

m1_phen_evenness_Litter_pred <- Effect("Litter_Cover", 
                                       m1_phen_evenness,
                                       xlevels = list(Litter_Cover = Litter_range)) %>% 
  as.data.frame()

plot_PhenEvenness_Litter <- ggplot(m1_phen_evenness_Litter_pred, aes(Litter_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0.15))  +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Litter cover, %", 
       y = "Phenological evenness") 

plot_PhenEvenness_Litter

##  Bare_Ground_Cover  ---------

BareSoil_range = seq(min(Data_1m2$Bare_Ground_Cover),
                     max(Data_1m2$Bare_Ground_Cover), 
                     by = 0.001)

m1_phen_evenness_BareSoil_pred <- Effect("Bare_Ground_Cover", 
                                         m1_phen_evenness,
                                         xlevels = list(Bare_Ground_Cover = BareSoil_range)) %>% 
  as.data.frame() 

plot_PhenEvenness_BareSoil <- ggplot(m1_phen_evenness_BareSoil_pred, 
                                     aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Phenological evenness") 

plot_PhenEvenness_BareSoil


##  Slope effect   -------------------------------------------------------------

Slope_range = seq(min(Data_1m2$slope_degr),
                  max(Data_1m2$slope_degr), 
                  by = 0.001)

m1_phen_evenness_slope_pred <- Effect("slope_degr", 
                                      m1_phen_evenness,
                                      xlevels = list(slope_degr = Slope_range)) %>% 
  as.data.frame() 


plot_PhenEvenness_Slope <-ggplot(m1_phen_evenness_slope_pred, aes(slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope", 
       y = "Phenological evenness") 


plot_PhenEvenness_Slope



## Ttree distance  ----------------------------------------------------------------
Tree_dist_range = seq(min(Data_1m2$dist_tree),
                      max(Data_1m2$dist_tree), 
                      by = 0.001)

m1_phen_evenness_TreeDist_pred <- Effect("dist_tree", 
                                         m1_phen_evenness,
                                         xlevels = list(dist_tree = Tree_dist_range)) %>% 
  as.data.frame() 


plot_PhenEvenness_TreeDist <-
  ggplot(m1_phen_evenness_TreeDist_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Phenological evenness") 

plot_PhenEvenness_TreeDist


##  Sky view   -----------------------------------------------------------

Sky_range = seq(min(Data_1m2$sky_view_factor),
                max(Data_1m2$sky_view_factor), 
                by = 0.001)

m1_phen_evenness_Sky_pred <- Effect("sky_view_factor", 
                                    m1_phen_evenness,
                                    xlevels = list(sky_view_factor = Sky_range)) %>% 
  as.data.frame() 


plot_PhenEvenness_Sky <- ggplot(m1_phen_evenness_Sky_pred, 
                                aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Phenological evenness") 


plot_PhenEvenness_Sky

##  Biotop_richness_specific  --------------------------------------------------


Biotop_range = seq(min(Data_1m2$Biotop_richness_specific),
                   max(Data_1m2$Biotop_richness_specific), 
                   by = 0.001)

m1_phen_evenness_Biotop_pred <- Effect("Biotop_richness_specific", 
                                       m1_phen_evenness,
                                       xlevels = list(Biotop_richness_specific = Biotop_range)) %>% 
  as.data.frame()

plot_PhenEvenness_Biotop <- ggplot(m1_phen_evenness_Biotop_pred, 
                                   aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.2, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Phenological evenness") 

plot_PhenEvenness_Biotop


##  Protected areas  ----------------------------------------------------------


ProtecArea_range = seq(min(Data_1m2$protected_cover_pct),
                       max(Data_1m2$protected_cover_pct), 
                       by = 0.001)

m1_phen_evenness_Protect_pred <- Effect("protected_cover_pct", 
                                        m1_phen_evenness,
                                        xlevels = list(protected_cover_pct = ProtecArea_range))%>% 
  as.data.frame() 



plot_PhenEvenness_ProtcArea <- ggplot(m1_phen_evenness_Protect_pred, 
                                      aes(x = protected_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(protected_cover_pct, phen_evenness, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.1, height = 0.04)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Protected areas (%)",
       y = "Phenological evenness") 

plot_PhenEvenness_ProtcArea


## Management ---------------------------------------------------------------------

m1_phen_evenness_no_int <- update(m1_phen_evenness, . ~ . - MowFreq:Month) # we remove interaction for lotting main effcets
drop1(m1_phen_evenness_no_int)

emmeans(m1_phen_evenness_no_int, list(pairwise ~ MowFreq))

emmeans_m1_phen_evenness_MowFreq1 <- cld(emmeans(m1_phen_evenness_no_int, list(pairwise ~ MowFreq)), 
                                         Letters = letters) %>% 
  arrange(MowFreq)

phen_evenness_max1 <-  Data_1m2 %>% 
  summarise(max=max(phen_evenness), .by = c(MowFreq))

plot_PhenEvenness_Management <- Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = phen_evenness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0.1)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Phenological evenness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m1_phen_evenness_MowFreq1 %>% 
              left_join(phen_evenness_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+0.5,
                label= .group 
            ),
            size=3.5, col="black") 

plot_PhenEvenness_Management

## Month ----------------------------------------------------------------------

emmeans(m1_phen_evenness_no_int, list(pairwise ~ Month))

emmeans_m1_phen_evenness_Month <- cld(emmeans(m1_phen_evenness_no_int, list(pairwise ~ Month)), 
                                      Letters = letters) %>% 
  arrange(Month)

phen_evenness_max2 <-  Data_1m2 %>% 
  summarise(max=max(phen_evenness), .by = c(Month))

plot_PhenEvenness_Month <- Data_1m2 %>% 
  ggplot(aes(x = Month, y = phen_evenness)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0.1)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Phenological evenness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m1_phen_evenness_Month %>% 
              left_join(phen_evenness_max2, by=c("Month")),
            aes(x=Month, y=max+1,
                # marginally signififcant difference, thus replace with custem letters
                #  label=  c("a", "ab", "b'" ),
                label=.group),
            size=3.5, col="black")


plot_PhenEvenness_Month

## MowFreq * Month  ------------------------------------------------------------

emmeans(m1_phen_evenness, list(pairwise ~ MowFreq | Month))

emmeans_m1_phen_evenness_MowFreq <- cld(emmeans(m1_phen_evenness, list(pairwise ~ MowFreq | Month)), 
                                        Letters = letters) %>% 
  arrange(MowFreq) %>% 
  as.tibble()

phen_evenness_max <-  Data_1m2 %>% 
  summarise(max=max(phen_evenness), .by = c(MowFreq, Month))

plot_PhenEvenness_Month_Manag_intr <- Data_1m2 %>% 
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
  geom_text(data=emmeans_m1_phen_evenness_MowFreq %>% 
              left_join(phen_evenness_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+1,
                label=.group),
            size=3.5, col="black") 

plot_PhenEvenness_Month_Manag_intr


## COMBINE PLOTS ---------------------------------------------------------------

library(patchwork)

combined_phen_evenness <- 
  plot_PhenEvenness_Mow_ev + plot_PhenEvenness_Litter +
  plot_PhenEvenness_BareSoil + plot_PhenEvenness_Slope + 
  plot_PhenEvenness_TreeDist + plot_PhenEvenness_Sky + 
  plot_PhenEvenness_Biotop + plot_PhenEvenness_ProtcArea + 
  # plot_PhenEvenness_Management + plot_PhenEvenness_Month +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = 'bold', size = 15),
    plot.tag.position = c(0.1, 1.05),
    plot.margin = margin(t = 22, r = 10, b = 10, l = 10)
  )

print(combined_phen_evenness)


ggsave("results/plots/phen_evenness.png", combined_phen_evenness, width = 7, height = 10, dpi = 450)




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
                            #      scale(road_density_km_per_ha) + 
                            #      scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(log1p(Biotop_richness_specific)) +
                              scale(protected_cover_pct) + 
                           # scale(green_cover_pct) + 
                            (1|PlotNo),
                          data = Data_1m2)


# check model assumptions
check_convergence(m1_FRic)
check_model(m1_FRic)
check_collinearity(m1_FRic)
# check interactions
drop1(m1_FRic)

# interaction is  marginally significant


# check predictor effects
drop1(m1_FRic)
# Anova(m1_FRic, type = "II")

## R2 ---------------------------------------------------------------
# R2 for the entire model
MuMIn::r.squaredGLMM(m1_FRic)
# Partial R2 for fixed effects
r2glmm::r2beta(m1_FRic,  partial = T)

Mod_results_FRic <- drop1(m1_FRic) %>% as.data.frame() %>% 
  rownames_to_column("Driver") %>% select(-"Sum Sq", -"Mean Sq") %>% 
  # relocate raw "MowFreq:Month" in column "Driver" to the top
  arrange(Driver != "MowFreq:Month") %>% 
  left_join(
    r2glmm::r2beta(m1_FRic,  partial = T) %>% as.data.frame() %>% 
      rename(Driver="Effect") %>% 
      select(Driver,  Rsq), by = "Driver") %>% 
  mutate(Responce = "FRic",.before= Driver)

Mod_results_FRic %>% 
  write_csv("results/LMM_Functional_richness.csv")
## Plots ------------------------------------------------------------------------

library(effects)
plot(allEffects(m1_FRic))

##  Mowing events before sampling  ---------

mow_events_range = seq(min(Data_1m2$n_mow_events_befre_sampling),
                       max(Data_1m2$n_mow_events_befre_sampling), 
                       by = 0.001)

m1_FRic_Mow_pred <- Effect("n_mow_events_befre_sampling", m1_FRic,
                           xlevels = list(n_mow_events_befre_sampling = mow_events_range)) %>% 
  as.data.frame()

plot_FRic_Mow_ev <- ggplot(m1_FRic_Mow_pred, 
                           aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.2, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Functional richness") 

plot_FRic_Mow_ev

##  Litter_Cover  ---------

Litter_range = seq(min(Data_1m2$Litter_Cover),
                   max(Data_1m2$Litter_Cover), 
                   by = 0.01)

m1_FRic_Litter_pred <- Effect("Litter_Cover", 
                              m1_FRic,
                              xlevels = list(Litter_Cover = Litter_range)) %>% 
  as.data.frame()

plot_FRic_Litter <- ggplot(m1_FRic_Litter_pred, aes(Litter_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0))  +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Litter cover, %", 
       y = "Functional richness") 

plot_FRic_Litter

##  Bare_Ground_Cover  ---------

BareSoil_range = seq(min(Data_1m2$Bare_Ground_Cover),
                     max(Data_1m2$Bare_Ground_Cover), 
                     by = 0.001)

m1_FRic_BareSoil_pred <- Effect("Bare_Ground_Cover", 
                                m1_FRic,
                                xlevels = list(Bare_Ground_Cover = BareSoil_range)) %>% 
  as.data.frame() 

plot_FRic_BareSoil <- ggplot(m1_FRic_BareSoil_pred, 
                             aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.05, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Functional richness") 

plot_FRic_BareSoil


##  Slope effect   -------------------------------------------------------------

Slope_range = seq(min(Data_1m2$slope_degr),
                  max(Data_1m2$slope_degr), 
                  by = 0.001)

m1_FRic_slope_pred <- Effect("slope_degr", 
                             m1_FRic,
                             xlevels = list(slope_degr = Slope_range)) %>% 
  as.data.frame() 


plot_FRic_Slope <-ggplot(m1_FRic_slope_pred, aes(slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope", 
       y = "Functional richness") 


plot_FRic_Slope



## Ttree distance  ----------------------------------------------------------------
Tree_dist_range = seq(min(Data_1m2$dist_tree),
                      max(Data_1m2$dist_tree), 
                      by = 0.001)

m1_FRic_TreeDist_pred <- Effect("dist_tree", 
                                m1_FRic,
                                xlevels = list(dist_tree = Tree_dist_range)) %>% 
  as.data.frame() 


plot_FRic_TreeDist <-
  ggplot(m1_FRic_TreeDist_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Functional richness") 

plot_FRic_TreeDist


##  Sky view   -----------------------------------------------------------

Sky_range = seq(min(Data_1m2$sky_view_factor),
                max(Data_1m2$sky_view_factor), 
                by = 0.001)

m1_FRic_Sky_pred <- Effect("sky_view_factor", 
                           m1_FRic,
                           xlevels = list(sky_view_factor = Sky_range)) %>% 
  as.data.frame() 


plot_FRic_Sky <- ggplot(m1_FRic_Sky_pred, 
                        aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Functional richness") 


plot_FRic_Sky

##  Biotop_richness_specific  --------------------------------------------------


Biotop_range = seq(min(Data_1m2$Biotop_richness_specific),
                   max(Data_1m2$Biotop_richness_specific), 
                   by = 0.001)

m1_FRic_Biotop_pred <- Effect("Biotop_richness_specific", 
                              m1_FRic,
                              xlevels = list(Biotop_richness_specific = Biotop_range)) %>% 
  as.data.frame()

plot_FRic_Biotop <- ggplot(m1_FRic_Biotop_pred, 
                           aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.15, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Functional richness") 

plot_FRic_Biotop


##  Protected areas  ----------------------------------------------------------


ProtecArea_range = seq(min(Data_1m2$protected_cover_pct),
                       max(Data_1m2$protected_cover_pct), 
                       by = 0.001)

m1_FRic_Protect_pred <- Effect("protected_cover_pct", 
                               m1_FRic,
                               xlevels = list(protected_cover_pct = ProtecArea_range))%>% 
  as.data.frame() 



plot_FRic_ProtcArea <- ggplot(m1_FRic_Protect_pred, 
                              aes(x = protected_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(protected_cover_pct, FRic, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.1, height = 0.04)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Protected areas (%)",
       y = "Functional richness") 

plot_FRic_ProtcArea


## Management ---------------------------------------------------------------------

m1_FRic_no_int <- update(m1_FRic, . ~ . - MowFreq:Month) # we remove interaction for lotting main effcets
drop1(m1_FRic_no_int)

emmeans(m1_FRic_no_int, list(pairwise ~ MowFreq))

emmeans_m1_FRic_MowFreq1 <- cld(emmeans(m1_FRic_no_int, list(pairwise ~ MowFreq)), 
                                Letters = letters) %>% 
  arrange(MowFreq)

FRic_max1 <-  Data_1m2 %>% 
  summarise(max=max(FRic), .by = c(MowFreq))

plot_FRic_Management <- Data_1m2 %>% 
  ggplot(aes(x = MowFreq, y = FRic)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Management", y = "Functional richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m1_FRic_MowFreq1 %>% 
              left_join(FRic_max1, by=c("MowFreq")),
            aes(x=MowFreq, y=max+0.5,
                label= .group 
            ),
            size=3.5, col="black") 

plot_FRic_Management

## Month ----------------------------------------------------------------------

emmeans(m1_FRic_no_int, list(pairwise ~ Month))

emmeans_m1_FRic_Month <- cld(emmeans(m1_FRic_no_int, list(pairwise ~ Month)), 
                             Letters = letters) %>% 
  arrange(Month)

FRic_max2 <-  Data_1m2 %>% 
  summarise(max=max(FRic), .by = c(Month))

plot_FRic_Month <- Data_1m2 %>% 
  ggplot(aes(x = Month, y = FRic)) +
  theme_bw() +
  geom_point(aes(color=Month),
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.25, height = 0)) +
  geom_boxplot(outliers = F, alpha=0) +
  labs(x = "Month", y = "Functional richness") +
  scale_color_manual(values = Month_col) +
  # theme(legend.position = "none") +
  geom_text(data=emmeans_m1_FRic_Month %>% 
              left_join(FRic_max2, by=c("Month")),
            aes(x=Month, y=max+1,
                label=.group),
            size=3.5, col="black")


plot_FRic_Month

## MowFreq * Month  ------------------------------------------------------------

emmeans(m1_FRic, list(pairwise ~ MowFreq | Month))

emmeans_m1_FRic_MowFreq <- cld(emmeans(m1_FRic, list(pairwise ~ MowFreq | Month)), 
                               Letters = letters) %>% 
  arrange(MowFreq) %>% 
  as.tibble()

FRic_max <-  Data_1m2 %>% 
  summarise(max=max(FRic), .by = c(MowFreq, Month))

plot_FRic_Month_Manag_intr <- Data_1m2 %>% 
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
  geom_text(data=emmeans_m1_FRic_MowFreq %>% 
              left_join(FRic_max, by=c("MowFreq", "Month")),
            aes(x=MowFreq, y=max+5,
                label=.group),
            size=3.5, col="black") 

plot_FRic_Month_Manag_intr


## COMBINE PLOTS ---------------------------------------------------------------

library(patchwork)

combined_FRic <- 
  plot_FRic_Mow_ev + plot_FRic_Litter +
  plot_FRic_BareSoil + plot_FRic_Slope + 
  plot_FRic_TreeDist + plot_FRic_Sky + 
  plot_FRic_Biotop + plot_FRic_ProtcArea + 
  # plot_FRic_Management + plot_FRic_Month +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = 'bold', size = 15),
    plot.tag.position = c(0.1, 1.05),
    plot.margin = margin(t = 22, r = 10, b = 10, l = 10)
  )

print(combined_FRic)


ggsave("results/plots/FRic.png", combined_FRic, width = 7, height = 10, dpi = 450)



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
                            MowFreq * Month +  
                            scale(n_mow_events_befre_sampling) +
                            scale(Litter_Cover) + 
                            #      scale(road_density_km_per_ha) + 
                            #      scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(Biotop_richness_specific) + 
                               scale(protected_cover_pct) + 
                            # scale(green_cover_pct) + 
                            (1|PlotNo),
                          data = Data_1m2)


# check model assumptions
check_convergence(m1_FEve)
check_model(m1_FEve)
check_collinearity(m1_FEve)
# check interactions
drop1(m1_FEve)

# interaction is not significant
m2_FEve <- lmerTest::lmer(FEve ~ 
                              MowFreq + Month +  
                              scale(n_mow_events_befre_sampling) +
                              scale(Litter_Cover) + 
                              #      scale(road_density_km_per_ha) + 
                              #      scale(patch_size_m2) +
                              Bare_Ground_Cover + 
                              scale(slope_degr) + scale(dist_tree) + 
                              scale(sky_view_factor) +
                              scale(Biotop_richness_specific) + 
                                 scale(protected_cover_pct) + 
                            # scale(green_cover_pct) + 
                              (1|PlotNo),
                            data = Data_1m2)


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

##  Mowing events before sampling  ---------

mow_events_range = seq(min(Data_1m2$n_mow_events_befre_sampling),
                       max(Data_1m2$n_mow_events_befre_sampling), 
                       by = 0.001)

m2_FEve_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_FEve,
                           xlevels = list(n_mow_events_befre_sampling = mow_events_range)) %>% 
  as.data.frame()

plot_FEve_Mow_ev <- ggplot(m2_FEve_Mow_pred, 
                           aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.2, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Functional evenness") 

plot_FEve_Mow_ev

##  Litter_Cover  ---------

Litter_range = seq(min(Data_1m2$Litter_Cover),
                   max(Data_1m2$Litter_Cover), 
                   by = 0.01)

m2_FEve_Litter_pred <- Effect("Litter_Cover", 
                              m2_FEve,
                              xlevels = list(Litter_Cover = Litter_range)) %>% 
  as.data.frame()

plot_FEve_Litter <- ggplot(m2_FEve_Litter_pred, aes(Litter_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0))  +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Litter cover, %", 
       y = "Functional evenness") 

plot_FEve_Litter

##  Bare_Ground_Cover  ---------

BareSoil_range = seq(min(Data_1m2$Bare_Ground_Cover),
                     max(Data_1m2$Bare_Ground_Cover), 
                     by = 0.001)

m2_FEve_BareSoil_pred <- Effect("Bare_Ground_Cover", 
                                m2_FEve,
                                xlevels = list(Bare_Ground_Cover = BareSoil_range)) %>% 
  as.data.frame() 

plot_FEve_BareSoil <- ggplot(m2_FEve_BareSoil_pred, 
                             aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.05, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Functional evenness") 

plot_FEve_BareSoil


##  Slope effect   -------------------------------------------------------------

Slope_range = seq(min(Data_1m2$slope_degr),
                  max(Data_1m2$slope_degr), 
                  by = 0.001)

m2_FEve_slope_pred <- Effect("slope_degr", 
                             m2_FEve,
                             xlevels = list(slope_degr = Slope_range)) %>% 
  as.data.frame() 


plot_FEve_Slope <-ggplot(m2_FEve_slope_pred, aes(slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope", 
       y = "Functional evenness") 


plot_FEve_Slope



## Ttree distance  ----------------------------------------------------------------
Tree_dist_range = seq(min(Data_1m2$dist_tree),
                      max(Data_1m2$dist_tree), 
                      by = 0.001)

m2_FEve_TreeDist_pred <- Effect("dist_tree", 
                                m2_FEve,
                                xlevels = list(dist_tree = Tree_dist_range)) %>% 
  as.data.frame() 


plot_FEve_TreeDist <-
  ggplot(m2_FEve_TreeDist_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Functional evenness") 

plot_FEve_TreeDist


##  Sky view   -----------------------------------------------------------

Sky_range = seq(min(Data_1m2$sky_view_factor),
                max(Data_1m2$sky_view_factor), 
                by = 0.001)

m2_FEve_Sky_pred <- Effect("sky_view_factor", 
                           m2_FEve,
                           xlevels = list(sky_view_factor = Sky_range)) %>% 
  as.data.frame() 


plot_FEve_Sky <- ggplot(m2_FEve_Sky_pred, 
                        aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Functional evenness") 


plot_FEve_Sky

##  Biotop_richness_specific  --------------------------------------------------


Biotop_range = seq(min(Data_1m2$Biotop_richness_specific),
                   max(Data_1m2$Biotop_richness_specific), 
                   by = 0.001)

m2_FEve_Biotop_pred <- Effect("Biotop_richness_specific", 
                              m2_FEve,
                              xlevels = list(Biotop_richness_specific = Biotop_range)) %>% 
  as.data.frame()

plot_FEve_Biotop <- ggplot(m2_FEve_Biotop_pred, 
                           aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.15, height = 0)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Functional evenness") 

plot_FEve_Biotop


##  Protected areas  ----------------------------------------------------------


ProtecArea_range = seq(min(Data_1m2$protected_cover_pct),
                       max(Data_1m2$protected_cover_pct), 
                       by = 0.001)

m2_FEve_Protect_pred <- Effect("protected_cover_pct", 
                               m2_FEve,
                               xlevels = list(protected_cover_pct = ProtecArea_range))%>% 
  as.data.frame() 



plot_FEve_ProtcArea <- ggplot(m2_FEve_Protect_pred, 
                              aes(x = protected_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(protected_cover_pct, FEve, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.1, height = 0.04)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Protected areas (%)",
       y = "Functional evenness") 

plot_FEve_ProtcArea


## Management ---------------------------------------------------------------------

emmeans(m2_FEve, list(pairwise ~ MowFreq))

emmeans_m2_FEve_MowFreq1 <- cld(emmeans(m2_FEve, list(pairwise ~ MowFreq)), 
                                Letters = letters) %>% 
  arrange(MowFreq)

FEve_max1 <-  Data_1m2 %>% 
  summarise(max=max(FEve), .by = c(MowFreq))

plot_FEve_Management <- Data_1m2 %>% 
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
            aes(x=MowFreq, y=max+0.05,
                label=c("a", "b'", "a") 
                # label= .group 
            ),
            size=3.5, col="black") 

plot_FEve_Management

## Month ----------------------------------------------------------------------

emmeans(m2_FEve, list(pairwise ~ Month))

emmeans_m2_FEve_Month <- cld(emmeans(m2_FEve, list(pairwise ~ Month)), 
                             Letters = letters) %>% 
  arrange(Month)

FEve_max2 <-  Data_1m2 %>% 
  summarise(max=max(FEve), .by = c(Month))

plot_FEve_Month <- Data_1m2 %>% 
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


plot_FEve_Month

## MowFreq * Month  ------------------------------------------------------------

emmeans(m1_FEve, list(pairwise ~ MowFreq | Month))

emmeans_m2_FEve_MowFreq <- cld(emmeans(m1_FEve, list(pairwise ~ MowFreq | Month)), 
                               Letters = letters) %>% 
  arrange(MowFreq) %>% 
  as.tibble()

FEve_max <-  Data_1m2 %>% 
  summarise(max=max(FEve), .by = c(MowFreq, Month))

plot_FEve_Month_Manag_intr <- Data_1m2 %>% 
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
            aes(x=MowFreq, y=max+0.1,
                label=.group),
            size=3.5, col="black") 

plot_FEve_Month_Manag_intr


## COMBINE PLOTS ---------------------------------------------------------------

library(patchwork)

combined_FEve <- 
  plot_FEve_Mow_ev + plot_FEve_Litter +
  plot_FEve_BareSoil + plot_FEve_Slope + 
  plot_FEve_TreeDist + plot_FEve_Sky + 
  plot_FEve_Biotop + plot_FEve_ProtcArea + 
  # plot_FEve_Management + plot_FEve_Month +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = 'bold', size = 15),
    plot.tag.position = c(0.1, 1.05),
    plot.margin = margin(t = 22, r = 10, b = 10, l = 10)
  )

print(combined_FEve)


ggsave("results/plots/FEve.png", combined_FEve, width = 7, height = 10, dpi = 450)





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
                            #      scale(road_density_km_per_ha) + 
                            #      scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(Biotop_richness_specific) + 
                            scale(protected_cover_pct) + 
                            #  scale(green_cover_pct) + 
                           # scale(impervious_pct) +
                            (1|PlotNo),
                          data = Data_1m2)


# check model assumptions
check_convergence(m1_FDis)
check_model(m1_FDis)
check_collinearity(m1_FDis)
# check interactions
drop1(m1_FDis)

# interaction is  not significant
m2_FDis <- lmerTest::lmer(FDis ~ 
                            MowFreq + Month +  
                            scale(log1p(n_mow_events_befre_sampling)) +
                            scale(Litter_Cover) + 
                            #      scale(road_density_km_per_ha) + 
                            #      scale(patch_size_m2) +
                            Bare_Ground_Cover + 
                            scale(slope_degr) + scale(dist_tree) + 
                            scale(sky_view_factor) +
                            scale(log1p(Biotop_richness_specific)) + 
                            scale(protected_cover_pct) + 
                            #  scale(green_cover_pct) + 
                          #  scale(impervious_pct) +
                            (1|PlotNo),
                          data = Data_1m2)

summary(m2_FDis)
# check model assumptions

check_convergence(m2_FDis)
check_model(m2_FDis)
check_collinearity(m2_FDis)


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

##  Mowing events before sampling  ---------

mow_events_range = seq(min(Data_1m2$n_mow_events_befre_sampling),
                       max(Data_1m2$n_mow_events_befre_sampling), 
                       by = 0.001)

m2_FDis_Mow_pred <- Effect("n_mow_events_befre_sampling", m2_FDis,
                           xlevels = list(n_mow_events_befre_sampling = mow_events_range)) %>% 
  as.data.frame()

plot_FDis_Mow_ev <- ggplot(m2_FDis_Mow_pred, 
                           aes(x = n_mow_events_befre_sampling, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(n_mow_events_befre_sampling, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.2, height = 0.15)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Mowing events before sampling", 
       y = "Functional dispersion") 

plot_FDis_Mow_ev

##  Litter_Cover  ---------

Litter_range = seq(min(Data_1m2$Litter_Cover),
                   max(Data_1m2$Litter_Cover), 
                   by = 0.01)

m2_FDis_Litter_pred <- Effect("Litter_Cover", 
                              m2_FDis,
                              xlevels = list(Litter_Cover = Litter_range)) %>% 
  as.data.frame()

plot_FDis_Litter <- ggplot(m2_FDis_Litter_pred, aes(Litter_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Litter_Cover, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0, height = 0))  +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Litter cover, %", 
       y = "Functional dispersion") 

plot_FDis_Litter

##  Bare_Ground_Cover  ---------

BareSoil_range = seq(min(Data_1m2$Bare_Ground_Cover),
                     max(Data_1m2$Bare_Ground_Cover), 
                     by = 0.001)

m2_FDis_BareSoil_pred <- Effect("Bare_Ground_Cover", 
                                m2_FDis,
                                xlevels = list(Bare_Ground_Cover = BareSoil_range)) %>% 
  as.data.frame() 

plot_FDis_BareSoil <- ggplot(m2_FDis_BareSoil_pred, 
                             aes(x = Bare_Ground_Cover, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Bare_Ground_Cover, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.05, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Bare ground cover (%)",
       y = "Functional dispersion") 

plot_FDis_BareSoil


##  Slope effect   -------------------------------------------------------------

Slope_range = seq(min(Data_1m2$slope_degr),
                  max(Data_1m2$slope_degr), 
                  by = 0.001)

m2_FDis_slope_pred <- Effect("slope_degr", 
                             m2_FDis,
                             xlevels = list(slope_degr = Slope_range)) %>% 
  as.data.frame() 


plot_FDis_Slope <-ggplot(m2_FDis_slope_pred, aes(slope_degr, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(slope_degr, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Slope", 
       y = "Functional dispersion") 


plot_FDis_Slope



## Ttree distance  ----------------------------------------------------------------
Tree_dist_range = seq(min(Data_1m2$dist_tree),
                      max(Data_1m2$dist_tree), 
                      by = 0.001)

m2_FDis_TreeDist_pred <- Effect("dist_tree", 
                                m2_FDis,
                                xlevels = list(dist_tree = Tree_dist_range)) %>% 
  as.data.frame() 


plot_FDis_TreeDist <-
  ggplot(m2_FDis_TreeDist_pred, aes(x = dist_tree, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(dist_tree, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Tree distance (m)",
       y = "Functional dispersion") 

plot_FDis_TreeDist


##  Sky view   -----------------------------------------------------------

Sky_range = seq(min(Data_1m2$sky_view_factor),
                max(Data_1m2$sky_view_factor), 
                by = 0.001)

m2_FDis_Sky_pred <- Effect("sky_view_factor", 
                           m2_FDis,
                           xlevels = list(sky_view_factor = Sky_range)) %>% 
  as.data.frame() 


plot_FDis_Sky <- ggplot(m2_FDis_Sky_pred, 
                        aes(x = sky_view_factor, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(sky_view_factor, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.04, height = 0.15)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Sky-view factor",
       y = "Functional dispersion") 


plot_FDis_Sky

##  Biotop_richness_specific  --------------------------------------------------


Biotop_range = seq(min(Data_1m2$Biotop_richness_specific),
                   max(Data_1m2$Biotop_richness_specific), 
                   by = 0.001)

m2_FDis_Biotop_pred <- Effect("Biotop_richness_specific", 
                              m2_FDis,
                              xlevels = list(Biotop_richness_specific = Biotop_range)) %>% 
  as.data.frame()

plot_FDis_Biotop <- ggplot(m2_FDis_Biotop_pred, 
                           aes(Biotop_richness_specific, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(Biotop_richness_specific, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.15, height = 0)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Landscape heterogeneity", 
       y = "Functional dispersion") 

plot_FDis_Biotop


##  Protected areas  ----------------------------------------------------------


ProtecArea_range = seq(min(Data_1m2$protected_cover_pct),
                       max(Data_1m2$protected_cover_pct), 
                       by = 0.001)

m2_FDis_Protect_pred <- Effect("protected_cover_pct", 
                               m2_FDis,
                               xlevels = list(protected_cover_pct = ProtecArea_range))%>% 
  as.data.frame() 



plot_FDis_ProtcArea <- ggplot(m2_FDis_Protect_pred, 
                              aes(x = protected_cover_pct, y = fit)) +
  # add CI across all dataset
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_point(data=Data_1m2,
             aes(protected_cover_pct, FDis, color=Month), 
             pch=19, size=1.5, alpha=0.6, stroke = 0.8,
             position = position_jitter(width = 0.1, height = 0.04)) +
  geom_line(linewidth = 0.8, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  scale_fill_manual(values = Month_col) +
  labs(x = "Protected areas (%)",
       y = "Functional dispersion") 

plot_FDis_ProtcArea


## Management ---------------------------------------------------------------------

emmeans(m2_FDis, list(pairwise ~ MowFreq))

emmeans_m2_FDis_MowFreq1 <- cld(emmeans(m2_FDis, list(pairwise ~ MowFreq)), 
                                Letters = letters) %>% 
  arrange(MowFreq)

FDis_max1 <-  Data_1m2 %>% 
  summarise(max=max(FDis), .by = c(MowFreq))

plot_FDis_Management <- Data_1m2 %>% 
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
            aes(x=MowFreq, y=max+0.2,
                label=c("a", "b'", "a") 
                # label= .group 
            ),
            size=3.5, col="black") 

plot_FDis_Management

## Month ----------------------------------------------------------------------

emmeans(m2_FDis, list(pairwise ~ Month))

emmeans_m2_FDis_Month <- cld(emmeans(m2_FDis, list(pairwise ~ Month)), 
                             Letters = letters) %>% 
  arrange(Month)

FDis_max2 <-  Data_1m2 %>% 
  summarise(max=max(FDis), .by = c(Month))

plot_FDis_Month <- Data_1m2 %>% 
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
            aes(x=Month, y=max+0.2,
                label=.group),
            size=3.5, col="black")


plot_FDis_Month

## MowFreq * Month  ------------------------------------------------------------

emmeans(m1_FDis, list(pairwise ~ MowFreq | Month))

emmeans_m2_FDis_MowFreq <- cld(emmeans(m1_FDis, list(pairwise ~ MowFreq | Month)), 
                               Letters = letters) %>% 
  arrange(MowFreq) %>% 
  as.tibble()

FDis_max <-  Data_1m2 %>% 
  summarise(max=max(FDis), .by = c(MowFreq, Month))

plot_FDis_Month_Manag_intr <- Data_1m2 %>% 
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
            aes(x=MowFreq, y=max+0.25,
                label=.group),
            size=3.5, col="black") 

plot_FDis_Month_Manag_intr


## COMBINE PLOTS ---------------------------------------------------------------

library(patchwork)

combined_FDis <- 
  plot_FDis_Mow_ev + plot_FDis_Litter +
  plot_FDis_BareSoil + plot_FDis_Slope + 
  plot_FDis_TreeDist + plot_FDis_Sky + 
  plot_FDis_Biotop + plot_FDis_ProtcArea + 
  # plot_FDis_Management + plot_FDis_Month +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = 'bold', size = 15),
    plot.tag.position = c(0.1, 1.05),
    plot.margin = margin(t = 22, r = 10, b = 10, l = 10)
  )

print(combined_FDis)


ggsave("results/plots/FDis.png", combined_FDis, width = 7, height = 10, dpi = 450)





# 5) R2 figure -----------------------------------------------------------------

# Function to extract coefficient sign
get_coefs_sign <- function(model, response_name) {
  tryCatch({
    coefs_raw <- fixef(model)
    coefs <- if (is.list(coefs_raw)) unlist(coefs_raw) else coefs_raw
    coef_df <- data.frame(
      Driver = names(coefs),
      sign   = ifelse(coefs > 0, "positive", "negative"),
      stringsAsFactors = FALSE
    )
    coef_df %>% mutate(Response = response_name)
  }, error = function(e) {
    message(paste("Error with", response_name, ":", e$message))
    NULL
  })
}

# call for one model
read_csv("results/LMM_biomass.csv") %>% mutate(Response = "Biomass") %>% 
  select(Driver, p_value=`Pr(>F)`, Rsq, Response) %>% 
  left_join(
  get_coefs_sign(m1_biomass, "Biomass"),
  by = c("Driver", "Response"))
  

coef_sign <- bind_rows(
          get_coefs_sign(m1_biomass,"Biomass"),
          get_coefs_sign(m2_SR_dummy, "Species richness"),
          get_coefs_sign(m2_evenness, "Hill-Simpson diversity"),
          get_coefs_sign(m1_phen_Richness, "Phenological richness"),
          get_coefs_sign(m1_phen_evenness, "Phenological evenness"),
          get_coefs_sign(m1_FRic, "Functional richness"),
          get_coefs_sign(m2_FEve, "Functional evenness"),
          get_coefs_sign(m2_FDis, "Functional dispersion")) %>% 
  select( Driver, sign, Response) 

# Combine all models
all_r2 <- bind_rows(
  read_csv("results/LMM_biomass.csv") %>% mutate(Response = "Biomass") %>% 
    select(Driver, p_value=`Pr(>F)`, Rsq, Response), 
  read_csv("results/GLMM_SpecRich.csv")%>% mutate(Response = "Species richness")%>% 
    select(Driver, p_value=`Pr(Chi)`, Rsq, Response), 
  read_csv("results/LMM_evenness.csv")%>% mutate(Response = "Hill-Simpson diversity")%>% 
    select(Driver, p_value=`Pr(>F)`, Rsq, Response),
  read_csv("results/LMM_Phenol_Richness.csv")%>% mutate(Response = "Phenological richness")%>% 
    select(Driver, p_value=`Pr(>F)`, Rsq, Response),
  read_csv("results/LMM_Phen_evenness.csv")%>% mutate(Response = "Phenological evenness")%>% 
    select(Driver, p_value=`Pr(>F)`, Rsq, Response),
  read_csv("results/LMM_Functional_richness.csv")%>% mutate(Response = "Functional richness")%>% 
    select(Driver, p_value=`Pr(>F)`, Rsq, Response),
  read_csv("results/LMM_Funct_evenness.csv")%>% mutate(Response = "Functional evenness")%>% 
    select(Driver, p_value=`Pr(>F)`, Rsq, Response),
  read_csv("results/LMM_Functional_dispersion.csv")%>% mutate(Response = "Functional dispersion")%>% 
    select(Driver, p_value=`Pr(>F)`, Rsq, Response)) %>% 
  left_join(coef_sign, by = c("Driver", "Response")) %>%
  mutate(Driver_clean = case_when(
    grepl("MowFreq:Month|Month:MowFreq", Driver) ~ "Management × Season",
    grepl("MowFreq", Driver) ~ "Management",
    grepl("Month", Driver) ~ "Season",
    grepl("Litter_Cover", Driver) ~ "Litter cover",
    grepl("n_mow_events", Driver) ~ "Mowing events",
    grepl("protected_cover_pct", Driver) ~ "Protected areas",
    grepl("protected_cover_scaled", Driver) ~ "Protected areas",
    grepl("Biotop_richness", Driver) ~ "Landscape heterogeneity",
    grepl("slope_degr", Driver) ~ "Slope",
    grepl("dist_tree", Driver) ~ "Tree distance",
    grepl("sky_view", Driver) ~ "Sky-view factor",
    grepl("Bare_Ground", Driver) ~ "Bare ground",
    TRUE ~ Driver)) %>% 
  # Set categorical variables
  mutate(sign = case_when(
           Driver_clean %in% c("Management × Season", "Management", "Season") ~ "categorical",
           TRUE ~ sign
         )) %>% 
  rename(partial_R2 = Rsq) %>% 
  mutate(significance=ifelse(p_value < 0.05, "significant", "non-significant")) 
  
all_r2 %>% 
  filter(is.na(sign))

all_r2 %>% pull(Driver_clean) %>% unique()
  
# Set order of responses
all_r2 <- all_r2 |>
  mutate(Response = factor(Response, levels = c(
    "Biomass", "Species richness", "Hill-Simpson diversity",
    "Phenological richness", "Phenological evenness",
    "Functional richness", "Functional evenness", "Functional dispersion"
  )))

# Set order of drivers
all_r2 <- all_r2 |>
  mutate(Driver_clean = factor(Driver_clean, levels = c(
    "Management × Season",
    "Management",
    "Season",
    "Mowing events",
    "Litter cover",
    "Bare ground",
    "Slope",
    "Tree distance",
    "Sky-view factor",
    "Patch size",
    "Landscape heterogeneity",
    "Protected areas"
  )))


all_r2 %>% 
  select(Driver_clean, partial_R2, sign, Response) %>% 
  pivot_longer(cols = c(Response), names_to = "metric", values_to = "value") 

  
  # Plot with color by sign
R2_plot <- ggplot(all_r2 %>% 
         mutate(sign=case_when(
           sign == "positive" ~ "positive effect",
           sign == "negative" ~ "negative effect",
           TRUE ~ "categorical driver")), 
       aes(x = Response, y = Driver_clean, 
           fill = sign, size = partial_R2, alpha=significance)) +
  geom_point(shape = 21, stroke = 0.5) +
  scale_fill_manual(
    values = c("negative effect" = "#B2182B", "positive effect" = "#2166AC", 
               "categorical driver" = "olivedrab"),
    name = "Effect type"
  ) +
  scale_alpha_manual(values = c("significant" = 1, "non-significant" = 0.5), 
                     name = "Effect significance") +
  scale_size_continuous(
    range = c(1, 8),
    name = expression(Partial~R^2)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color = "black", size=10),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "right") +
  guides(
    fill = guide_legend(override.aes = list(size = 5)),
    alpha = guide_legend(override.aes = list(
      size   = 5,
      fill   = c("grey69", "grey75"),   # dark = significant, light = non-significant
      colour = "black",
      alpha  = c(1, 0.5)
    ))
  ) +
  labs(x = "Plant community variables", y = "Drivers")

R2_plot
ggsave("results/plots/R2_plot.png", R2_plot, width = 8, height = 6, dpi = 450)
