# Explorations of possible drivers of diversity patterns in the 1m2 data, including cover data and GIS data.

library(tidyverse)

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



# assigning colors for the plots
MowFreq_col <- c("#F8766D", "#00B0F6","#00BA38")
Month_col <-  c( "orange", "#287271", "#6D326D","brown") 


# Correlations among response variables------------------------------------------------------------------


# Functional diversity:
cor <- round(cor(Funct_div %>% 
                   dplyr::select(-c(PlotNo, Subplot, Month)),
                 method = c("pearson"), use = "pairwise.complete.obs"), 2)

ggcorrplot::ggcorrplot(cor,
                       hc.order = F, type = "lower",
                       lab = TRUE, lab_size = 5,
                       colors = c("red", "white", "blue"))



corl1 <- round(cor(Dat1_1m2 %>% 
                     dplyr::select(cover, biomass, SR, 
                                   evenness, shannon,
                                   phen_Richness, phen_Shannon, phen_evenness,
                                   FRic,   FEve,  FDiv,  FDis) %>% 
                     rename("Total cover" =cover,
                            "Total biomass" =biomass,
                            "Species richness" =SR,
                            "Hill-Shannon diversity (specieis)" =shannon,
                            "Phenological richness" =phen_Richness,
                            "Hill-Shannon diversity (phenological)" =phen_Shannon,
                            "Hill-Simpson (phenological)" =phen_evenness,
                            "Hill-Simpson (specieis)" =evenness,
                            "Functional richness" =FRic,
                            "Functional evenness" =FEve,
                            "Functional divergence" =FDiv,
                            "Functional dispersion" =FDis),
                   method = c("pearson"), use = "pairwise.complete.obs"), 2)

corl1

ggcorrplot::ggcorrplot(corl1,
                       hc.order = F, type = "lower",
                       lab = TRUE, lab_size = 4,
                       colors = c("red", "white", "blue"))


# Exploratory statistics -------------------------------------------------------

## 1) With Cover data: ----------------------------------------------------


#### SR -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, SR) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = SR)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Species Richness")



#### Biomass -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = biomass)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Biomass")

#### Phenological Rich -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Richness) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Richness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Richness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Richness")

#### Phenological evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_evenness) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_evenness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_evenness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological evenness")


#### Phenological Shannon  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Shannon) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Shannon), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Shannon)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Shannon")



#### Functional richness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FRic) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FRic), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FRic)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional richness")




#### Functional evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FEve) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FEve), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FEve)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional evenness")


#### Functional dispersion  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FDis) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FDis), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FDis)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional dispersion")



#### Neophytes_SR_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, SR, Neophytes_SR_propr) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR, Neophytes_SR_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_SR_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = binomial),
              aes(weight = SR),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on species richness)")



#### Neophytes_mass_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass, Neophytes_mass_propr) %>% 
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass, Neophytes_mass_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_mass_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", aes(weight = biomass),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on biomass)")



## 2) With Carolin's GIS data: ----------------------------------------------------

names(Dat1_1m2)

### 2.1) Field measured:  ---------------------------------------------------------------

#### SR -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, SR) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = SR)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Species Richness")


#### Biomass -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = biomass)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Biomass")

#### Phenological Rich -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Richness) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Richness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Richness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Richness")

#### Phenological evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_evenness) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_evenness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_evenness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological evenness")


#### Phenological Shannon  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Shannon) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Shannon), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Shannon)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Shannon")



#### Functional richness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FRic) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FRic), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FRic)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional richness")




#### Functional evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FEve) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FEve), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FEve)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional evenness")


#### Functional dispersion  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FDis) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FDis), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FDis)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional dispersion")



#### Neophytes_SR_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, SR, Neophytes_SR_propr) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR, Neophytes_SR_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_SR_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = binomial),
              aes(weight = SR),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on species richness)")



#### Neophytes_mass_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass, Neophytes_mass_propr) %>% 
  left_join(Fild_data %>% 
              dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass, Neophytes_mass_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_mass_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", aes(weight = biomass),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on biomass)")


### 2.2) GIS main (patch type, size, distances) :  ---------------------------------------------------------------
names(GIS_main)

#### SR -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month,SR) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = SR)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Species Richness")


#### Biomass -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = log(biomass))) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Biomass")

#### Phenological Rich -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Richness) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Richness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Richness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Richness")

#### Phenological evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_evenness) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_evenness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_evenness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological evenness")


#### Phenological Shannon  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Shannon) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Shannon), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Shannon)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Shannon")



#### Functional richness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FRic) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FRic), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FRic)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional richness")




#### Functional evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FEve) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FEve), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FEve)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional evenness")


#### Functional dispersion  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FDis) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FDis), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FDis)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional dispersion")



#### Neophytes_SR_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, SR, Neophytes_SR_propr) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR, Neophytes_SR_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_SR_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = binomial),
              aes(weight = SR),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on species richness)")



#### Neophytes_mass_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass, Neophytes_mass_propr) %>% 
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass, Neophytes_mass_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_mass_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", aes(weight = biomass),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on biomass)")





### 2.3) GIS in 500m radius:  ---------------------------------------------------------------

#### SR -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month,SR) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = SR)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Species Richness", title = "500m buffer") +
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme(legend.position = "bottom")


#### Biomass -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = log(biomass))) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Biomass", title = "500m buffer") +   
  guides(color = guide_legend(override.aes = list(size = 4)))+     
  theme(legend.position = "bottom")

#### Phenological Rich -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Richness) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Richness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Richness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Richness", title = "500m buffer") +   
  guides(color = guide_legend(override.aes = list(size = 4)))+     
  theme(legend.position = "bottom")

#### Phenological evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_evenness) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_evenness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_evenness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological evenness", title = "500m buffer") +  
  guides(color = guide_legend(override.aes = list(size = 4)))+    
  theme(legend.position = "bottom")


#### Phenological Shannon  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Shannon) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Shannon), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Shannon)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Shannon", title = "500m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")



#### Functional richness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FRic) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FRic), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FRic)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional richness", title = "500m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")




#### Functional evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FEve) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FEve), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FEve)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional evenness", title = "500m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")


#### Functional dispersion  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FDis) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FDis), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FDis)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional dispersion", title = "500m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")



#### Neophytes_SR_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, SR, Neophytes_SR_propr) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR, Neophytes_SR_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_SR_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = binomial),
              aes(weight = SR),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on species richness)", 
       title = "500m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")



#### Neophytes_mass_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass, Neophytes_mass_propr) %>% 
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass, Neophytes_mass_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_mass_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", aes(weight = biomass),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on biomass)", title = "500m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")

### 2.4) GIS in 1000m radius:  ---------------------------------------------------------------

#### SR -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month,SR) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = SR)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Species Richness", title = "1000m buffer") +
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme(legend.position = "bottom")


#### Biomass -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = log(biomass))) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Biomass", title = "1000m buffer") +   
  guides(color = guide_legend(override.aes = list(size = 4)))+     
  theme(legend.position = "bottom")

#### Phenological Rich -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Richness) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Richness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Richness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Richness", title = "1000m buffer") +   
  guides(color = guide_legend(override.aes = list(size = 4)))+     
  theme(legend.position = "bottom")

#### Phenological evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_evenness) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_evenness), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_evenness)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological evenness", title = "1000m buffer") +  
  guides(color = guide_legend(override.aes = list(size = 4)))+    
  theme(legend.position = "bottom")


#### Phenological Shannon  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, phen_Shannon) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, phen_Shannon), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = phen_Shannon)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Phenological Shannon", title = "1000m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")



#### Functional richness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FRic) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FRic), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FRic)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional richness", title = "1000m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")




#### Functional evenness  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FEve) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FEve), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FEve)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional evenness", title = "1000m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")


#### Functional dispersion  -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, FDis) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, FDis), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = FDis)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", 
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Functional dispersion", title = "1000m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")



#### Neophytes_SR_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, SR, Neophytes_SR_propr) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, SR, Neophytes_SR_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_SR_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = binomial),
              aes(weight = SR),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on species richness)", 
       title = "1000m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")



#### Neophytes_mass_propr -----
Dat1_1m2 %>%
  dplyr::select(PlotNo, Subplot,  MowFreq, Month, biomass, Neophytes_mass_propr) %>% 
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, biomass, Neophytes_mass_propr), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = Neophytes_mass_propr)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "lm", aes(weight = biomass),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Neophytes proportion (based on biomass)", title = "1000m buffer") +   guides(color = guide_legend(override.aes = list(size = 4)))+     theme(legend.position = "bottom")



# 3) Predictor correlations: ----------------------------------------------------

plot_data <- Mowing_data %>% 
  mutate(MowFreq=ifelse(MowFreq == "reduced_sown", "reduced & sowing", MowFreq)) %>% 
  mutate(MowFreq=fct_relevel(MowFreq,"regular", 
                             "reduced", 
                             "reduced & sowing")) %>% 
  mutate(Month=fct_relevel(Month,"March", "May", "July", "September")) %>% 
  mutate(PlotNo=factor(PlotNo))



plot_data %>%
  ggplot(aes(x = MowFreq, y = n_mow_events_befre_sampling, color=Month)) +
  geom_point(aes(color=Month), alpha=0.6,
             position=position_jitterdodge(jitter.width=0.4, jitter.height=0.07,
                                           dodge.width=0.7)) +
  geom_boxplot(outliers = F, alpha=0) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = "Management", y = "Number of mowing events before sampling")



# 3.1) with Cover data -------------------------------------------------------------
#### Mowing events number  -----
plot_daa %>%
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, n_mow_events_befre_sampling), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = n_mow_events_befre_sampling)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Number of mowing events before sampling")


#### Management  -----
plot_data %>%
  left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, n_mow_events_befre_sampling), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = MowFreq, y = value, color=Month)) +
  geom_point(aes(color=Month), alpha=0.6,
             position=position_jitterdodge(jitter.width=0.4, jitter.height=0.07,
                                           dodge.width=0.7)) +
  geom_boxplot(outliers = F, alpha=0) +
  facet_wrap(~ variable, scales = "free",
             ncol=3) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = "Management", y = "Driver") +
  theme(legend.position = "bottom")





## 3) Predictor correlations: ----------------------------------------------------

# 3.2) with Carolin's Field measured -------------------------------------------------------------

#### Mowing events number  -----
plot_data %>%
  left_join(Fild_data %>% dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, n_mow_events_befre_sampling), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = n_mow_events_befre_sampling)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Number of mowing events before sampling")


#### Management  -----
plot_data %>%
  group_by(PlotNo, MowFreq, Subplot) %>%
  count() %>% dplyr::select(-n) %>% ungroup() %>%
  left_join(Fild_data %>% dplyr::select(-xcoord, -ycoord), by = c("PlotNo", "Subplot")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = MowFreq, y = value, color=MowFreq)) +
  geom_point(aes(color=MowFreq), alpha=0.6,
             position=position_jitterdodge(jitter.width=0.4, jitter.height=0.07,
                                           dodge.width=0.7)) +
  geom_boxplot(outliers = F, alpha=0) +
  facet_wrap(~ variable, scales = "free",
             ncol=3) +
  theme_bw() +
  scale_color_manual(values = MowFreq_col) +
  labs(x = "Management", y = "Driver") +
  theme(legend.position = "bottom")




# 3.3) with Carolin's GIS main (patch type, size, distances) :  ----------------

#### Mowing events number  -----
plot_data %>%
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>% 
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, n_mow_events_befre_sampling), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = n_mow_events_befre_sampling)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Number of mowing events before sampling")


#### Management  -----
plot_data %>%
  group_by(PlotNo, MowFreq) %>%
  count() %>% dplyr::select(-n) %>% ungroup() %>%
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  left_join(Sky_view, by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>% 
  left_join(Biotop_diversity, by = c("PlotNo")) %>% 
  pivot_longer(-c(PlotNo, MowFreq), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = MowFreq, y = value, color=MowFreq)) +
  geom_point(aes(color=MowFreq), alpha=0.6,
             position=position_jitterdodge(jitter.width=0.4, jitter.height=0.07,
                                           dodge.width=0.7)) +
  geom_boxplot(outliers = F, alpha=0) +
  facet_wrap(~ variable, scales = "free",
             ncol=3) +
  theme_bw() +
  scale_color_manual(values = MowFreq_col) +
  labs(x = "Management", y = "Driver") +
  theme(legend.position = "bottom")


plot_data %>%
  group_by(PlotNo, MowFreq) %>%
  count() %>% dplyr::select(-n) %>% ungroup() %>%
  left_join(GIS_main %>% 
              dplyr::select(PlotNo, patch_biotope_area_sqm, min_dist_to_street_m,
                            min_dist_to_water_m, min_dist_to_business_district_km_no_buffer), 
            by = c("PlotNo")) %>%
  rename(patch_size_m2=patch_biotope_area_sqm) %>%
  ggplot(aes(x = MowFreq, y = patch_size_m2, color=MowFreq)) +
  geom_point(aes(color=MowFreq), alpha=0.6,
             position=position_jitterdodge(jitter.width=0.4, jitter.height=0.07,
                                           dodge.width=0.7)) +
  geom_boxplot(outliers = F, alpha=0) +
  theme_bw() +
  scale_color_manual(values = MowFreq_col) +
  labs(x = "Management", y = "Patch size (m2)") +
  theme(legend.position = "bottom")

# 3.4) 500m: with Carolin's GIS data in 500m buffer :  ----------------

#### Mowing events number  -----
plot_data %>%
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, n_mow_events_befre_sampling), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = n_mow_events_befre_sampling)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Number of mowing events before sampling", title = "500m buffer") +   
  guides(color = guide_legend(override.aes = list(size = 4)))+     
  theme(legend.position = "bottom")

#### Management  -----
plot_data %>%
  group_by(PlotNo, MowFreq) %>%
  count() %>% dplyr::select(-n) %>% ungroup() %>%
  left_join(GIS_500, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, MowFreq), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = MowFreq, y = value, color=MowFreq)) +
  geom_point(aes(color=MowFreq), alpha=0.6,
             position=position_jitterdodge(jitter.width=0.4, jitter.height=0.07,
                                           dodge.width=0.7)) +
  geom_boxplot(outliers = F, alpha=0) +
  facet_wrap(~ variable, scales = "free",
             ncol=3) +
  theme_bw() +
  scale_color_manual(values = MowFreq_col) +
  labs(x = "Management", y = "Driver", title = "500m buffer") +
  theme(legend.position = "bottom")


#### among the GIS data -----

GIS_500 %>% 
   ggplot(aes(impervious_pct, green_cover_pct)) +
  geom_point(size=2, color="gray44") +
  geom_smooth(method = "lm", 
              se = TRUE, color = "black") +
 labs(x = "Impervious surfaces cover, %", 
      y = "Green surfaces cover, %") +
  theme_bw()


GIS_500 %>% 
  ggplot(aes(protected_biotopes_polygons_cover_pct, green_cover_pct)) +
  geom_point(size=2, color="gray44") +
  geom_smooth(method = "lm", 
              se = TRUE, color = "black") +
  labs(y = "Protected biotopes cover, %", 
       x = "Green surfaces cover, %") +
  theme_bw()


# 3.5) 1000m: with Carolin's GIS data in 1000m buffer :  ----------------

#### Mowing events number  -----
plot_data %>%
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, Subplot, MowFreq, Month, n_mow_events_befre_sampling), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = n_mow_events_befre_sampling)) +
  geom_point(aes(color=Month), alpha=0.6) +
  geom_smooth(method = "glm", method.args = list(family = poisson),
              se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = NULL, y = "Number of mowing events before sampling", title = "1000m buffer") +   
  guides(color = guide_legend(override.aes = list(size = 4)))+     
  theme(legend.position = "bottom")

#### Management  -----
plot_data %>%
  group_by(PlotNo, MowFreq) %>%
  count() %>% dplyr::select(-n) %>% ungroup() %>%
  left_join(GIS_1000, by = c("PlotNo")) %>%
  pivot_longer(-c(PlotNo, MowFreq), 
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = MowFreq, y = value, color=MowFreq)) +
  geom_point(aes(color=MowFreq), alpha=0.6,
             position=position_jitterdodge(jitter.width=0.4, jitter.height=0.07,
                                           dodge.width=0.7)) +
  geom_boxplot(outliers = F, alpha=0) +
  facet_wrap(~ variable, scales = "free",
             ncol=3) +
  theme_bw() +
  scale_color_manual(values = MowFreq_col) +
  labs(x = "Management", y = "Driver", title = "1000m buffer") +
  theme(legend.position = "bottom")
