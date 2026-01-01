# Purpose: Calculate community diversity and phenology metrics at 1m2 plot level


library(tidyverse)

dat <- read_csv("data/processed_data/vegetation2025_1m2.csv")


# Community composition: species & phenological ----------------------------------------------------------
Composition_1m2 <- dat %>%
  # calculate cover-weighted phenology for each species for each plot and month
  mutate(
    Seedling = cover * Seedling,
    Juvenile = cover * Juvenile,
    FlowerBud     = cover * FlowerBud,
    Flowering     = cover * Flowering,
    Fruiting = cover * Fruiting,
    PostFruiting  = cover * PostFruiting
  ) 
# readme for Composition_1m2: 

Composition_1m2

write_csv(Composition_1m2, "data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv")

# Diversity metrics ----------------------------------------------------------
# calculate diversity and community phenology weighted by cover 

Diversity_phenology_1m2 <- Composition_1m2 %>% 
  summarise(SR  = n_distinct(Taxon),
            evenness = vegan::diversity(cover, index = "invsimpson"),
            shannon = exp(vegan::diversity(cover, index = "shannon")),
            cover = sum(cover, na.rm = TRUE),
            biomass = sum(Biomass, na.rm = TRUE),
            height_mean=mean(height, na.rm = TRUE),
            height_max=max(height, na.rm = TRUE),
            Seedling = sum(Seedling, na.rm = TRUE),
            Juvenile = sum(Juvenile, na.rm = TRUE),
            FlowerBud = sum(FlowerBud, na.rm = TRUE),
            Flowering = sum(Flowering, na.rm = TRUE),
            Fruiting = sum(Fruiting, na.rm = TRUE),
            PostFruiting = sum(PostFruiting, na.rm = TRUE),
            .by=c("PlotNo", "Month", "Subplot")) %>% 
  mutate(biomass=ifelse(is.na(height_mean), NA, biomass), # if height is NA, biomass should be NA and not 0
         height_max=ifelse(is.na(height_mean), NA, height_max)) # if height is NA, height_max should be NA and not Inf

Diversity_phenology_1m2

write_csv(Diversity_phenology_1m2, "data/processed_data/Diversity_phenology_1m2.csv")
