# Phenology

library(tidyverse)
library(vegan)

# data
phen <- read_csv("data/Phenology1m2.csv") 

Phenophase_compos <- phen %>%
    mutate(
    Seedling = cover * Seedling,
    Juvenile = cover * Juvenile,
    FlowerBud     = cover * FlowerBud,
    Flowering     = cover * Flowering,
    Fruiting = cover * Fruiting,
    PostFruiting  = cover * PostFruiting
  ) %>%
  group_by(PlotNo, MowFreq, MonthLabel, Subplot) %>%
  summarise(
    Seedling = sum(Seedling, na.rm = TRUE),
    Juvenile = sum(Juvenile, na.rm = TRUE),
    FlowerBud = sum(FlowerBud, na.rm = TRUE),
    Flowering = sum(Flowering, na.rm = TRUE),
    Fruiting = sum(Fruiting, na.rm = TRUE),
    PostFruiting = sum(PostFruiting, na.rm = TRUE)
  ) %>%
  rowwise() %>%
  mutate(
    phen_Richness = sum(c_across(c(Seedling, Juvenile, FlowerBud, 
                                   Flowering, Fruiting, PostFruiting)) > 0),
    
    phen_Shannon = vegan::diversity(
      c_across(c(Seedling, Juvenile, FlowerBud, 
                 Flowering, Fruiting, PostFruiting)), index = "shannon"),
    phen_evenness = vegan::diversity(
      c_across(c(Seedling, Juvenile, FlowerBud, 
                 Flowering, Fruiting, PostFruiting)), index = "invsimpson")
    ) %>%
  ungroup()

Phenophase_compos 

write_csv(Phenophase_compos, "data/Phenological_composition.csv")
