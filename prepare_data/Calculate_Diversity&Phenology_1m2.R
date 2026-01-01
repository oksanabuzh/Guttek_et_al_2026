# Purpuse: Calculate community diversity and phenology metrics at 1m2 plot level


library(tidyverse)

dat <- read_csv("data/processed_data/vegetation2025_1m2.csv")


# Diversity metrics ----------------------------------------------------------
Diversity_phenology_1m2 <- dat %>%
  # calculate cover weighted phenology for each species for each plot and month
  mutate(
    Seedling = cover * Seedling,
    Juvenile = cover * Juvenile,
    FlowerBud     = cover * FlowerBud,
    Flowering     = cover * Flowering,
    Fruiting = cover * Fruiting,
    PostFruiting  = cover * PostFruiting
  ) %>%
  # calculate diversity and community phenology weighted by cover 
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
            .by=c("PlotNo", "Month", "Subplot"))

Diversity_phenology_1m2

write_csv(Diversity_phenology_1m2, "data/processed_data/Diversity_phenology_1m2.csv")
