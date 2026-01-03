# Purpose: Calculate community diversity and phenology metrics at 1m2 plot level


library(tidyverse)

dat <- read_csv("data/processed_data/vegetation2025_1m2.csv")

traits <- read_csv("data/processed_data/Traits_Dist.Ind.Values.csv")
names(traits)

# Community composition: species & phenological ----------------------------------------------------------
Composition_1m2 <- dat %>%
   left_join(traits %>% select(-EuroMed, -status_detailed),
             by="Taxon") %>%
  # calculate cover-weighted phenology for each species for each plot and month
  mutate(Neophyte=ifelse(status=="Neophyte", 1, 0),
         Endangered=ifelse(status=="Endangered", 1, 0))%>%
  mutate(
    Seedling = cover * Seedling,
    Juvenile = cover * Juvenile,
    FlowerBud     = cover * FlowerBud,
    Flowering     = cover * Flowering,
    Fruiting = cover * Fruiting,
    PostFruiting  = cover * PostFruiting,
    Annuals = cover * lifespan_annual,
    ShortLived = cover * lifespan_biennialOrShortLived,
    Perennials = cover * lifespan_perennial,
    Neophytes = cover * Neophyte,
    Endangereds = cover * Endangered) 

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
            Seedling_cover = sum(Seedling, na.rm = TRUE),
            Juvenile_cover = sum(Juvenile, na.rm = TRUE),
            FlowerBud_cover = sum(FlowerBud, na.rm = TRUE),
            Flowering_cover = sum(Flowering, na.rm = TRUE),
            Fruiting_cover = sum(Fruiting, na.rm = TRUE),
            PostFruiting_cover = sum(PostFruiting, na.rm = TRUE),
            Annuals_cover = sum(Annuals, na.rm = TRUE),
            ShortLived_cover = sum(ShortLived, na.rm = TRUE),
            Perennials_cover = sum(Perennials, na.rm = TRUE),
            Neophytes_cover = sum(Neophytes, na.rm = TRUE),
            Endangereds_cover = sum(Endangereds, na.rm = TRUE),
            Seedling_SR  = n_distinct(ifelse(Seedling>0, Taxon, NA), na.rm = TRUE),
            Juvenile_SR  = n_distinct(ifelse(Juvenile>0, Taxon, NA), na.rm = TRUE),
            FlowerBud_SR  = n_distinct(ifelse(FlowerBud>0, Taxon, NA), na.rm = TRUE),
            Flowering_SR  = n_distinct(ifelse(Flowering>0, Taxon, NA), na.rm = TRUE),
            Fruiting_SR  = n_distinct(ifelse(Fruiting>0, Taxon, NA), na.rm = TRUE),
            PostFruiting_SR  = n_distinct(ifelse(PostFruiting>0, Taxon, NA), na.rm = TRUE),
            Annuals_SR  = n_distinct(ifelse(Annuals>0, Taxon, NA), na.rm = TRUE),
            ShortLived_SR  = n_distinct(ifelse(ShortLived>0, Taxon, NA), na.rm = TRUE),
            Perennials_SR  = n_distinct(ifelse(Perennials>0, Taxon, NA), na.rm = TRUE),
            Neophytes_SR  = n_distinct(ifelse(Neophytes>0, Taxon, NA), na.rm = TRUE),
            Endangereds_SR  = n_distinct(ifelse(Endangereds>0, Taxon, NA), na.rm = TRUE),
            .by=c("PlotNo", "Month", "Subplot")) %>% 
  mutate(biomass=ifelse(is.na(height_mean), NA, biomass), # if height is NA, biomass should be NA and not 0
         height_max=ifelse(is.na(height_mean), NA, height_max)) %>%  # if height is NA, height_max should be NA and not Inf
  mutate(Seedling_cover_propor = Seedling_cover / biomass,
         Juvenile_cover_propor = Juvenile_cover / biomass,
         FlowerBud_cover_propor = FlowerBud_cover / biomass,
         Flowering_cover_propor = Flowering_cover / biomass,
         Fruiting_cover_propor = Fruiting_cover / biomass,
         PostFruiting_cover_propor = PostFruiting_cover / biomass,
         Annuals_cover_propor = Annuals_cover / biomass,
         ShortLived_cover_propor = ShortLived_cover / biomass,
         Perennials_cover_propor = Perennials_cover / biomass,
         Neophytes_cover_propor = Neophytes_cover / biomass,
         Endangereds_cover_propor = Endangereds_cover / biomass,
         Seedling_SR_propor = Seedling_SR / SR,
         Juvenile_SR_propor = Juvenile_SR / SR,
         FlowerBud_SR_propor = FlowerBud_SR / SR,
         Flowering_SR_propor = Flowering_SR / SR,
         Fruiting_SR_propor = Fruiting_SR / SR,
         PostFruiting_SR_propor = PostFruiting_SR / SR,
         Annuals_SR_propor = Annuals_SR / SR,
         ShortLived_SR_propor = ShortLived_SR / SR,
         Perennials_SR_propor = Perennials_SR / SR,
         Neophytes_SR_propor = Neophytes_SR / SR,
         Endangereds_SR_propor = Endangereds_SR / SR)

Diversity_phenology_1m2

write_csv(Diversity_phenology_1m2, "data/processed_data/Diversity_phenology_1m2.csv")
