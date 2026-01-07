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
    Vegetative = cover * Vegetative,
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
            Vegetative_cover = sum(Vegetative, na.rm = TRUE),
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
            Vegetative_SR  = n_distinct(ifelse(Vegetative>0, Taxon, NA), na.rm = TRUE),
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
  mutate(Seedling_cover_propor = Seedling_cover / cover,
         Juvenile_cover_propor = Juvenile_cover / cover,
         FlowerBud_cover_propor = FlowerBud_cover / cover,
         Flowering_cover_propor = Flowering_cover / cover,
         Fruiting_cover_propor = Fruiting_cover / cover,
         PostFruiting_cover_propor = PostFruiting_cover / cover,
         Annuals_cover_propor = Annuals_cover / cover,
         ShortLived_cover_propor = ShortLived_cover / cover,
         Perennials_cover_propor = Perennials_cover / cover,
         Neophytes_cover_propor = Neophytes_cover / cover,
         Endangereds_cover_propor = Endangereds_cover / cover,
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
         Endangereds_SR_propor = Endangereds_SR / SR) %>% 
  # Phenological diversity 
  rowwise() %>%
  mutate(
    phen_Richness = sum(c_across(c(Vegetative_cover,
                                   Seedling_cover,	Juvenile_cover,FlowerBud_cover,
                                   Flowering_cover,	Fruiting_cover,	PostFruiting_cover)) > 0),
    
    phen_Shannon = exp(vegan::diversity(
      c_across(c(Vegetative_cover, Seedling_cover,	Juvenile_cover,
                 FlowerBud_cover, Flowering_cover,	Fruiting_cover,	
                 PostFruiting_cover)), index = "shannon")),
    phen_evenness = vegan::diversity(
      c_across(c(Vegetative_cover, Seedling_cover,	Juvenile_cover,
                 FlowerBud_cover, Flowering_cover,	Fruiting_cover,	
                 PostFruiting_cover)), index = "invsimpson"),
    .after = "height_max") %>% 
  ungroup()

Diversity_phenology_1m2

#  Warning message concerns height_max, it get -Inf as no measurements were present for the plot


write_csv(Diversity_phenology_1m2, "data/processed_data/Diversity_phenology_1m2.csv")



# Functional composition & diversity ----------------------------------------------------------

## species data ----
species_data <- read_csv("data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv") %>% 
  filter(!(Taxon %in% c("Asteracea", "Plantae (rosette)", "Rosaceae"))) %>% # no any traits
  unite(Plot, PlotNo, Subplot, Month, remove = TRUE) %>% 
  select(Plot, Taxon, cover) %>%
  arrange(Taxon) %>%
  pivot_wider(names_from = Taxon, values_from = cover, values_fill = 0) %>% 
  column_to_rownames("Plot")

str(species_data)

# check if there are species with zero appearance
which(rowSums(species_data) == 0)
which(colSums(species_data) == 0)


## functional groups data ----
Func_groups <- read_csv("data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv") %>% 
  filter(!(Taxon %in% c("Asteracea", "Plantae (rosette)", "Rosaceae"))) %>% # no traits for it
  select(Taxon) %>% distinct(Taxon) %>% arrange(Taxon) %>% 
  left_join(read_csv("data/processed_data/Functional_composition.csv") %>% 
              select(-EuroMed, -status_detailed),
            by="Taxon") %>% 
  column_to_rownames("Taxon") 

str(Func_groups)
names(Func_groups)

# any species with NA across all traits?
Func_groups %>% filter(if_all(everything(), is.na))


## 1) CWM for functional groups ----
dim(species_data)
dim(Func_groups)

FuncComp <- FD::functcomp(Func_groups %>% 
                            rename(FG=functional_group,
                                   lifespan_shortlived=lifespan_biennialOrShortLived,
                                   lifeform_tree.shrub=lifeform_tree_schrub),  
                          as.matrix(species_data), 
                          CWM.type = "all", 
                          bin.num=c(  #  bin.num - indicates binary traits to be treated as continuous  
                            "raunkiaer_phanerophyte", "lifeform_tree.shrub",
                            "lifeform_herbPoli", "lifeform_herbMono")) 

FuncComp

dim(FuncComp)

write_csv(as.data.frame(FuncComp) %>% 
            rownames_to_column("Plot") %>% 
            separate(Plot, into = c("PlotNo", "Subplot", "Month"), sep = "_",
                     remove = TRUE),
          "data/processed_data/CWM_FunctCompos_1m2.csv")


## 2) Functional diversity ----

# When calculating functional diversity, FD::dbFD Error occurs when mixed categorical and numeric traits
# thus, convert categorical traits to wide format for each category
Func_groups_modif <- Func_groups %>% 
  select(-starts_with("raunkiaer_"), # correlate strongly with lifespan
         -starts_with("lifeform_"),  # correlate strongly with lifespan
         # exclude all disturb indic for calculating functional diversity
         -Disturbance.Frequency,-Disturbance.Severity,  
         -Mowing.Frequency, -Grazing.Pressure, -Soil.Disturbance) %>%
  rownames_to_column("Taxon") %>%
  # for column status, create wide format 
  mutate(present = 1) %>%
  pivot_wider(names_from = status, values_from=present, values_fill = 0,
              names_prefix = "status_") %>%
  mutate(across(starts_with("status_"), ~ replace(.x, status_NA==1, NA))) %>%
  select(-status_NA) %>% 
  # for column functional_group, create wide format 
  mutate(present = 1) %>%
  pivot_wider(names_from = functional_group, values_from=present, values_fill = 0,
              names_prefix = "FunGr_") %>%
  column_to_rownames("Taxon")

str(Func_groups_modif)

FuncDiv <- FD::dbFD(Func_groups_modif, as.matrix(species_data), 
                    calc.CWM = F,
                    corr = c("cailliez")) 
FuncDiv

str(FuncDiv)

FuncDiv <- as.data.frame(FuncDiv) %>% 
  rownames_to_column("Plot") %>%
  select(-nbsp, -qual.FRic, -sing.sp) %>% 
  separate(Plot, into = c("PlotNo", "Subplot", "Month"), sep = "_",
           remove = TRUE)


# save Functional diversity indices
write_csv(FuncDiv, "data/processed_data/Functional_Diversity_1m2.csv")

