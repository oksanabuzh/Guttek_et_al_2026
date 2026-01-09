# Purpose: Calculate community diversity and phenology metrics at 1m2 plot level


library(tidyverse)

dat <- read_csv("data/processed_data/vegetation2025_1m2.csv")

traits <- read_csv("data/processed_data/Traits_Dist.Ind.Values.csv")
names(traits)

# Community composition: species & phenological ----------------------------------------------------------
Composition_1m2 <- dat %>%
   left_join(traits %>% select(-EuroMed, -status_detailed),
             by="Taxon") %>%
  # calculate Biomass-weighted phenology for each species for each plot and month
  mutate(Neophyte=ifelse(status=="Neophyte", 1, 0),
         Endangered=ifelse(status=="Endangered", 1, 0))%>%
  mutate(
    Seedling = Biomass * Seedling,
    Juvenile = Biomass * Juvenile,
    Vegetative = Biomass * Vegetative,
    FlowerBud     = Biomass * FlowerBud,
    Flowering     = Biomass * Flowering,
    Fruiting = Biomass * Fruiting,
    PostFruiting  = Biomass * PostFruiting,
    Annuals = Biomass * lifespan_annual,
    ShortLived = Biomass * lifespan_biennialOrShortLived,
    Perennials = Biomass * lifespan_perennial,
    Neophytes = Biomass * Neophyte,
    Endangereds = Biomass * Endangered) 

Composition_1m2

Composition_1m2 %>% 
  filter(is.na(Biomass))

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
            Seedling_mass = sum(Seedling, na.rm = TRUE),
            Juvenile_mass = sum(Juvenile, na.rm = TRUE),
            Vegetative_mass = sum(Vegetative, na.rm = TRUE),
            FlowerBud_mass = sum(FlowerBud, na.rm = TRUE),
            Flowering_mass = sum(Flowering, na.rm = TRUE),
            Fruiting_mass = sum(Fruiting, na.rm = TRUE),
            PostFruiting_mass = sum(PostFruiting, na.rm = TRUE),
            Annuals_mass = sum(Annuals, na.rm = TRUE),
            ShortLived_mass = sum(ShortLived, na.rm = TRUE),
            Perennials_mass = sum(Perennials, na.rm = TRUE),
            Neophytes_mass = sum(Neophytes, na.rm = TRUE),
            Endangereds_mass = sum(Endangereds, na.rm = TRUE),
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
 # mutate(biomass=ifelse(is.na(height_mean), NA, biomass), # if height is NA, biomass should be NA and not 0
 #        height_max=ifelse(is.na(height_mean), NA, height_max)) %>%  # if height is NA, height_max should be NA and not Inf
  mutate(Seedling_mass_propr = Seedling_mass / biomass,
         Juvenile_mass_propr = Juvenile_mass / biomass,
         FlowerBud_mass_propr = FlowerBud_mass / biomass,
         Flowering_mass_propr = Flowering_mass / biomass,
         Fruiting_mass_propr = Fruiting_mass / biomass,
         PostFruiting_mass_propr = PostFruiting_mass / biomass,
         Annuals_mass_propr = Annuals_mass / biomass,
         ShortLived_mass_propr = ShortLived_mass / biomass,
         Perennials_mass_propr = Perennials_mass / biomass,
         Neophytes_mass_propr = Neophytes_mass / biomass,
         Endangereds_mass_propr = Endangereds_mass / biomass,
         Seedling_SR_propr = Seedling_SR / SR,
         Juvenile_SR_propr = Juvenile_SR / SR,
         FlowerBud_SR_propr = FlowerBud_SR / SR,
         Flowering_SR_propr = Flowering_SR / SR,
         Fruiting_SR_propr = Fruiting_SR / SR,
         PostFruiting_SR_propr = PostFruiting_SR / SR,
         Annuals_SR_propr = Annuals_SR / SR,
         ShortLived_SR_propr = ShortLived_SR / SR,
         Perennials_SR_propr = Perennials_SR / SR,
         Neophytes_SR_propr = Neophytes_SR / SR,
         Endangereds_SR_propr = Endangereds_SR / SR) %>% 
  # Phenological diversity 
  rowwise() %>%
  mutate(
    phen_Richness = sum(c_across(c(Vegetative_mass,
                                   Seedling_mass,	Juvenile_mass,FlowerBud_mass,
                                   Flowering_mass,	Fruiting_mass,	PostFruiting_mass)) > 0),
    
    phen_Shannon = exp(vegan::diversity(
      c_across(c(Vegetative_mass, Seedling_mass,	Juvenile_mass,
                 FlowerBud_mass, Flowering_mass,	Fruiting_mass,	
                 PostFruiting_mass)), index = "shannon")),
    phen_evenness = vegan::diversity(
      c_across(c(Vegetative_mass, Seedling_mass,	Juvenile_mass,
                 FlowerBud_mass, Flowering_mass,	Fruiting_mass,	
                 PostFruiting_mass)), index = "invsimpson"),
    .after = "height_max") %>% 
  ungroup()

Diversity_phenology_1m2


#  Warning message concerns height_max, it get -Inf as no measurements were present for the plot


write_csv(Diversity_phenology_1m2, "data/processed_data/Diversity_phenology_1m2.csv")



# Functional composition & diversity ----------------------------------------------------------


Functtraits <- read_csv("data/processed_data/Functional_composition.csv") %>% 
  select(-EuroMed, -status_detailed) 

# any species with NA across all traits?
Functtraits %>% 
  column_to_rownames("Taxon") %>%
  filter(if_all(everything(), is.na))
# three species with all traits=NA
# c("Asteracea", "Plantae (rosette)", "Rosaceae")

## species data ----

species_data <- read_csv("data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv") %>% 
  filter(!(Taxon %in% c("Asteracea", "Plantae (rosette)", "Rosaceae"))) %>% # species with all traits=NA
  unite(Plot, PlotNo, Subplot, Month, remove = TRUE) %>% 
  select(Plot, Taxon, Biomass) %>%
  arrange(Taxon) %>%
  pivot_wider(names_from = Taxon, values_from = Biomass, values_fill = 0) %>% 
  column_to_rownames("Plot")

str(species_data)

# check if there are species with zero appearance
which(rowSums(species_data) == 0)
which(colSums(species_data) == 0)


## functional groups data ----
Func_groups <- read_csv("data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv") %>% 
  filter(!(Taxon %in% c("Asteracea", "Plantae (rosette)", "Rosaceae"))) %>% # species with all traits=NA
  select(Taxon) %>% distinct(Taxon) %>% arrange(Taxon) %>%  # to ensure we use same species list and order as in species_data
  left_join(read_csv("data/processed_data/Functional_composition.csv") %>% 
              select(-EuroMed, -status_detailed),
            by="Taxon") %>% 
  column_to_rownames("Taxon") 

str(Func_groups)
names(Func_groups)

# check again if any species with NA across all traits?
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
  # for column status, create wide a format 
  mutate(present = 1) %>%
  pivot_wider(names_from = status, values_from=present, values_fill = 0,
              names_prefix = "status_") %>%
  mutate(across(starts_with("status_"), ~ replace(.x, status_NA==1, NA))) %>%
  select(-status_NA) %>% 
  # for column functional_group, create a wide format 
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

