# Purpose: Data wrangling for vegetation data


library(tidyverse)

# metadata -----------

# file Sampling5.0Data.csv:

# 1) X in cover columns means that we forgot to give the cover during the sample 
# to the species that were present on the plot. 

#  2) for 10m2 in the 3.16 column we gave cover for every species 
# by itself (not mean from two 1m2 columns). When we had a species in both corner 
# we only entered it once and left it blank for the other

# 3) phenology and height were measured for each species only on 10m2 scale

# data ---------------

# mowing data
Mowing_data <- read_csv("data/raw_data/mowing_events_2025.csv") %>% 
  pivot_longer(cols = c(September,	July,	May,	March),
               names_to = "Month",
               values_to = "n_mow_events_befre_sampling") 

# litter and other cover data
Cover_data <- read_csv("data/raw_data/BC_2025_Cover_Data.csv") %>% 
  filter(Scale_m2 == 1) %>%
  select(-Date, -Scale_m2, -Veg_Total_Cover, -"10m_Max_Cryptogam_Height", -Remarks)



# prepare phenology data -------------
# phenology and height were measured only on 10m2 scale
data_10m2 <- read_csv("data/raw_data/Sampling5.0Data.csv") %>%
  mutate(Date = as.Date(Date, "%d/%m/%Y")) %>% 
  mutate(Month = lubridate::month(Date, label = TRUE, abbr = FALSE),
         .after =Date, .keep = "unused") %>% 
  mutate(Month = case_when(Month == "April" ~ "March",
                           Month == "August" ~ "July",
                           Month == "October" ~ "September",
                           .default = Month)) %>%
  rename("cover10m2" = "3.16m") %>% # 10m2 scale
  mutate(cover10m2=ifelse(cover10m2=="X", 0.001, cover10m2)) %>% # X in cover columns means that we forgot to give the cover during the sample 
  # to the species that were present on the plot. We replace "X" with a very low cover value to account for the species in species richness
  mutate(cover10m2=as.numeric(cover10m2),
         Juvenile = as.numeric(Juvenile),
         PostFruiting = as.numeric(PostFruiting)) %>%
  filter(Layer!="B", Layer!="L", # remove bryophytes & lichens
         !is.na(cover10m2)) %>%  # there is one NA in cover that is that species is not present
  mutate(height=ifelse(height=="X", NA, height)) %>%
  mutate(height = vapply(strsplit(height, ",\\s*"),
                         function(vals) {
                           nums <- as.numeric(vals)
                           nums <- nums[!is.na(nums)]
                           if (length(nums) == 0) NA_real_ else mean(nums)
                         }, numeric(1))) %>% 
  select(PlotNo, Subplot, Month, Layer, Taxon, cover10m2, height, phen,	Seedling, Juvenile, 
         FlowerBud, Flowering, Fruiting, PostFruiting) %>% 
  summarise(cover10m2=mean(cover10m2, na.rm = TRUE),  # When species was in both 1m2 corners, it was entered in 10m2 only once (e.g. NE subplot) and left it blank for the other corner
            height=mean(height, na.rm = TRUE),
            Seedling=mean(Seedling, na.rm = TRUE), 
            Juvenile=mean(Juvenile, na.rm = TRUE), 
            FlowerBud=mean(FlowerBud, na.rm = TRUE), 
            Flowering=mean(Flowering, na.rm = TRUE), 
            Fruiting=mean(Fruiting, na.rm = TRUE),
            PostFruiting=mean(PostFruiting, na.rm = TRUE),
            .by=c("PlotNo", "Month", "Taxon"))  %>% 
  mutate(across(c(Seedling, Juvenile, FlowerBud, Flowering, Fruiting, PostFruiting),
                ~ replace(., is.na(.), 0))) %>% 
  mutate(Vegetative=ifelse((Seedling + Juvenile + FlowerBud + Flowering + Fruiting + PostFruiting) > 0, 1, 0))
str(data_10m2)


# vegetation data 1m2 ----------------
data_1m2 <- read_csv("data/raw_data/Sampling5.0Data.csv") %>%
  mutate(Date = as.Date(Date, "%d/%m/%Y")) %>% 
  mutate(Month = lubridate::month(Date, label = TRUE, abbr = FALSE),
         .after =Date, .keep = "unused") %>% 
  mutate(Month = case_when(Month == "April" ~ "March",
                           Month == "August" ~ "July",
                           Month == "October" ~ "September",
                           .default = Month)) %>%
  rename("cover" = "1m",
         EuroMed = "Euro+Med Taxon") %>%
  mutate(Layer=ifelse(Layer=="S", "Seedling", Layer),
         cover=ifelse(cover=="X", 0.001, cover)) %>% # X in cover columns means that we forgot to give the cover during the sample 
                                                     # to the species that were present on the plot. We replace "X" with a very low cover value to account for the species in species richness
  mutate(cover=as.numeric(cover)) %>%
  filter(Layer!="B", Layer!="L", # remove bryophytes & lichens
         Subplot!="EXT", # EXT means that we found species in the big 10^2m plot that are not inside the 1^2m nested subplots. Such species always have 0% cover in NW and SE corners
         !is.na(cover)) %>%  # there is one NA in cover that is that species is not present
  select(PlotNo, Subplot, Month, Layer, Taxon, EuroMed, cover) %>% 
  summarise(cover=sum(cover, na.rm = TRUE),  # sum cover for repeted species to get unique species per plot
            .by=c("PlotNo", "Subplot", "Month", "Taxon", "EuroMed")) %>% 
  # merge with phenology and height data
  left_join(data_10m2 %>% 
              select(-cover10m2),
            by=c("PlotNo", "Month", "Taxon")) %>%
  mutate(Biomass = height * cover, .after=height) %>% 
  mutate(Month = factor(Month, levels = c("March", "May", "July", "September"))) 


str(data_1m2)

# check if all species from 1m2 are present in 10m2
data_10m2 %>% 
  select(PlotNo, Month, Taxon, cover10m2) %>%
  left_join(data_1m2 %>% 
              summarise(cover1m2 = mean(cover, na.rm = TRUE),
                        .by=c("PlotNo", "Month", "Taxon")),
            by=c("PlotNo", "Month", "Taxon")
            ) %>% 
  filter(is.na(cover10m2))


write_csv(data_1m2, "data/processed_data/vegetation2025_1m2.csv")

# trait data ---------------------------------------------------------------

# traits data
trait_data <- read_csv("data/raw_data/traits.csv")

disturbance_indic <- read_csv("data/raw_data/disturbance_indicator_values_Midolo_et_al_2022.csv") %>%
  select(species_corrected_OB, Disturbance.Severity, Disturbance.Frequency, 
                  Mowing.Frequency, Grazing.Pressure, Soil.Disturbance) %>% 
  rename("EuroMed"=species_corrected_OB)

traits <- trait_data %>% 
  left_join(disturbance_indic, by = "EuroMed") 

write_csv(traits, "data/processed_data/Traits_Dist.Ind.Values.csv")

# check merged data
missing_Taxa <- traits%>% 
  select(Taxon, EuroMed, functional_group, Disturbance.Severity, Disturbance.Frequency, 
         Mowing.Frequency, Grazing.Pressure, Soil.Disturbance) %>%
  filter(if_any(all_of(
    c("Disturbance.Severity", "Disturbance.Frequency",
      "Mowing.Frequency", "Grazing.Pressure", "Soil.Disturbance")),
    is.na)) %>% print(n=Inf)

missing_Taxa

# write_csv(missing_Taxa, "data/missing_taxa_for_disturbance_indicators.csv")

