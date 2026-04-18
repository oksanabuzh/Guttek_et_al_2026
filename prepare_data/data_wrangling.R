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
Mowing_data <- read_csv("data/raw_data/mowing_events_2025_DB.csv") %>% 
  # mowing_events_2025 - total number of mowing events in 2025 (including after September)
  select(-mowing_events_2025) %>%
  pivot_longer(cols = c(September,	July,	May,	March),
               names_to = "Month",
               values_to = "n_mow_events_befre_sampling") 

# litter and other cover data
Cover_data <- read_csv("data/raw_data/BC_2025_Cover_Data.csv") %>% 
  filter(Scale_m2 == 1) %>%
  select(-Date, -Scale_m2, -Veg_Total_Cover, -"10m_Max_Cryptogam_Height", -Remarks)

names(Cover_data)

# prepare phenology data -------------
# phenology and height were measured only on 10m2 scale
data_10m2 <- read_csv("data/raw_data/Sampling7.0Data.csv") %>%
 # mutate(Date = as.Date(Date, "%d/%m/%Y")) %>% 
  mutate(Month = lubridate::month(Date, label = TRUE, abbr = FALSE),
         .after =Date, .keep = "unused") %>% 
  mutate(Month = case_when(Month == "April" ~ "March",
                           Month == "August" ~ "July",
                           Month == "October" ~ "September",
                           .default = Month)) %>%
  rename("cover10m2" = "3.16") %>%     # 1m2 scale
  filter(Layer!="B", Layer!="L" ) %>% # , # remove bryophytes & lichens
  mutate(height=ifelse(height=="X", NA, height)) %>%
  mutate(height = vapply(strsplit(height, ",\\s*"),
                         function(vals) {
                           nums <- as.numeric(vals)
                           nums <- nums[!is.na(nums)]
                           if (length(nums) == 0) NA_real_ else mean(nums)
                         }, numeric(1))) %>% 
  select(PlotNo, Month, Taxon, cover10m2, height, phen,	Seedling, Juvenile, 
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
  mutate(Vegetative=ifelse(
    (Seedling + Juvenile + FlowerBud + Flowering + Fruiting + PostFruiting) == 0, 
    1, 0)) 



str(data_10m2)

data_10m2 %>% 
  filter(is.na(cover10m2)) %>% 
  print(n=Inf)

data_10m2 %>% 
  filter(is.na(height)) %>% 
  print(n=Inf)

## fill missing height data -----
data_10m2_filled <- data_10m2 %>%
  left_join(Mowing_data %>% 
              summarise(n_mow_events_befre_sampling=mean(n_mow_events_befre_sampling),
                        .by=c("PlotNo", "MowFreq", "Month")) %>% 
              relocate(c("MowFreq", "n_mow_events_befre_sampling"), .after=Month),
            by=c("PlotNo", "Month")) %>%
# Step 1: based on Month, MowFreq, n_mow_events_befre_sampling, Taxon
  group_by(Month, MowFreq, n_mow_events_befre_sampling, Taxon) %>%
  mutate(hight_filled_step1 = ifelse(is.na(height),
                                     mean(height, na.rm = TRUE), height),
         .after=height) %>% 
  ungroup() %>% 
# Step 2: based on Month, MowFreq, Taxon
group_by(Month, MowFreq, Taxon) %>%
  mutate(hight_filled_step2 = ifelse(is.na(hight_filled_step1),
                                     mean(height, na.rm = TRUE), hight_filled_step1),
         .after=hight_filled_step1) %>%
  ungroup() %>% 
# Step 3: based on Month, n_mow_events_befre_sampling, Taxon
  group_by(Month, n_mow_events_befre_sampling, Taxon) %>%
  mutate(hight_filled_step3 = ifelse(is.na(hight_filled_step2),
                                      mean(height, na.rm = TRUE), hight_filled_step2),
           .after=hight_filled_step2) %>%
    ungroup() %>% 

    # Step 4: based on Month, Mowing events, Taxon
group_by(MowFreq, n_mow_events_befre_sampling, Taxon) %>%
  mutate(hight_filled_step4 = ifelse(is.na(hight_filled_step3),
                                     mean(height, na.rm = TRUE), hight_filled_step3),
         .after=hight_filled_step3) %>%
  ungroup() %>% 

# Step 5: based on MowFreq, Taxon
  # at step 5 I tried to group by "Month & Taxon" and by "n_mow_events_befre_sampling & Taxon"but it didn't solve NAs in addition to Step 4 
  group_by(MowFreq, Taxon) %>%
  mutate(hight_filled_step5 = ifelse(is.na(hight_filled_step4),
                                     mean(height, na.rm = TRUE), hight_filled_step4),
         .after=hight_filled_step4) %>%
  ungroup()


data_10m2_filled %>% 
  filter(is.na(hight_filled_step1)) %>% 
  print(n=Inf)

data_10m2_filled %>% 
  filter(is.na(hight_filled_step2)) %>% 
  print(n=Inf)

data_10m2_filled %>% 
  filter(is.na(hight_filled_step3)) %>% 
  print(n=Inf)

data_10m2_filled %>% 
  filter(is.na(hight_filled_step4)) %>% 
  print(n=Inf)


data_10m2_filled %>% 
  filter(is.na(hight_filled_step5)) %>% 
  print(n=Inf)

# vegetation data 1m2 ----------------
data_1m2 <- read_csv("data/raw_data/Sampling7.0Data.csv") %>%
  mutate(Date = as.Date(Date, "%d/%m/%Y")) %>% 
  mutate(Month = lubridate::month(Date, label = TRUE, abbr = FALSE),
         .after =Date, .keep = "unused") %>% 
  mutate(Month = case_when(Month == "April" ~ "March",
                           Month == "August" ~ "July",
                           Month == "October" ~ "September",
                           .default = Month)) %>%
  rename("cover" = "1.00",
         EuroMed = "Euro+Med Taxon") %>%
  mutate(Layer=ifelse(Layer=="S", "Seedling", Layer)) %>% 
  filter(Layer!="B", Layer!="L", # remove bryophytes & lichens
         Subplot!="EXT") %>% #, # EXT means that we found species in the big 10^2m plot that are not inside the 1^2m nested subplots. Such species always have 0% cover in NW and SE corners
  select(PlotNo, Subplot, Month, Layer, Taxon, EuroMed, cover) %>% 
    # merge with phenology and height data
  left_join(data_10m2_filled %>% 
              select(-cover10m2, -height, -hight_filled_step1, 
                     -hight_filled_step2, -hight_filled_step3, -hight_filled_step4) %>% 
              rename("height"="hight_filled_step5"),
            by=c("PlotNo", "Month", "Taxon")) %>%
  mutate(Biomass = height * cover, .after=height) %>% 
  mutate(Month = factor(Month, levels = c("March", "May", "July", "September"))) 

# check for missing cover
data_1m2 %>% 
  filter(is.na(cover)) %>% 
  select(Taxon)

# check for repeated species 
data_1m2 %>% 
  group_by(PlotNo, Subplot, Month) %>% 
  count(Taxon) %>% 
  arrange(desc(n))


# check if all species from 1m2 are present in 10m2
data_10m2 %>% 
  select(PlotNo, Month, Taxon, cover10m2) %>%
  left_join(data_1m2 %>% 
              summarise(cover1m2 = mean(cover, na.rm = TRUE),
                        .by=c("PlotNo", "Month", "Taxon")),
            by=c("PlotNo", "Month", "Taxon")
            ) %>% 
  filter(is.na(cover10m2))

data_1m2 %>% 
  filter(is.na(Biomass))

write_csv(data_1m2, "data/processed_data/vegetation2025_1m2.csv")

# trait data ---------------------------------------------------------------

trait_data <- read_csv("data/raw_data/trait_joined_finCharly.csv")

disturbance_indic <- read_csv("data/raw_data/disturbance_indicator_values_Midolo_et_al_2022.csv") %>%
  select(species_corrected_OB, Disturbance.Severity, Disturbance.Frequency, 
                  Mowing.Frequency, Grazing.Pressure, Soil.Disturbance) %>% 
  rename("EuroMed"=species_corrected_OB)

# some taxa are missing in Midolo_et_al_2022. Below some of the taxa are added manually by Dariia Borovyk:
Missing_taxa <- read_csv("data/raw_data/missing_taxa_for_disturbance_indicators_addedd_DB.csv") %>% 
  select(EuroMed, Disturbance.Severity, Disturbance.Frequency, 
         Mowing.Frequency, Grazing.Pressure, Soil.Disturbance)
names(Missing_taxa)

disturbance_indic_new <- disturbance_indic %>% 
  bind_rows(Missing_taxa) 


traits <- trait_data %>% 
  left_join(disturbance_indic_new, by = "EuroMed") 

# check merged data
traits %>% 
  filter(is.na(Disturbance.Severity) | is.na(Disturbance.Frequency) | 
           is.na(Mowing.Frequency) | is.na(Grazing.Pressure) | 
           is.na(Soil.Disturbance)) %>% 
  print(n=Inf)

write_csv(traits, "data/processed_data/Traits_Dist.Ind.Values.csv")


# GIS data ----------------------------------------------------------------


# Data measured in the field:
Fild_data <- read_csv("data/raw_data/GISdata/BC_Plots_Carolin_field_measurements.csv") %>% 
  # split column "name" into "PlotNo" and  "Subplot"
  separate(name, into = c("PlotNo", "Subplot"), sep = 4, remove = FALSE) %>% 
  dplyr::select(PlotNo, Subplot, ele, xcoord, ycoord, slope_aspe, slope_degr,
                mole_holes, trampling, dist_shrub, dist_tree, 
                shrub_cov_buff20, tree_cov_buff20, dist_trail) 
Fild_data %>% 
  print(n=Inf)

write_csv(Fild_data, "data/processed_data/GIS_data/GIS_field_measurements.csv")

# GIS main data (patch type, area, distances):
GIS_main <- read_csv("data/raw_data/GISdata/BC_Plots_Carolin_GIS_patch_type_area_distances.csv") %>% 
  filter(!str_detect(name, "raw")) %>%
  mutate(PlotNo = str_remove_all(name, "established"), 
         .before=name) %>%
  dplyr::select(-name)

write_csv(GIS_main, "data/processed_data/GIS_data/GISmain_PatchTypeArea_Distances.csv")

# Sky-view data:
Sky_view <- read_csv("data/raw_data/GISdata/BC_Plots_Carolin_Sky_view_factor.csv") 
write_csv(Sky_view, "data/processed_data/GIS_data/Sky_view_factor.csv")



# GIS data in buffers:
GIS_250 <- read_csv("data/raw_data/GISdata/BC_Plots_Carolin_GIS_250m.csv") %>% 
  filter(!str_detect(name, "raw")) %>%
  mutate(PlotNo = str_remove_all(name, "250m_|established"), 
         .before=name) %>%
  dplyr::select(-name) 

names(GIS_250)
write_csv(GIS_250, "data/processed_data/GIS_data/GIS_buffer250m.csv")


GIS_500 <- read_csv("data/raw_data/GISdata/BC_Plots_Carolin_GIS_500m.csv") %>% 
  filter(!str_detect(name, "raw")) %>%
  mutate(PlotNo = str_remove_all(name, "500m_|established"), 
         .before=name) %>%
  dplyr::select(-name) 

names(GIS_500)
write_csv(GIS_500, "data/processed_data/GIS_data/GIS_buffer500m.csv")


GIS_1000 <- read_csv("data/raw_data/GISdata/BC_Plots_Carolin_GIS_1000m.csv") %>%
  filter(!str_detect(name, "raw")) %>%
  mutate(PlotNo = str_remove_all(name, "1000m_|established"), 
         .before=name) %>%
  dplyr::select(-name) 

names(GIS_1000)
write_csv(GIS_1000, "data/processed_data/GIS_data/GIS_buffer1000m.csv")

#  Biotops data:
GIS_biotops <- read_csv("data/raw_data/GISdata/BC_Plots_Carolin_GIS_biotops_FU.BiotopMaps_25m.csv") %>% 
  filter(!str_detect(name, "raw")) %>%
  mutate(PlotNo = str_remove_all(name, "established"), 
         .before=name) %>%
  dplyr::select(-name)

names(GIS_biotops)

Biotop_specific_diversity <- GIS_biotops %>% 
  group_by(PlotNo) %>% 
  summarise(Biotop_richness_specific = n_distinct(biotop_code_specific),
            Biotop_Shannon_specific = vegan::diversity(biotope_specific_percent, index = "shannon"),
            .groups = "drop")

Biotop_coarse_diversity <- GIS_biotops %>% 
  group_by(PlotNo, biotop_code_general, biotop_name_general) %>% 
  summarise(biotope_general_percent = sum(biotope_specific_percent),
            .groups = "drop") %>% 
  summarise(Biotop_richness_coarse = n_distinct(biotop_code_general),
            Biotop_Shannon_coarse = vegan::diversity(biotope_general_percent, index = "shannon"),
            .by=PlotNo)

Biotop_diversity <- Biotop_specific_diversity %>% 
  left_join(Biotop_general_diversity, by = "PlotNo")

write_csv(Biotop_diversity, "data/processed_data/GIS_data/Biotop_Diversity.csv")
