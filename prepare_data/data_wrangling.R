# data wrangling


library(tidyverse)


# data ---------------

# mowing data
Mowing_data <- read_csv("data/mowing_events_2025.csv") %>% 
  pivot_longer(cols = c(September,	July,	May,	March),
               names_to = "Month",
               values_to = "n_mow_events_befre_sampling") 

# litter and other cover data
Cover_data <- read_csv("data/BC_2025_Cover_Data.csv") %>% 
  filter(Scale_m2 == 1) %>%
  select(-Date, -Scale_m2, -Veg_Total_Cover, -"10m_Max_Cryptogam_Height", -Remarks)

# traits data

traits <- read_csv("data/traits.csv")

# vegetation data
dat <- read_csv("data/Sampling5.0Data.csv") %>%
  rename("cover" = `1.00`,
         EuroMed = "Euro+Med Taxon") %>%
  mutate(Layer=ifelse(Layer=="S", "Seedling", Layer),
         cover=ifelse(cover=="X", 0.001, cover)) %>% 
  # X in cover columns means that we forgot to give the cover during the sample 
  # to the species that were present on the plot. We replace "X" with a very low cover value to account for the species in species richness
  mutate(Date = as.Date(Date, "%d/%m/%Y"),
         Juvenile = as.numeric(Juvenile),
         PostFruiting = as.numeric(PostFruiting),
         cover=as.numeric(cover)) %>%
  filter(Layer!="B", Layer!="L", # remove bryophytes & lichens
         Subplot!="EXT") %>%  # EXT means that we found species in the big 10^2m plot that are not inside the 1^2m nested subplots. Such species always have 0% cover in NW and SE corners
 # filter(!Subplot=="EXT") %>%  # remove data for 10m2
  mutate(height=ifelse(height=="X", NA, height)) %>%
  mutate(height = vapply(strsplit(height, ",\\s*"),
                             function(vals) {
                               nums <- as.numeric(vals)
                               nums <- nums[!is.na(nums)]
                               if (length(nums) == 0) NA_real_ else mean(nums)
                             }, numeric(1))) %>% 
  select(PlotNo, Subplot, Date, Layer, Taxon, EuroMed, cover, height, phen,	Seedling, Juvenile, 
         FlowerBud, Flowering, Fruiting, PostFruiting) %>% 
  mutate(across(c(Seedling, Juvenile, FlowerBud, Flowering, Fruiting, PostFruiting),
                ~ replace(., is.na(.), 0))) %>% 
  mutate(Month = lubridate::month(Date, label = TRUE, abbr = FALSE),
         .after =Date, .keep = "unused") %>% 
  mutate(Month = case_when(Month == "April" ~ "March",
                               Month == "August" ~ "July",
                               Month == "October" ~ "September",
                               .default = Month)) %>% 
 # left_join(Mowing_data, by = c("PlotNo", "Subplot", "Month")) %>%  
 # left_join(Cover_data, by = c("PlotNo", "Subplot", "Month")) %>%
  mutate(Biomass = height * cover, .after=height) %>% 
  mutate(Month = factor(Month, levels = c("March", "May", "July", "September"))) 


# dplyr::last_dplyr_warnings() # NAs are introduced where empty cells

str(dat)

# write_csv(traits, "data/processed_data/traits.csv")

# check data  -----
dat %>% pull(Layer) %>% unique()

dat %>% filter(is.na(cover))

dat %>% select(Taxon,  EuroMed) %>% 
  group_by(Taxon,  EuroMed) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  print(n=Inf) %>% 
  left_join(traits, by = c("Taxon", "EuroMed")) %>% 
  filter(is.na(lifeform_lianaHerb))

# Community composition --------


CommunCompos <- dat %>% 
  group_by(PlotNo, Month, Taxon) %>%
  summarise(cover=sum(cover, na.rm = TRUE),  # or mean?
            .groups = "drop") %>%  
  mutate(cover=ifelse(is.na(cover), 0.01, cover))

write_csv(CommunCompos, "data/processed_data/Community_Composition.csv")  

?summarise

