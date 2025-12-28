# data wrangling


library(tidyverse)


# data ---------------

dat <- read_csv("data/Sampling5.0Data.csv") %>%
  mutate(Date = as.Date(Date, "%d/%m/%Y"),
         Juvenile = as.numeric(Juvenile),
         PostFruiting = as.numeric(PostFruiting)) %>%
  rename("abundance" = `1.00`,
         EuroMed = "Euro+Med Taxon") %>%
  filter(Layer!="B") %>%
  select(PlotNo, Subplot, Date, Taxon, EuroMed, phen,	Seedling, Juvenile, 
         FlowerBud, Flowering, Fruiting, PostFruiting, height) %>% 
  mutate(across(c(phen, Seedling, Juvenile, FlowerBud, Flowering, Fruiting, PostFruiting),
                ~ replace(., is.na(.), 0)))






names <- read_csv("data/raw_data/Sampling5.0Data.csv") %>% 
  rename(EuroMed = `Euro+Med Taxon`) %>%  
  select(Taxon,  EuroMed) %>% 
  group_by(Taxon,  EuroMed) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  print(n=Inf)


traits <- read_csv("data/traits.csv") %>% 
  left_join(names, by = "Taxon") %>% 
  select(-count) %>% 
  relocate(EuroMed, .after = Taxon)

traits %>% 
  filter(is.na(EuroMed)) 

# write_csv(traits, "data/traits.csv")
