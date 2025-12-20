# data wtrangling

names <- read_csv("data/raw_data/Sampling4.2Data.csv") %>% 
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
