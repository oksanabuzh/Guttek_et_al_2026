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
  mutate(height=ifelse(height=="X", NA, height)) %>%
  mutate(height = vapply(strsplit(height, ",\\s*"),
                             function(vals) {
                               nums <- as.numeric(vals)
                               nums <- nums[!is.na(nums)]
                               if (length(nums) == 0) NA_real_ else mean(nums)
                             }, numeric(1))) %>% 
  select(PlotNo, Subplot, Date, Taxon, EuroMed, height, phen,	Seedling, Juvenile, 
         FlowerBud, Flowering, Fruiting, PostFruiting) %>% 
  mutate(across(c(phen, Seedling, Juvenile, FlowerBud, Flowering, Fruiting, PostFruiting),
                ~ replace(., is.na(.), 0))) %>% 
  mutate(Month = lubridate::month(Date, label = TRUE, abbr = FALSE))


str(dat)


dat$Month


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
