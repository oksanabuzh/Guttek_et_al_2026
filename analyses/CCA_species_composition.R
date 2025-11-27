# Purpose:  Canonical correspondence analysis for species composition sampled at 1m2 plots

dev.off()


library(tidyverse)
library(vegan)
library(ggrepel)

citation()

# read & wrangle data ----- 
species_data <- read_csv("data/Compos1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, MonthLabel, remove = FALSE) %>% 
  column_to_rownames("Plot") %>% 
  dplyr::select(-c( "PlotNo", "MowFreq", "MonthLabel", "Subplot"))

str(species_data)

plot_data <- read_csv("data/Compos1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, MonthLabel, remove = FALSE) %>%  
  unite(Random_effect, PlotNo, MonthLabel, remove = FALSE) %>% 
  dplyr::select(Plot, Random_effect, PlotNo, MowFreq, MonthLabel, Subplot) 

str(plot_data)

Trait_data <- read_csv("data/traits.csv") %>% 
  rename(species=Taxon) %>% 
  mutate(trait=
           case_when(trait=="NotThreatened"~ "Not threatened",
                     trait=="ForeWarning"~ "Threatened",
                     is.na(trait) ~ "Data insufficient / unknown",
                     .default=trait))  

Trait_data %>% 
  distinct(trait)

# Data exploration -----
# check both data sets for NA‘s
anyNA(species_data) # no NA's
anyNA(plot_data) # no NA's

## Linear or nonlinear methods to use? ----
# check gradient length of first DCA axis (optional)
# if axis lengths for DCA1 is 
# <3 -> linear methods (PCA)
# >3 -> nonlinear methods (CCA)
# in any case non metric distance based methods can be used (NMDS or PCoA)
decorana(species_data) 
#  >3 -> linear methods  are not applicable as axis lengths for DCA1 is >3

# Canonical correspondence analysis
set.seed(1)
ord_mod1 <- cca(species_data ~ MowFreq:MonthLabel + MowFreq + MonthLabel, data = plot_data,
             scale = FALSE) # scale data to have the same units
ord_mod1
anova(ord_mod1, strata = as.factor(plot_data$PlotNo), # random effects
                                   by= "terms") # each term (sequentially from first to last), depends on the order

set.seed(1)
ord_mod <- cca(species_data ~ MowFreq + MonthLabel , data = plot_data,
             scale = FALSE) # scale data to have the same units
ord_mod

anova(ord_mod, strata = as.factor(plot_data$PlotNo),
      by= "margin") # test for conditional effects (does not depend on the order)


# multicolinearity check
vif.cca(ord_mod)
# proportion variance explained by RDA axes
summary(eigenvals(ord_mod))
# adjusted R2
RsquareAdj(ord_mod)

# Permutation tests ------
## --- model fit statistics ---
set.seed(1)
anova(ord_mod, strata = as.factor(plot_data$PlotNo)) # model fit statistics

plot(ord_mod)

# extract species scores
sp.scrs <- scores(ord_mod, display = "species",
                  scaling = "species") %>% 
  as_tibble(rownames = "species") %>% 
  left_join(Trait_data, by="species") %>% 
  mutate(species_full_name=species,
         species = if_else(
           str_count(species, "\\S+") == 1,      # If only one word (non-space sequence)
           paste(species, "sp."),                # add "sp."
           species                               # else keep as is
         ),
         species = str_c(str_split_i(species, '\\s+', 1) %>%    # splits the species_name at each empty space in the species name and extracts the first word (the genus)
                           str_sub(.,  1, 5 ),          #  in this string (".") subtracts first 4 letters of genus (start, end 
                         str_split_i( species, '\\s+', 2) %>%   # gets the second part of the species name after the first empty space (species)
                           str_sub(., 1, 3),            #  subtracts first 5 letters of from that second part (species)
                         sep = '.')) %>% 
  mutate(trait=fct_relevel(trait,"Threatened",
                           "Not threatened", 
                           "Neophyte", 
                           "Data insufficient / unknown"))

sp.scrs

sp.scrs %>% 
  filter(is.na(species))

# extract plot scores 
plot.scrs <- scores(ord_mod, display = "sites",
                    scaling = "sites") %>% 
  as_tibble(rownames = "Plot") %>% 
  left_join(plot_data, by="Plot") 

plot.scrs

# calculate centroid for  Grazing_season
centroid_mowing <- scores(ord_mod, 
                    display="cn",  
                    scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "MowFreq")) %>% 
  mutate(MowFreq=stringr::str_sub(treatment, 8)) %>% 
  dplyr::select(-treatment) %>% 
  rename(CCA1_mowing=CCA1,
         CCA2_mowing=CCA2)

centroid_mowing

centroid_month <- scores(ord_mod, 
                         display="cn",  
                         scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "MonthLabel")) %>% 
  mutate(MonthLabel=stringr::str_sub(treatment, 11)) %>% 
  dplyr::select(-treatment) %>% 
  rename(CCA1_month=CCA1,
         CCA2_month=CCA2)

centroid_month

# centroid for interaction from raw data
centroids <- plot.scrs %>% 
  group_by(MowFreq, MonthLabel) %>% 
  summarise(CCA1_centroid=mean(CCA1),
            CCA2_centroid=mean(CCA2)) %>% 
  ungroup() %>% 
  left_join(centroid_mowing, by=c("MowFreq")) %>% 
  left_join(centroid_month, by=c("MonthLabel")) %>%
  mutate(Mowing=case_when(
    MowFreq == "reduced_sown" ~ "reduced mowing & sowing",
    MowFreq == "regular" ~ "regular mowing",
    MowFreq == "reduced" ~ "reduced mowing",
    TRUE ~ as.character(MowFreq)))

centroids

# merge with site scores, order levels of categorical predictors
plot.scrs <- plot.scrs %>%
   left_join(centroids, by=c("MowFreq", "MonthLabel")) %>%
    mutate(Mowing=fct_relevel(Mowing,"regular mowing", "reduced mowing", "reduced mowing & sowing")) %>% 
   mutate(MonthLabel=fct_relevel(MonthLabel,"March", "May", "July", "September")) 

plot.scrs

set.seed(11)
# plot for plots data
plot1 <- ggplot(data=plot.scrs, 
                aes(x=CCA1, y=CCA2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # spiders

  geom_segment(data = plot.scrs,        
               mapping = aes(xend = CCA1_centroid, yend = CCA2_centroid, 
                             color=Mowing),
               alpha=0.5) +
  # add plot scores as point:
  geom_point(data=plot.scrs, 
             aes(x=CCA1, y=CCA2, 
                 color=Mowing),
             size=1.5, pch=21) + 
  # add centroids as point:
  geom_point(data=plot.scrs, 
             aes(x=CCA1_centroid, y=CCA2_centroid, 
                 color=Mowing),
             size=3, pch=18) + 
  # centroids as text
  geom_text_repel(data=centroids, 
                  #geom_text(data=centroids, 
                  aes(x=CCA1_centroid, y=CCA2_centroid, 
                      color=Mowing, label = MonthLabel), 
                  size=5, fontface="bold", show_guide = F) +
  theme_bw()+
  scale_color_manual(values = c("#F8766D", "#00B0F6","#00BA38"))+
  labs(color="Management", x="CCA1 (3.6 %)", y="CCA2 (3.2 %)")


print(plot1)


# ggsave("CCA_plot1.png", plot1, width = 6, height = 6, dpi = 350)
# ggsave("CCA_plot1.jpeg", plot1, width = 6, height = 6, dpi = 350)

# plot for species data
set.seed(11)
plot2 <- ggplot(data=plot.scrs, 
                aes(x=CCA1, y=CCA2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # ellipse 
  stat_ellipse(aes(fill=Mowing), alpha=0.1,
               type='t', # type = 't' means the ellipses are calculated assuming a multivariate t-distribution (robust to outliers)
               linewidth =0.0001, geom="polygon",
               level=0.95, # ellipses represent a 95% confidence interval for the multivariate mean of each group) +
               color="gray88") +
  # species
  geom_point(data=sp.scrs, 
             aes(x=CCA1, y=CCA2, color=trait), 
             size = 0.5,  
             alpha=0.8, pch=19)+
  geom_text_repel(data=sp.scrs, 
                  aes(x=CCA1, y=CCA2, label = species,
                      color=trait), #Genus , Sociality ,  life_strategy, Nest_location
                  size=3, fontface="bold", show_guide = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
  scale_color_manual(values = c(
    "Not threatened" = "#1b9e77",   
   # "Endangered"    = "#e41a1c",   
    "Threatened"     = "#e41a1c",   
    "Neophyte"       = "#3b4cc0",   
    "Data insufficient / unknown"  = "#999999")) +
  scale_fill_manual(values = c(
    "regular mowing" = "#F8766D",
    "reduced mowing" = "yellow3",
    "reduced mowing & sowing" = "#00B0F6"
  )) +
  labs(color="Red list species and neophytes", fill="Management",
       x="CCA1 (3.6 %)", y="CCA2 (3.2 %)")

print(plot2)
