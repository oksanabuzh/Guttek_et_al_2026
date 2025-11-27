# Purpose: to run ordination for species composition sampled at 1m2 plots

library(tidyverse)
library(vegan)
library(ggrepel)

citation()

# read & wrangle data ----- 
species_data <- read_csv("data/Compos1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, MonthLabel, remove = FALSE) %>% 
  column_to_rownames("Plot") %>% 
  select(-c( "PlotNo", "MowFreq", "MonthLabel", "Subplot"))

str(species_data)

plot_data <- read_csv("data/Compos1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, MonthLabel, remove = FALSE) %>%  
  unite(Random_effect, PlotNo, MonthLabel, remove = FALSE) %>% 
  dplyr::select(Plot, Random_effect, PlotNo, MowFreq, MonthLabel, Subplot) %>% 
  mutate(Random_effect=as.factor(Random_effect))

str(plot_data)

Trait_data <- read_csv("data/traits.csv") %>% 
  rename(species=Taxon) %>% 
  mutate(trait=
           case_when(trait=="NotThreatened"~ "Not threatened",
                     trait=="ForeWarning"~ "Endangered",
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
#  >3 -> linear methods methods are not applicable as axis lengths for DCA1 is >3

# anyway we do NMDS 

## Graphical data exploration -----
species_data %>%
  pivot_longer(everything(), names_to = "Species", 
               values_to = "Abundance") %>% 
  ggplot(aes(x = Abundance, y = Species)) +
  geom_boxplot() 
# we see large differences in abundances
# therefore the ordination can be dominated by dominant taxa

# wisconsin transformation in NMDS  removes the influence of dominant abundance, so that dominant species don't dominate the ordination.
species_data %>%
  wisconsin() %>% 
  pivot_longer(everything(), names_to = "Species", 
               values_to = "Abundance") %>% 
  ggplot(aes(x = Abundance, y = Species)) +
  geom_boxplot() 


# check also plots:
species_data %>% 
  rownames_to_column("plot_ID") %>% 
  pivot_longer(- plot_ID, values_to = "abund", names_to = "species") %>% 
  group_by( plot_ID) %>% 
  summarise(sum=sum(abund))%>% 
  ggplot(aes(x = sum, y = plot_ID)) +
  geom_bar(stat = "identity") 

# wisconsin transformation in NMDS  removes the influence of overall abundance at a plot, so that sites with higher total species counts don't dominate the ordination.

species_data %>% 
  wisconsin() %>% 
  rownames_to_column("plot_ID") %>%  
  pivot_longer(- plot_ID, values_to = "abund", names_to = "species") %>% 
  group_by( plot_ID) %>% 
  summarise(sum=sum(abund))%>% 
  ggplot(aes(x = sum, y = plot_ID)) +
  geom_bar(stat = "identity") 



set.seed(2435)
ord_mod <- metaMDS(wisconsin(species_data), 
                   scale = FALSE, distance = "bray") 

ord_mod

# NMDS fit
ord_mod$stress
# fit
stressplot(ord_mod, main = "Shepard plot")


# Permutation test:
set.seed(10)
PERM1 <- vegan::adonis2(species_data ~ 
                         MowFreq * MonthLabel, 
                            data=plot_data,
                            permutations = 1000, method = "bray",
                       strata=as.factor(plot_data$PlotNo),
                       by = "terms")

PERM1

PERM2 <- vegan::adonis2(species_data ~ 
                          MowFreq + MonthLabel, 
                        data=plot_data,
                        permutations = 1000, method = "bray",
                        strata=as.factor(plot_data$PlotNo),
                        by = "terms")

PERM2

# variable fitting for posthoc plotting
set.seed(1259)
fit2 <- vegan::envfit(ord_mod   ~  
                        MowFreq + MonthLabel, 
                      data=plot_data,
                    #  strata=as.factor(plot_data$PlotNo),
                      perm=1000) #


fit2



# exploratory plot
plot(ord_mod, main = "NMDS plot", display = "sites")
plot(ord_mod, main = "NMDS plot", display = "species")
plot(ord_mod, main = "NMDS plot")
plot(fit2)

### Plotting NMDS results using the ggplot --------------------------------------

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
  mutate(trait=fct_relevel(trait,"Endangered", "Threatened",
                           "Not threatened", 
                           "Neophyte", 
                           "Data insufficient / unknown"))

sp.scrs %>% 
  filter(is.na(species))

# extract plot cores 
plot.scrs <- scores(ord_mod, display = "sites",
                    scaling = "sites") %>% 
  as_tibble(rownames = "Plot") %>% 
  left_join(plot_data, by="Plot") 

plot.scrs


# calculate centroid for  Grazing_season
centroid_mowing <- scores(fit2, 
                          display="cn",  
                          scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "MowFreq")) %>% 
  mutate(MowFreq=stringr::str_sub(treatment, 8)) %>% 
  dplyr::select(-treatment) %>% 
  rename(NMDS1_mowing=NMDS1,
         NMDS2_mowing=NMDS2)

centroid_mowing

centroid_month <- scores(fit2, 
                         display="cn",  
                         scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "MonthLabel")) %>% 
  mutate(MonthLabel=stringr::str_sub(treatment, 11)) %>% 
  dplyr::select(-treatment) %>% 
  rename(NMDS1_month=NMDS1,
         NMDS2_month=NMDS2)

centroid_month

# centroid for interaction from raw data
centroids <- plot.scrs %>% 
  group_by(MowFreq, MonthLabel) %>% 
  summarise(NMDS1_centroid=mean(NMDS1),
            NMDS2_centroid=mean(NMDS2)) %>% 
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
  mutate(MonthLabel=fct_relevel(MonthLabel,"March", "May", "July")) 

plot.scrs


# plot for plots data
plot1 <- ggplot(data=plot.scrs, 
                aes(x=NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # spiders
  
  geom_segment(data = plot.scrs,        
               mapping = aes(xend = NMDS1_mowing, yend = NMDS2_mowing, 
                             color=Mowing),
               alpha=0.7) +
  # add plot scores as point:
  geom_point(data=plot.scrs, 
             aes(x=NMDS1, y=NMDS2, 
                 color=Mowing),
             size=1.5, pch=21) + 
  # add centroids as point:
  geom_point(data=plot.scrs, 
             aes(x=NMDS1_centroid, y=NMDS2_centroid, 
                 color=Mowing),
             size=3, pch=18) + 
  # centroids as text
  geom_text_repel(data=centroids, 
                  #geom_text(data=centroids, 
                  aes(x=NMDS1_centroid, y=NMDS2_centroid, 
                      color=Mowing, label = MonthLabel), 
                  size=5, fontface="bold", show_guide = F) +
  theme_bw()+
  scale_color_manual(values = c("#F8766D", "#00B0F6","#00BA38"))+
  labs(color="Management")


print(plot1)


# ggsave("NMDS_plot1.png", plot1, width = 6, height = 6, dpi = 350)
# ggsave("NMDS_plot1.jpeg", plot1, width = 6, height = 6, dpi = 350)

# plot for species data
set.seed(11)
plot2 <- ggplot(data=plot.scrs, 
                aes(x=NMDS1, y=NMDS2))+
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
             aes(x=NMDS1, y=NMDS2, color=trait), 
             size = 0.5,  
             alpha=0.8, pch=19)+
  geom_text_repel(data=sp.scrs, 
                  aes(x=NMDS1, y=NMDS2, label = species,
                      color=trait), #Genus , Sociality ,  life_strategy, Nest_location
                  size=3, fontface="bold", show_guide = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
  scale_color_manual(values = c(
    "Not threatened" = "#1b9e77",   
    "Endangered"    = "#e41a1c",   
    "Threatened"     = "#ffb000",   
    "Neophyte"       = "#3b4cc0",   
    "Data insufficient / unknown"  = "#999999")) +
  scale_fill_manual(values = c(
    "regular mowing" = "#F8766D",
    "reduced mowing" = "yellow3",
    "reduced mowing & sowing" = "#00B0F6"
  )) +
  labs(color="Red list species and neophytes", fill="Management")

plot2
