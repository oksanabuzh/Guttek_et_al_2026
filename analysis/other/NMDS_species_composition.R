# Purpose: to run ordination for species composition sampled at 1m2 plots
dev.off()


library(tidyverse)
library(vegan)
library(ggrepel)

citation()

# read & wrangle data ----- 

## species data ----
species_data <- read_csv("data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, Month, remove = TRUE) %>% 
  rename(Taxon_EuroMed=EuroMed) %>%
  select(Plot, Taxon_EuroMed, cover) %>%
  arrange(Taxon_EuroMed) %>%
  pivot_wider(names_from = Taxon_EuroMed, values_from = cover, values_fill = 0) %>% 
  column_to_rownames("Plot")

str(species_data)

## predictor data ----

# mowing data
Mowing_data <- read_csv("data/raw_data/mowing_events_2025_DB.csv") %>% 
  select(-mowing_events_2025) %>% 
  pivot_longer(cols = c(September,	July,	May,	March),
               names_to = "Month",
               values_to = "n_mow_events_befre_sampling") %>% 
  relocate(Month, .after=Subplot)

# litter and other cover data
Cover_data <- read_csv("data/raw_data/BC_2025_Cover_Data.csv") %>% 
  filter(Scale_m2 == 1) %>%
  mutate(Litter_Cover_litte_from_2025_mowing=
           ifelse(is.na(Litter_Cover_litte_from_2025_mowing), 0, 
                  Litter_Cover_litte_from_2025_mowing)) %>%
  select(-Date, -Scale_m2, -Veg_Total_Cover, -"10m_Max_Cryptogam_Height", -Remarks)

str(Cover_data)
summary(Cover_data)

plot_data <- Mowing_data %>% 
  left_join(Cover_data, by=c("PlotNo", "Subplot", "Month")) %>% 
  unite(Plot, PlotNo, Subplot, Month, remove = FALSE) 


predictor_data <- species_data %>% # make same row order as in species data 
  rownames_to_column("Plot") %>% select(Plot) %>%
  left_join(plot_data, by="Plot") 

str(predictor_data)

# species trait (for coloring species in ordination plot)
Trait_data <- read_csv("data/processed_data/Traits_Dist.Ind.Values_fuzzy.csv") %>% 
  rename(Taxon_EuroMed=EuroMed) %>%
  mutate(status=
           case_when(is.na(status) ~ "data insufficient",
                     status == "least_concerned" ~ "least-concerned",
                     status == "neophyte" ~ "neophytes",
                     .default=status))  


Trait_data %>% 
  distinct(status)

# Data exploration -----
# check both data sets for NA‘s
anyNA(species_data) # no NA's
anyNA(plot_data) # no NA's

## Linear or nonlinear methods to use? ----
# check gradient length of first DCA axis (optional)
# if axis lengths for DCA1 is 
# <3 -> linear methods (PCA)
# >3 -> nonlinear methods (NMDS)
# in any case non metric distance based methods can be used (NMDS or PCoA)
decorana(species_data) 
#  >3 -> linear methods  are not applicable as axis lengths for DCA1 is >3
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


# Permutation test:  --------------------------------------
set.seed(10)
PERM1 <- vegan::adonis2(species_data ~ 
                          MowFreq + Month + 
                          n_mow_events_befre_sampling, 
                            data=plot_data,
                            permutations = 1000, method = "bray",
                       strata=as.factor(plot_data$PlotNo),
                       by = "terms")

PERM1



# variable fitting for posthoc plotting  ------------------
set.seed(1259)
fit2 <- vegan::envfit(ord_mod   ~  
                        MowFreq + Month + 
                        n_mow_events_befre_sampling, 
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
  as_tibble(rownames = "Taxon_EuroMed") %>% 
  left_join(Trait_data, by="Taxon_EuroMed") %>% 
  mutate(species_full_name=Taxon_EuroMed,
         Taxon_EuroMed = if_else(
           str_count(Taxon_EuroMed, "\\S+") == 1,      # If only one word (non-space sequence)
           paste(Taxon_EuroMed, "sp."),                # add "sp."
           Taxon_EuroMed                               # else keep as is
         ),
         Taxon_EuroMed = str_c(str_split_i(Taxon_EuroMed, '\\s+', 1) %>%    # splits the species_name at each empty space in the species name and extracts the first word (the genus)
                                 str_sub(.,  1, 5 ),          #  in this string (".") subtracts first 4 letters of genus (start, end 
                               str_split_i( Taxon_EuroMed, '\\s+', 2) %>%   # gets the second part of the species name after the first empty space (species)
                                 str_sub(., 1, 3),            #  subtracts first 5 letters of from that second part (species)
                               sep = '.')) %>% 
  mutate(Taxon_EuroMed = ifelse(Taxon_EuroMed=="Plant.(ro", "Plantae", Taxon_EuroMed)) %>% 
  mutate(trait=fct_relevel(status,"endangered",
                           "vulnerable", 
                           "least-concerned",
                           "neophytes", 
                           "data insufficient"))



sp.scrs %>% 
  pull(status) %>% 
  unique()


# extract plot scores 
plot.scrs <- scores(ord_mod, display = "sites",
                    scaling = "sites") %>% 
  as_tibble(rownames = "Plot") %>% 
  left_join(predictor_data, by="Plot") 

plot.scrs

names(plot.scrs)

# calculate centroid for  Grazing_season
centroid_mowing <- scores(fit1, 
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
  filter(str_detect(treatment, "Month")) %>% 
  mutate(Month=stringr::str_sub(treatment, 6)) %>% 
  dplyr::select(-treatment) %>% 
  rename(NMDS1_month=NMDS1,
         NMDS2_month=NMDS2)

centroid_month

# centroid for interaction from raw data
centroids <- plot.scrs %>% 
  group_by(MowFreq, Month) %>% 
  summarise(NMDS1_centroid=mean(NMDS1),
            NMDS2_centroid=mean(NMDS2)) %>% 
  ungroup() %>% 
  left_join(centroid_mowing, by=c("MowFreq")) %>% 
  left_join(centroid_month, by=c("Month")) %>%
  mutate(Mowing=case_when(
    MowFreq == "reduced_sown" ~ "reduced mowing & sowing",
    MowFreq == "regular" ~ "regular mowing",
    MowFreq == "reduced" ~ "reduced mowing",
    TRUE ~ as.character(MowFreq)))

centroids


# merge with site scores, order levels of categorical predictors
plot.scrs <- plot.scrs %>%
  left_join(centroids, by=c("MowFreq", "Month")) %>%
  mutate(Mowing=fct_relevel(Mowing,"regular mowing", "reduced mowing", "reduced mowing & sowing")) %>% 
  mutate(Month=fct_relevel(Month,"March", "May", "July", "September")) 


plot.scrs

# plot for plots data
set.seed(11)
# plot for plots data
plot1 <- ggplot(data=plot.scrs, 
                aes(x=NMDS1, y=NMDS2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # spiders
  
  geom_segment(data = plot.scrs,        
               mapping = aes(xend = NMDS1_centroid, yend = NMDS2_centroid, 
                             color=Mowing),
               alpha=0.5) +
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
                      color=Mowing, label = Month), 
                  size=5, fontface="bold", show.legend = F) +
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
               color=NA) +
  # species
  geom_point(data=sp.scrs, 
             aes(x=NMDS1, y=NMDS2, color=trait), 
             size = 0.5,  
             alpha=0.8, pch=19)+
  geom_text_repel(data=sp.scrs, 
                  aes(x=NMDS1, y=NMDS2, label = Taxon_EuroMed,
                      color=status ), #Genus , Sociality ,  life_strategy, Nest_location
                  size=3, fontface="bold", show.legend = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
  scale_color_manual(values = c(
    "least-concerned" = "forestgreen",#"#1b9e77",   
    "endangered"    = "red4",   
    "vulnerable"     = "red",   
    "neophytes"       = "#3b4cc0",   
    "data insufficient"  = "#999999")) +
  scale_fill_manual(values = c(
    "regular mowing" = "red", ##F8766D",
    "reduced mowing" = "#00B0F8", # "#00B0F6",  #"yellow3",
    "reduced mowing & sowing" = "green" #"#00BA38" # "#00B0F6"
  )) +
  labs(color="Red-list species and neophytes", fill="Grassland management",
       x="NMDS1", y="NMDS2")

print(plot2)

