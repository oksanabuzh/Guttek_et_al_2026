
# Purpose:  Canonical correspondence analysis for species composition sampled at 1m2 plots

dev.off()


library(tidyverse)
library(vegan)
library(ggrepel)

citation()

# read & wrangle data ----- 

## species data ----
species_data <- read_csv("data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, Month, remove = TRUE) %>% 
  select(Plot, Taxon, cover) %>%
  arrange(Taxon) %>%
  pivot_wider(names_from = Taxon, values_from = cover, values_fill = 0) %>% 
  column_to_rownames("Plot")

str(species_data)

## predictor data ----

# mowing data
Mowing_data <- read_csv("data/raw_data/mowing_events_2025.csv") %>% 
  pivot_longer(cols = c(September,	July,	May,	March),
               names_to = "Month",
               values_to = "n_mow_events_befre_sampling")%>% 
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
Trait_data <- read_csv("data/processed_data/Traits_Dist.Ind.Values.csv") %>% 
  mutate(status=
           case_when(is.na(status) ~ "Data insufficient",
                     status  =="NotEndangered" ~ "Not endangered",
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
# >3 -> nonlinear methods (CCA)
# in any case non metric distance based methods can be used (NMDS or PCoA)
decorana(species_data) 
#  >3 -> linear methods  are not applicable as axis lengths for DCA1 is >3

# Canonical correspondence analysis
set.seed(1)
ord_mod1 <- cca(species_data ~ MowFreq:Month + MowFreq + Month +
                  +   n_mow_events_befre_sampling, data = predictor_data,
             scale = FALSE) # scale data to have the same units
ord_mod1
anova(ord_mod1, strata = as.factor(predictor_data$PlotNo), # random effects
                                   by= "terms") # each term (sequentially from first to last), depends on the order

set.seed(1)
ord_mod2 <- cca(species_data ~ MowFreq + Month + 
                 n_mow_events_befre_sampling, data = predictor_data,
               scale = FALSE) # scale data to have the same units
ord_mod2
anova(ord_mod2, strata = as.factor(predictor_data$PlotNo), # random effects
                     by= "terms") # each term (sequentially from first to last),

ord_effects <- anova(ord_mod2, strata = as.factor(predictor_data$PlotNo), # random effects
                     by= "terms") # each term (sequentially from first to last), depends on the order
ord_effects

vif.cca(ord_mod2)
# proportion variance explained by CCA axes
summary(eigenvals(ord_mod2))
# adjusted R2
RsquareAdj(ord_mod2)

# Permutation tests ------
## --- model fit ---
set.seed(1)
Mod_sign <- anova(ord_mod2, strata = as.factor(predictor_data$PlotNo)) # model fit statistics
Mod_sign

# save results ------

write_csv(Mod_sign %>% 
            as_tibble(rownames = "Predictors") %>% 
            filter(Predictors!="Residual") %>%
            mutate(Model_R2=RsquareAdj(ord_mod2)[[2]],
                   CCA1.Prop.Explained=summary(eigenvals(ord_mod2))[[2,1]],
                   CCA2.Prop.Explained=summary(eigenvals(ord_mod2))[[2,2]]) %>% 
            bind_rows(ord_effects %>% 
                        as_tibble(rownames = "Predictors")),
          "results/CCA_species_results.csv")
plot(ord_mod)


# remove non significant effect of mowing for plotting
set.seed(1)
ord_mod <- cca(species_data ~ MowFreq + Month, 
               data = predictor_data,
               scale = FALSE) # scale data to have the same units

ord_mod

anova(ord_mod, strata = as.factor(predictor_data$PlotNo), # random effects
      by= "terms") # each term (sequentially from first to last),


# extract species scores
sp.scrs <- scores(ord_mod, display = "species",
                  scaling = "species") %>% 
  as_tibble(rownames = "Taxon") %>% 
  left_join(Trait_data, by="Taxon") %>% 
  mutate(species_full_name=Taxon,
         Taxon = if_else(
           str_count(Taxon, "\\S+") == 1,      # If only one word (non-space sequence)
           paste(Taxon, "sp."),                # add "sp."
           Taxon                               # else keep as is
         ),
         Taxon = str_c(str_split_i(Taxon, '\\s+', 1) %>%    # splits the species_name at each empty space in the species name and extracts the first word (the genus)
                           str_sub(.,  1, 5 ),          #  in this string (".") subtracts first 4 letters of genus (start, end 
                         str_split_i( Taxon, '\\s+', 2) %>%   # gets the second part of the species name after the first empty space (species)
                           str_sub(., 1, 3),            #  subtracts first 5 letters of from that second part (species)
                         sep = '.')) %>% 
  mutate(Taxon = ifelse(Taxon=="Plant.(ro", "Plantae", Taxon)) %>% 
  mutate(trait=fct_relevel(status,"Endangered",
                           "Not endangered", 
                           "Warning",
                           "Neophyte", 
                           "Data insufficient"))

sp.scrs


# extract plot scores 
plot.scrs <- scores(ord_mod, display = "sites",
                    scaling = "sites") %>% 
  as_tibble(rownames = "Plot") %>% 
  left_join(predictor_data, by="Plot") 

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
  filter(str_detect(treatment, "Month")) %>% 
  mutate(Month=stringr::str_sub(treatment, 6)) %>% 
  dplyr::select(-treatment) %>% 
  rename(CCA1_month=CCA1,
         CCA2_month=CCA2)

centroid_month

# centroid for interaction from raw data
centroids <- plot.scrs %>% 
  group_by(MowFreq, Month) %>% 
  summarise(CCA1_centroid=mean(CCA1),
            CCA2_centroid=mean(CCA2)) %>% 
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
                      color=Mowing, label = Month), 
                  size=5, fontface="bold", show.legend = F) +
  theme_bw()+
  scale_color_manual(values = c("#F8766D", "#00B0F6","#00BA38"))+
  labs(color="Management", x="CCA1 (3.8 %)", y="CCA2 (3.4 %)")

print(plot1)


ggsave("results/plots/CCA_plot1.png", plot1, width = 8, height = 6, dpi = 350)
# ggsave("CCA_plot1.jpeg", plot1, width = 8, height = 6, dpi = 350)

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
               color=NA) +
  # species
  geom_point(data=sp.scrs, 
             aes(x=CCA1, y=CCA2, color=trait), 
             size = 0.5,  
             alpha=0.8, pch=19)+
  geom_text_repel(data=sp.scrs, 
                  aes(x=CCA1, y=CCA2, label = Taxon,
                      color=status ), #Genus , Sociality ,  life_strategy, Nest_location
                  size=3, fontface="bold", show_guide = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
  scale_color_manual(values = c(
    "Not endangered" = "forestgreen",#"#1b9e77",   
    "Endangered"    = "red4",   
    "Warning"     = "red",   
    "Neophyte"       = "#3b4cc0",   
    "Data insufficient"  = "#999999")) +
  scale_fill_manual(values = c(
    "regular mowing" = "red", ##F8766D",
    "reduced mowing" = "#00B0F8", # "#00B0F6",  #"yellow3",
    "reduced mowing & sowing" = "green" #"#00BA38" # "#00B0F6"
  )) +
  labs(color="Red list species and neophytes", fill="Management",
       x="CCA1 (3.6 %)", y="CCA2 (3.1 %)")


print(plot2)
ggsave("results/plots/CCA_plot2.png", plot2, width = 10, height = 8, dpi = 350)
