# Functional composition, ordination 
dev.off

library(tidyverse)
library(vegan)
library(ggrepel)

# data

## Functional composition data (CWM for each trait modality) ----
FuncComp <-read.csv("data/processed_data/CWM_FunctCompos_1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, Month, remove = TRUE) %>% 
  select(-starts_with("raunkiaer_"),
         -starts_with("lifeform_"),
         # -Mowing.Frequency, -Grazing.Pressure, -Soil.Disturbance,
         -Disturbance.Severity, -Disturbance.Frequency) %>%
  column_to_rownames("Plot")

names(FuncComp)
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


predictor_data <- FuncComp %>% # make same row order as in FuncComp 
  rownames_to_column("Plot") %>% select(Plot) %>%
  left_join(plot_data, by="Plot") 

str(predictor_data)


# Data exploration -----

## Linear or nonlinear methods to use? ----
# check gradient length of first DCA axis (optional)
# if axis lengths for DCA1 is 
# <3 -> linear methods (PCA)
# >3 -> nonlinear methods (CCA)
# in any case non metric distance based methods can be used (NMDS or PCoA)
decorana((FuncComp)) 
#  <3 -> linear methods  are  applicable as axis lengths for DCA1 is <3

# we would use CCA to be consistent with all other analysis


set.seed(1)
ord_mod <-  rda(log1p(FuncComp) ~ # MowFreq:Month + 
                  MowFreq + Month + 
                  n_mow_events_befre_sampling, data = predictor_data,
                scale = FALSE) # scale data to have the same units
ord_mod
ord_effects <- anova(ord_mod, strata = as.factor(predictor_data$PlotNo), # random effects
                     by= "terms") # each term (sequentially from first to last), depends on the order

ord_effects


vif.cca(ord_mod)
# proportion variance explained by CCA axes
summary(eigenvals(ord_mod))
# adjusted R2
RsquareAdj(ord_mod)

# Permutation tests ------
## --- model fit ---
set.seed(1)
Mod_sign <- anova(ord_mod, strata = as.factor(predictor_data$PlotNo)) # model fit statistics
Mod_sign

# save results ------

write_csv(Mod_sign %>% 
            as_tibble(rownames = "Predictors") %>% 
            filter(Predictors!="Residual") %>%
            mutate(Model_R2=RsquareAdj(ord_mod)[[2]],
                   RDA1.Prop.Explained=summary(eigenvals(ord_mod))[[2,1]],
                   RDA2.Prop.Explained=summary(eigenvals(ord_mod))[[2,2]]) %>% 
            bind_rows(ord_effects %>% 
                        as_tibble(rownames = "Predictors")),
          "results/RDA_FunctComp_results.csv")


# extract species scores
sp.scrs <- scores(ord_mod, display = "species",
                  scaling = "species") %>% 
  as_tibble(rownames = "Traits") %>% 
  separate(Traits,
           into = c("Trait_group", "Trait_modality"),
           remove = FALSE,
           sep = "_",
           extra = "merge",  # keep the rest after the first "_" together
           fill = "left") %>%    # if no "_" -> Trait_modality = is Traits
  mutate(Trait_group=case_when(
    Trait_group == "FG" ~ "Functional group",
    Trait_group == "lifespan" ~ "Lifespan",
    Trait_group == "lifeform" ~ "Life form",
    Trait_group == "raunkiaer" ~ "Raunkiaer life form",
    is.na(Trait_group) ~ "Disturbance indicator",
    Trait_group =="status" ~ "Red-list / aliens")) %>% 
  mutate(Trait_modality=case_when(
    Trait_modality == "NotEndangered" ~ "least-concern",
    Trait_modality == "Endangered" ~ "endangered",
    Trait_modality == "Warning" ~ "vulnerable",
    Trait_modality == "Neophyte" ~ "neophyte",
    Trait_modality == "tree.shrub" ~ "tree/shrub",
    Trait_modality == "herbPoli" ~ "herb.polyc",
    Trait_modality == "herbMono" ~ "herb.monoc",
    Trait_modality == "shortlived" ~ "short-lived",
    Trait_modality == "other_forb" ~ "other forb",
    Trait_modality == "Disturbance.Severity" ~ "Distr.Severity",
    Trait_modality == "Disturbance.Frequency" ~ "Distr.Frequency",
    Trait_modality == "Mowing.Frequency" ~ "Mowing.Distr",
    Trait_modality == "Grazing.Pressure" ~ "Grazing.Distr",
    Trait_modality == "Soil.Disturbance" ~ "Soil.Distr",
    .default=Trait_modality)) 

sp.scrs


# extract plot scores 
plot.scrs <- scores(ord_mod, display = "sites",
                    scaling = "sites") %>% 
  as_tibble(rownames = "Plot") %>% 
  left_join(predictor_data, by="Plot") 

plot.scrs

# vector 
vector.scrs <- scores(ord_mod, display = "bp", # vector
                      scaling = "species") %>% 
  as_tibble(rownames = "Plot") %>% 
  filter(Plot=="n_mow_events_befre_sampling")  

vector.scrs


# calculate centroid for  Grazing_season
centroid_mowing <- scores(ord_mod, 
                          display="cn",  
                          scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "MowFreq")) %>% 
  mutate(MowFreq=stringr::str_sub(treatment, 8)) %>% 
  dplyr::select(-treatment) %>% 
  rename( RDA1_mowing= RDA1,
          RDA2_mowing= RDA2)

centroid_mowing

centroid_month <- scores(ord_mod, 
                         display="cn",  
                         scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "Month")) %>% 
  mutate(Month=stringr::str_sub(treatment, 6)) %>% 
  dplyr::select(-treatment) %>% 
  rename( RDA1_month= RDA1,
          RDA2_month= RDA2)

centroid_month

# centroid for interaction from raw data
centroids <- plot.scrs %>% 
  group_by(MowFreq, Month) %>% 
  summarise( RDA1_centroid=mean( RDA1),
             RDA2_centroid=mean( RDA2)) %>% 
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
                aes(x= RDA1, y= RDA2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # spiders
  
  geom_segment(data = plot.scrs,        
               mapping = aes(xend =  RDA1_centroid, yend =  RDA2_centroid, 
                             color=Mowing),
               alpha=0.9) +
  # add plot scores as point:
  geom_point(data=plot.scrs, 
             aes(x= RDA1, y= RDA2, 
                 color=Mowing),
             size=1.5, pch=21) + 
  # add centroids as point:
  geom_point(data=plot.scrs, 
             aes(x= RDA1_centroid, y= RDA2_centroid, 
                 color=Mowing),
             size=3, pch=18) + 
  # centroids as text
  geom_text_repel(data=centroids, 
                  #geom_text(data=centroids, 
                  aes(x= RDA1_centroid, y= RDA2_centroid, 
                      color=Mowing, label = Month), 
                  size=5, fontface="bold", show_guide = F) +
  theme_bw()+
  scale_color_manual(values = c("#F8766D", "#00B0F6","#00BA38"))+
  labs(color="Management",  x=" RDA1 (9.3 %)", y=" RDA2 (5.8 %)")

print(plot1)


# ggsave(" RDA_plot1.png", plot1, width = 6, height = 6, dpi = 350)
# ggsave(" RDA_plot1.jpeg", plot1, width = 6, height = 6, dpi = 350)

# plot for species data
set.seed(11)
plot2 <- ggplot(data=plot.scrs, 
                aes(x= RDA1, y= RDA2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # ellipse 
  stat_ellipse(aes(fill=Mowing), alpha=0.2,
               type='t', # type = 't' means the ellipses are calculated assuming a multivariate t-distribution (robust to outliers)
               linewidth =0.0001, geom="polygon",
               level=0.95, # ellipses represent a 95% confidence interval for the multivariate mean of each group) +
               color="gray88") +
  # vector
  geom_segment(data=vector.scrs, 
               aes(x=0, y=0, xend=RDA1, yend=RDA2), 
               arrow=arrow(length=unit(0.3,"cm")), 
               color="gray23", linewidth=1) +
  geom_text_repel(data=vector.scrs, 
                  aes(RDA1, RDA2, label="Mowing"), 
                  color="black", fontface="bold", 
                  size=5, max.overlaps = Inf) +
  # species
  geom_point(data=sp.scrs, 
             aes(x= RDA1, y= RDA2, color=Trait_group), 
             size = 0.5,  
             alpha=0.8, pch=19)+
  geom_text_repel(data=sp.scrs, 
                  aes(x= RDA1, y= RDA2, color=Trait_group,
                      label = Trait_modality), 
                  size=4, fontface="bold", show_guide = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
  scale_fill_manual(values = c(
    "regular mowing" = "red", ##F8766D",
    "reduced mowing" = "#00B0F8", # "#00B0F6",  #"yellow3",
    "reduced mowing & sowing" = "green3" #"#00BA38" # "#00B0F6"
  )) +
  labs(color="Functional category", fill="Management",
       x=" RDA1 (9.3 %)", y=" RDA2 (5.8 %)")


print(plot2)

# plot for species data
set.seed(11)
plot3 <- ggplot(data=plot.scrs, 
                aes(x= RDA1, y= RDA2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # ellipse 
  stat_ellipse(aes(fill=Month), alpha=0.3,
               type='t', # type = 't' means the ellipses are calculated assuming a multivariate t-distribution (robust to outliers)
               linewidth =0.0001, geom="polygon",
               level=0.95, # ellipses represent a 95% confidence interval for the multivariate mean of each group) +
               color="gray88") +
  # species
  geom_point(data=sp.scrs, 
             aes(x= RDA1, y= RDA2, color=Trait_group,), 
             size = 0.5,  
             alpha=0.8, pch=19)+
  geom_text_repel(data=sp.scrs, 
                  aes(x= RDA1, y= RDA2, label = Trait_modality, 
                      color=Trait_group), 
                  size=4, fontface="bold", show_guide = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
  scale_fill_manual(values = c(
    "March" = "orange",
    "May" = "#287271",
    "July" = "#6D326D",
    "September"="brown"
  )) +
  labs(fill="Month", color="Functional category",
       x=" RDA1 (9.3 %)", y=" RDA2 (5.8 %)")


print(plot3)

