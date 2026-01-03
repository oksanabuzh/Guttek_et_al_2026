# Phenology, ordination --------
dev.off

library(tidyverse)
library(vegan)
library(ggrepel)

# data
Phenophase_compos <- read_csv("data/processed_data/Diversity_phenology_1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, Month, remove = TRUE) %>% 
  select(Plot, Seedling_cover,	Juvenile_cover, Vegetative_cover, FlowerBud_cover,	Flowering_cover,	Fruiting_cover,	PostFruiting_cover) %>%
  column_to_rownames("Plot")


str(Phenophase_compos)
names(Phenophase_compos)



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


predictor_data <- Phenophase_compos %>% # make same row order as in Phenophase_compos 
  rownames_to_column("Plot") %>% select(Plot) %>%
  left_join(plot_data, by="Plot") 

str(predictor_data)


# Data exploration -----
# check both data sets for NA‘s
anyNA(Phenophase_compos) # no NA's

## Linear or nonlinear methods to use? ----
# check gradient length of first DCA axis (optional)
# if axis lengths for DCA1 is 
# <3 -> linear methods (PCA)
# >3 -> nonlinear methods (CCA)
# in any case non metric distance based methods can be used (NMDS or PCoA)
decorana((Phenophase_compos)) 
#  <3 -> linear methods  are  applicable as axis lengths for DCA1 is <3

# CCA


set.seed(1)
ord_mod <-  cca(Phenophase_compos ~ #MowFreq:Month + 
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
                   CCA1.Prop.Explained=summary(eigenvals(ord_mod))[[2,1]],
                   CCA2.Prop.Explained=summary(eigenvals(ord_mod))[[2,2]]) %>% 
            bind_rows(ord_effects %>% 
                as_tibble(rownames = "Predictors")),
          "results/CCA_phenology_results.csv")


# extract species scores
sp.scrs <- scores(ord_mod, display = "species",
                  scaling = "species") %>% 
  as_tibble(rownames = "Phenophase") %>% 
  mutate(Phenophase = str_remove(Phenophase, "_cover$")) %>% 
  mutate(Phenophase = case_when(
    Phenophase == "FlowerBud" ~ "Flower bud",
    Phenophase == "Flowering" ~ "Flowering",
    Phenophase == "Fruiting" ~ "Fruiting",
    Phenophase == "PostFruiting" ~ "Post fruiting",
    .default = Phenophase))

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
  rename( CCA1_mowing= CCA1,
          CCA2_mowing= CCA2)

centroid_mowing

centroid_month <- scores(ord_mod, 
                         display="cn",  
                         scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "Month")) %>% 
  mutate(Month=stringr::str_sub(treatment, 6)) %>% 
  dplyr::select(-treatment) %>% 
  rename( CCA1_month= CCA1,
          CCA2_month= CCA2)

centroid_month

# centroid for interaction from raw data
centroids <- plot.scrs %>% 
  group_by(MowFreq, Month) %>% 
  summarise( CCA1_centroid=mean( CCA1),
             CCA2_centroid=mean( CCA2)) %>% 
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
                aes(x= CCA1, y= CCA2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # spiders
  
  geom_segment(data = plot.scrs,        
               mapping = aes(xend =  CCA1_centroid, yend =  CCA2_centroid, 
                             color=Mowing),
               alpha=0.9) +
  # add plot scores as point:
  geom_point(data=plot.scrs, 
             aes(x= CCA1, y= CCA2, 
                 color=Mowing),
             size=1.5, pch=21) + 
  # add centroids as point:
  geom_point(data=plot.scrs, 
             aes(x= CCA1_centroid, y= CCA2_centroid, 
                 color=Mowing),
             size=3, pch=18) + 
  # centroids as text
  geom_text_repel(data=centroids, 
                  #geom_text(data=centroids, 
                  aes(x= CCA1_centroid, y= CCA2_centroid, 
                      color=Mowing, label = Month), 
                  size=5, fontface="bold", show_guide = F) +
  theme_bw()+
  scale_color_manual(values = c("#F8766D", "#00B0F6","#00BA38"))+
  labs(color="Management",  x=" CCA1 (14.8 %)", y=" CCA2 (3.4 %)")


print(plot1)


# ggsave(" CCA_plot1.png", plot1, width = 6, height = 6, dpi = 350)
# ggsave(" CCA_plot1.jpeg", plot1, width = 6, height = 6, dpi = 350)

# plot for species data
set.seed(11)
plot2 <- ggplot(data=plot.scrs, 
                aes(x= CCA1, y= CCA2))+
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
               aes(x=0, y=0, xend=CCA1, yend=CCA2), 
               arrow=arrow(length=unit(0.3,"cm")), 
               color="gray23", linewidth=1) +
  geom_text_repel(data=vector.scrs, 
                  aes(CCA1, CCA2, label="Mowing"), 
                  color="black", fontface="bold", 
                  size=5, max.overlaps = Inf) +
  # species
  geom_point(data=sp.scrs, 
             aes(x= CCA1, y= CCA2), 
             size = 0.5,  
             alpha=0.8, pch=19, color="red4")+
  geom_text_repel(data=sp.scrs, color="red4",
                  aes(x= CCA1, y= CCA2, label = Phenophase), 
                  size=4, fontface="bold", show_guide = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
  scale_fill_manual(values = c(
    "regular mowing" = "red", ##F8766D",
    "reduced mowing" = "#00B0F8", # "#00B0F6",  #"yellow3",
    "reduced mowing & sowing" = "green3" #"#00BA38" # "#00B0F6"
  )) +
   labs(color="Red list species and neophytes", fill="Management",
        x=" CCA1 (29.8 %)", y=" CCA2 (7.5 %)")

print(plot2)

# plot for species data
set.seed(11)
plot3 <- ggplot(data=plot.scrs, 
                aes(x= CCA1, y= CCA2))+
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
             aes(x= CCA1, y= CCA2), 
             size = 0.5,  
             alpha=0.8, pch=19)+
  geom_text_repel(data=sp.scrs, 
                  aes(x= CCA1, y= CCA2, label = Phenophase), 
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
  labs(fill="Month",
       x=" CCA1 (29.8 %)", y=" CCA2 (7.5 %)")


print(plot3)

