# Functional composition, ordination 
dev.off()

library(tidyverse)
library(vegan)
library(ggrepel)

# data

## Functional composition data (CWM for each trait modality) ----
FuncComp_all <-read.csv("data/processed_data/CWM_FunctCompos_1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, Month, remove = TRUE) 

names(FuncComp_all)

ggcorrplot::ggcorrplot(round(cor(FuncComp_all %>% select(-Plot),method = c("pearson"), use = "pairwise.complete.obs"), 2),
                       hc.order = F, type = "lower",
                       lab = TRUE, lab_size = 3,
                       colors = c("red", "white", "blue"))


FuncComp <-read.csv("data/processed_data/CWM_FunctCompos_1m2.csv") %>% 
  unite(Plot, PlotNo, Subplot, Month, remove = TRUE) %>% 
  select(-starts_with("raunkiaer_"),
         -starts_with("lifeform_"),
         # -status_NotEndangered,
         # -Mowing.Frequency, 
         -Disturbance.Frequency,-Disturbance.Severity, -Grazing.Pressure #, -Soil.Disturbance
        ) %>%
  column_to_rownames("Plot")

names(FuncComp)

ggcorrplot::ggcorrplot(round(cor(FuncComp, method = c("pearson"), use = "pairwise.complete.obs"), 2),
                       hc.order = F, type = "lower",
                       lab = TRUE, lab_size = 3,
                       colors = c("red", "white", "blue"))

## predictor data ----
# mowing data
Mowing_data <- read_csv("data/raw_data/mowing_events_2025_DB.csv") %>% 
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
# >3 -> nonlinear methods (NMDS)
# in any case non metric distance based methods can be used (NMDS or PCoA)
decorana((FuncComp)) 
#  <3 -> linear methods  are  applicable as axis lengths for DCA1 is <3

# we use NMDS

# NMDS -----
# NMDS -----------------------
set.seed(2435)
ord_mod <- metaMDS(wisconsin(FuncComp), 
                   scale = FALSE, distance = "bray") 

ord_mod

# NMDS fit
ord_mod$stress
# fit
stressplot(ord_mod, main = "Shepard plot")


# Permutation test:
set.seed(10)
PERM1 <- vegan::adonis2(wisconsin(FuncComp) ~ 
                          MowFreq + Month + 
                          n_mow_events_befre_sampling, 
                        data=plot_data,
                        permutations = 1000, method = "bray",
                        strata=as.factor(plot_data$PlotNo),
                        by = "terms")

PERM1



# variable fitting for posthoc plotting  --------------
set.seed(1259)
fit1 <- vegan::envfit(ord_mod   ~  
                        MowFreq + Month + 
                        n_mow_events_befre_sampling, 
                      data=plot_data,
                      #  strata=as.factor(plot_data$PlotNo),
                      perm=1000) #


fit1



# exploratory plot
plot(ord_mod, main = "NMDS plot", display = "sites")
plot(ord_mod, main = "NMDS plot", display = "species")
plot(ord_mod, main = "NMDS plot")
plot(fit1)

### Plotting NMDS results using the ggplot --------------------------------------


# extract species scores   ---------------------
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
    is.na(Trait_group) ~ "Disturbance indicator species",
    Trait_group =="status" ~ "Red-list / aliens")) %>% 
  mutate(Trait_modality=case_when(
    Trait_modality == "NotEndangered" ~ "least-concerned",
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
    Trait_modality == "Mowing.Frequency" ~ "mowing-disturbance species",
    Trait_modality == "Grazing.Pressure" ~ "Grazing.Distr",
    Trait_modality == "Soil.Disturbance" ~ "soil-disturbance species",
    .default=Trait_modality)) 

sp.scrs


# extract plot scores    ---------------------
plot.scrs <- scores(ord_mod, display = "sites",
                    scaling = "sites") %>% 
  as_tibble(rownames = "Plot") %>% 
  left_join(predictor_data, by="Plot") 

plot.scrs

# vector     ---------------------
vector.scrs <- scores(fit1, display = "bp", # vector
                    scaling = "species") %>% 
  as_tibble(rownames = "Plot") %>% 
  filter(Plot=="n_mow_events_befre_sampling") %>% 
  mutate(NMDS1=NMDS1*9.5, 
         NMDS2=NMDS2*9.5) # rescale for better visualization

vector.scrs


# calculate centroid for  Grazing_season
centroid_mowing <- scores(fit1, 
                          display="cn",  
                          scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "MowFreq")) %>% 
  mutate(MowFreq=stringr::str_sub(treatment, 8)) %>% 
  dplyr::select(-treatment) %>% 
  rename( NMDS1_mowing= NMDS1,
          NMDS2_mowing= NMDS2)

centroid_mowing

centroid_month <- scores(fit1, 
                         display="cn",  
                         scaling="species") %>%   
  as_tibble(rownames = "treatment")  %>%
  filter(str_detect(treatment, "Month")) %>% 
  mutate(Month=stringr::str_sub(treatment, 6)) %>% 
  dplyr::select(-treatment) %>% 
  rename( NMDS1_month= NMDS1,
          NMDS2_month= NMDS2)

centroid_month

# centroid for interaction from raw data
centroids <- plot.scrs %>% 
  group_by(MowFreq, Month) %>% 
  summarise( NMDS1_centroid=mean( NMDS1),
             NMDS2_centroid=mean( NMDS2)) %>% 
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

# plots -----

## plot for plots data ---------------------------------------------------------
set.seed(11)
plot1 <- ggplot(data=plot.scrs, 
                aes(x= NMDS1, y= NMDS2))+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  # spiders
  
  geom_segment(data = plot.scrs,        
               mapping = aes(xend =  NMDS1_centroid, yend =  NMDS2_centroid, 
                             color=Mowing),
               alpha=0.9) +
  # add plot scores as point:
  geom_point(data=plot.scrs, 
             aes(x= NMDS1, y= NMDS2, 
                 color=Mowing),
             size=1.5, pch=21) + 
  # add centroids as point:
  geom_point(data=plot.scrs, 
             aes(x= NMDS1_centroid, y= NMDS2_centroid, 
                 color=Mowing),
             size=3, pch=18) + 
  # centroids as text
  geom_text_repel(data=centroids, 
                  #geom_text(data=centroids, 
                  aes(x= NMDS1_centroid, y= NMDS2_centroid, 
                      color=Mowing, label = Month), 
                  size=5, fontface="bold", show.legend = F) +
  theme_bw()+
  scale_color_manual(values = c("#F8766D", "#00B0F6","#00BA38"))+
  labs(color="Management")

print(plot1)

# ggsave(" NMDS_plot1.png", plot1, width = 6, height = 6, dpi = 350)
# ggsave(" NMDS_plot1.jpeg", plot1, width = 6, height = 6, dpi = 350)

## plot for species data ----
set.seed(11)
plot2 <- ggplot(data=plot.scrs, 
                aes(x= NMDS1, y= NMDS2))+
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
               aes(x=0, y=0, xend=NMDS1, yend=NMDS2), 
               arrow=arrow(length=unit(0.3,"cm")), 
               color="gray23", linewidth=1) +

  geom_text(data=vector.scrs, 
                  aes(NMDS1, NMDS2, label="Mowing"), 
                  color="black", fontface="bold", 
                  size=5, hjust=0.5, vjust=1) +
  # species
  geom_point(data=sp.scrs, 
             aes(x= NMDS1, y= NMDS2, color=Trait_group), 
             size = 3,  
             alpha=1, pch=8)+
  geom_text_repel(data=sp.scrs, 
                  aes(x= NMDS1, y= NMDS2, color=Trait_group,
                      label = Trait_modality), 
                  size=4, fontface="bold", show.legend = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
  scale_fill_manual(values = c(
    "regular mowing" = "red", ##F8766D",
    "reduced mowing" = "#00B0F8", # "#00B0F6",  #"yellow3",
    "reduced mowing & sowing" = "green3" #"#00BA38" # "#00B0F6"
  )) +
   labs(color="Functional category", fill="Management")


print(plot2)

## month plot for traits data -----
set.seed(11)
plot3 <- ggplot(data=plot.scrs, 
                aes(x= NMDS1, y= NMDS2))+
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
             aes(x= NMDS1, y= NMDS2, color=Trait_group,), 
             size = 3,  
             alpha=1, pch=8)+
  geom_text_repel(data=sp.scrs, 
                  aes(x= NMDS1, y= NMDS2, label = Trait_modality, 
                      color=Trait_group), 
                  size=4, fontface="bold", show.legend = F,
                  max.overlaps=Inf) +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3))) + # makes legend dots large
    scale_fill_manual(values = c(
      "March" = "orange",
      "May" = "#287271",
      "July" = "#6D326D",
      "September"="brown"
    )) +
  labs(color="Functional category", fill="Month", 
       x=" NMDS1 (14.5 %)", y=" NMDS2 (5.6 %)")


print(plot3)


