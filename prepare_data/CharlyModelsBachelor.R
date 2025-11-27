
# Script for only vascular plants

# setwd("~/BIOLOGIE/Blooming Campus/Statistics") ----------------------------
library(ggplot2)
library(dplyr)
library(readr)
library(readxl)
library(treemap)
library(tidyr)
library(tidyverse)
library(vegan)
library(nlme)
library(lme4)
library(janitor)
library(ggcats)
library(labdsv)

## for citation, do i cite all packages as well? --------------------------------
# citation("base")
# RStudio.Version()
# citation("ggplot2")
# citation("dplyr")
# citation("readr")
# citation("readxl")
# citation("treemap")
# citation("tidyr")
# citation("tidyverse")
# citation("vegan")
# citation("nlme")
# citation("lme4")
# citation("janitor")
# citation("ggcats")
# citation("labdsv")


## read in data ----------------------------------------------------------------
Mayday <- read_csv("data/Sampling4.1CSVdata.csv", 
                   col_types = cols(Date = col_date(format = "%d.%m.%Y"), 
                                    `0.01` = col_number(), `0.03` = col_number(), 
                                    `0.10` = col_number(), `0.32` = col_number(), 
                                    `1.00` = col_number(), `3.16` = col_number(), 
                                    Seedling = col_number(), Juvenile = col_number(), 
                                    FlowerBud = col_number(), Flowering = col_number(), 
                                    Fruiting = col_number(), PostFruiting = col_number()))
# View(Mayday)


## problems ---------------------------------------------------------------------
Mayday$HeightSplit <- strsplit(as.character(Mayday$height), ",\\s*")
Mayday$HeightSplit <- lapply(Mayday$HeightSplit, as.numeric)
Mayday$HeightMean <- sapply(Mayday$HeightSplit, mean)
# Mayday$Date <- as.Date(Mayday$Date, format="%d.%m.%Y")
Mayday$Date <- as.character(Mayday$Date)


## read NAs in all phenotype columns as Zeros -----------------------------------
Mayday <- Mayday %>% 
  mutate(Seedling = ifelse(is.na(Seedling), 0, Seedling),
         Juvenile = ifelse(is.na(Juvenile), 0, Juvenile),
                            Flowering = ifelse(is.na(Flowering), 0, Flowering),
                            Fruiting = ifelse(is.na(Fruiting), 0, Fruiting),
                            PostFruiting = ifelse(is.na(PostFruiting), 0, PostFruiting),
                            FlowerBud = ifelse(is.na(FlowerBud), 0, FlowerBud)
         ) %>%
           mutate(Seedling = ifelse(Seedling>1, 1, Seedling),
                      Juvenile = ifelse(Juvenile>1, 1, Juvenile),
                      Flowering = ifelse(Flowering>1, 1, Flowering),
                      Fruiting = ifelse(Fruiting>1, 1, Fruiting),
                      PostFruiting = ifelse(PostFruiting>1, 1, PostFruiting),
                      FlowerBud = ifelse(FlowerBud>1, 1, FlowerBud)
                  )

Mayday <- Mayday %>% mutate(Layer = ifelse(is.na(Layer), "H", Layer))


# now with as.logical into True or False
#Mayday$Seedling <- as.logical(Mayday$Seedling)
#Mayday$Juvenile <- as.logical(Mayday$Juvenile)
#Mayday$FlowerBud <- as.logical(Mayday$FlowerBud)
#Mayday$Flowering <- as.logical(Mayday$Flowering)
#Mayday$Fruiting <- as.logical(Mayday$Fruiting)
#Mayday$PostFruiting <- as.logical(Mayday$PostFruiting)


## Sampling Date into categories -----------------------------------------------
Mayday <- Mayday %>%
  mutate(
    MonthLabel = case_when(
      substr(Date, 6, 7) %in% c("03", "04") ~ "March",
      substr(Date, 6, 7) %in% c("05", "06") ~ "May",
      substr(Date, 6, 7) %in% c("07", "08") ~ "July",
      substr(Date, 6, 7) %in% c("09", "10") ~ "September",
      TRUE ~ NA_character_
    )
  )


## add Mowing Frequency for each Plot -------------------------------------------
Mayday <- Mayday %>%
  mutate(
    MowFreq = case_when(
      substr(PlotNo, 1, 4) %in% c("BC01", "BC04", "BC06", 
                                  "BC08", "BC11", "BC13", 
                                  "BC14", "BC16", "BC19", "BC20") ~ "regular",
      substr(PlotNo, 1, 4) %in% c("BC02", "BC03", "BC07", 
                                  "BC09", "BC12", "BC15", "BC18") ~ "reduced",
      substr(PlotNo, 1, 4) %in% c("BC05", "BC10", "BC17") ~ "reduced_sown",
      TRUE ~ NA_character_
    )
  )



# Select Data for Bachelor -----------------------------------------------------
#'* only Vascular plants *
Vascular <- subset(Mayday, Layer == "H")
# with or without sown plots?

## Species biomass -------------------------------------------------------------
Vascular <- Vascular %>% 
  mutate(Biomass10m2= HeightMean * `3.16`,
         Biomass1m2 = HeightMean * `1.00`)


## Community composition big plot ----------------------------------------------
Compos10m2 <- Vascular %>% 
  group_by(PlotNo, MowFreq, MonthLabel, Taxon) %>%
  summarise(cover=mean(`3.16`, na.rm = TRUE)) %>% 
  mutate(cover=ifelse(is.na(cover), 0.01, cover)) %>% 
           ungroup() %>% 
    pivot_wider(names_from = "Taxon", values_from = "cover", values_fill = 0)
Compos10m2

write_csv(Compos10m2, "data/Compos10m2.csv")

### Phenology big plot ----------------------------------------------------------
Phenology10m2 <- Vascular %>% 
  group_by(PlotNo, MowFreq, MonthLabel, Taxon) %>%
  summarise(cover=mean(`3.16`, na.rm = TRUE),
            Seedling=sum(Seedling, na.rm = TRUE), 
            Juvenile=sum(Juvenile, na.rm = TRUE), 
            FlowerBud=sum(FlowerBud, na.rm = TRUE), 
            Flowering=sum(Flowering, na.rm = TRUE), 
            Fruiting=sum(Fruiting, na.rm = TRUE),
            PostFruiting=sum(PostFruiting, na.rm = TRUE)) %>% 
  mutate(cover=ifelse(is.na(cover), 0.01, cover)) %>% 
  ungroup()

write_csv(Phenology10m2, "data/Phenology10m2.csv")


## Composition subplot ----------------------------------------------------------
Compos1m2 <- Vascular %>%
  filter(!Subplot=="EXT") %>% 
  group_by(PlotNo, MowFreq, MonthLabel, Subplot, Taxon) %>%
  summarise(cover=mean(`1.00`, na.rm = TRUE)) %>% 
  mutate(cover=ifelse(is.na(cover), 0.01, cover)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Taxon", values_from = "cover", values_fill = 0)

Compos1m2

write_csv(Compos1m2, "data/Compos1m2.csv")

### Phenology subplot ----------------------------------------------------------

Phenology1m2 <- Vascular %>%
  filter(!Subplot=="EXT") %>% 
  group_by(PlotNo, MowFreq, MonthLabel, Subplot, Taxon) %>%
  summarise(cover=mean(`1.00`, na.rm = TRUE)) %>% 
  mutate(cover=ifelse(is.na(cover), 0.01, cover)) %>% 
  ungroup() %>% 
  left_join(Phenology10m2 %>% 
              dplyr::select(-MowFreq, -cover),
            by=c("PlotNo", "MonthLabel", "Taxon"))
  

write_csv(Phenology1m2, "data/Phenology1m2.csv")

## traits ----------------------------------------------------------------------
traits <- Vascular %>% 
  group_by(PlotNo, MowFreq, MonthLabel, Taxon) %>%
  summarise(cover=mean(`3.16`, na.rm = TRUE)) %>% 
  mutate(cover=ifelse(is.na(cover), 0.01, cover)) %>% 
  ungroup() %>% 
  group_by(Taxon) %>% 
  count() %>% 
  ungroup() %>% 
  dplyr::select(-n) %>% 
  mutate(trait = case_when(
    Taxon %in% 
      c("Arabis hirsuta", "Phleum phleoides", "Veronica polita", "Medicago minima"
      ) ~ "Threatened",
    Taxon %in% 
      c("Crepis biennis", "Rumex acetosa", "Saxifraga tridactylites", "Valerianella locusta", "Carduus nutans", "Carex caryophyllea", "Potentilla verna", "Potentilla neumanniana", "Salvia pratensis"
      ) ~ "ForeWarning",
    Taxon %in% 
      c("Arrhenatherum elatius", "Berteroa incana", "Claytonia perfoliata", "Helianthus tuberosus", "Oenothera biennis", "Oxalis dillenii", "	
Scilla bifolia", "Scilla siberica", "Sedum album", "Solidago canadensis", "Trisetum flavescens", "Veronica persica", "Viola odorata", "Cardamine hirsuta", "Conyza canadensis", "Erigeron canadensis", "Corydalis solida", "Erigeron annuus", "Juglans regia seedling", "Medicago sativa", "Robinia pseudoacacie seedling"
      ) ~ "Neophyte",
    Taxon %in% 
      c("Abies concolor", "Agrimonia", "Asteracea",  "Astragalus", "Betula tree seedling", "Rubus", "Prunus sp.", "Setaria", "Silene coronaria", "Taraxacum", "Tragopogon", "Trifolium incarnatum", "Ulmus seedling", "Valerianella", "Verbascum", "Vicia sp.", "Bromus hordeaceus", "Calamagrostis", "Carex muricata", "Carex sp.", "Castanea tree seedling", "Cerastium", "Corydalis", "Crataegus", "Crepis", "Crocus", "Digitaria", "Diplotaxis", "Elymus", "Erigeron", "Fragaria", "Galanthus", "Geum", "Helichrysum", "Hypochaeris", "Hypochaeris uniflora", "Leontodon", "Scorzoneroides", "Malva", "Medicago", "Myosotis", "Narcissus", "Oenothera", "Oxalis", "Plant sp.", "Prunus", "Prunus", "Quercus sp", "Rosaceae tree/shrub seedling", "Rosette sp.", "Arabis", "Cardamine", "Campanula", "Carduus", "Carex", "Grass sp.", "Holcus", "Prunus seedling", "Rumex", "Scabiosa", "Setaria/Digitaria", "Sorbus tree", "Vicia", "Acer seedling"
      ) ~ NA_character_,
    TRUE ~ "NotThreatened"
  ))

traits
write_csv(traits, "data/traits.csv")

# Species Richness per bigPlot and Season --------------------------------------


Diversity_10m2 <- Vascular %>%
  group_by(PlotNo, MowFreq, MonthLabel, Taxon) %>%
  summarise(cover=mean(`3.16`, na.rm = TRUE),
            HeightMean=mean(HeightMean, na.rm = TRUE)) %>% 
  left_join(traits %>% 
              mutate(neophyte=ifelse(trait=="Neophyte", 1, 0)),
            by="Taxon") %>% 
  mutate(cover=ifelse(is.na(cover), 0.01, cover),
         biomass=HeightMean*cover) %>% 
  summarise(SR  = n_distinct(Taxon),
            evenness = vegan::diversity(cover, index = "invsimpson"),
            shannon = vegan::diversity(cover, index = "shannon"),
            cover = sum(cover, na.rm = TRUE),
            biomass = sum(biomass, na.rm = TRUE),
            neophyte=sum(neophyte, na.rm = TRUE)) %>% 
  ungroup()




## Species Richness per Subplot of each Plot and Season -------------------------



Diversity_1m2 <- Vascular %>%
  filter(!Subplot=="EXT") %>% 
  group_by(PlotNo, MowFreq, MonthLabel, Subplot, Taxon) %>%
  summarise(cover=mean(`1.00`, na.rm = TRUE)) %>% 
 left_join(Vascular %>%
  filter(!Subplot=="EXT") %>% 
  group_by(PlotNo, MowFreq, MonthLabel, Taxon) %>%
    summarise(HeightMean=mean(HeightMean, na.rm = TRUE)) %>% 
    ungroup(),
  by=c("PlotNo", "MowFreq", "MonthLabel", "Taxon")
  ) %>% 
  left_join(traits %>% 
              mutate(neophyte=ifelse(trait=="Neophyte", 1, 0)),
            by="Taxon") %>% 
  mutate(cover=ifelse(is.na(cover), 0.01, cover),
         biomass=HeightMean*cover) %>% 
  summarise(SR  = n_distinct(Taxon),
            evenness = vegan::diversity(cover, index = "invsimpson"),
            shannon = vegan::diversity(cover, index = "shannon"),
            cover = sum(cover, na.rm = TRUE),
            biomass = sum(biomass, na.rm = TRUE),
            neophyte=sum(neophyte, na.rm = TRUE)) %>%
  ungroup()

Diversity_1m2

write_csv(Diversity_1m2, "data/Diversity_1m2.csv")
write_csv(Diversity_10m2, "data/Diversity_10m2.csv")


# Rote Liste Status and Neophyte category per Taxon ----------------------------
# NA means not enough data, species not specified, cultivated (like Narcissus) or not found on List
# ForVas <- data.frame(
#   Vascular %>%
#     group_by(Taxon) %>%
#     summarise(count = n()) #%>%
#   # print(n = 30)
# )
# 
# ForVas <- ForVas %>%
#   mutate(NeoOrThreatened = case_when(
#     Taxon %in% 
#       c("Arabis hirsuta", "Phleum phleoides", "Veronica polita", "Medicago minima"
#       ) ~ "Threatened",
#     Taxon %in% 
#       c("Crepis biennis", "Rumex acetosa", "Saxifraga tridactylites", "Valerianella locusta", "Carduus nutans", "Carex caryophyllea", "Potentilla verna", "Potentilla neumanniana", "Salvia pratensis"
#       ) ~ "ForeWarning",
#     Taxon %in% 
#       c("Arrhenatherum elatius", "Berteroa incana", "Claytonia perfoliata", "Helianthus tuberosus", "Oenothera biennis", "Oxalis dillenii", "	
# Scilla bifolia", "Scilla siberica", "Sedum album", "Solidago canadensis", "Trisetum flavescens", "Veronica persica", "Viola odorata", "Cardamine hirsuta", "Conyza canadensis", "Erigeron canadensis", "Corydalis solida", "Erigeron annuus", "Juglans regia seedling", "Medicago sativa", "Robinia pseudoacacie seedling"
#       ) ~ "Neophyte",
#     Taxon %in% 
#       c("Abies concolor", "Agrimonia", "Asteracea sp.",  "Astragalus sp.", "Betula tree seedling", "Rubus", "Prunus sp.", "Setaria", "Silene coronaria", "Taraxacum", "Tragopogon", "Trifolium incarnatum", "Ulmus seedling", "Valerianella", "Verbascum", "Vicia sp.", "Bromus hordeaceus", "Calamagrostis", "Carex muricata", "Carex sp.", "Castanea tree seedling", "Cerastium", "Corydalis", "Crataegus", "Crepis", "Crocus", "Digitaria", "Diplotaxis", "Elymus", "Erigeron", "Fragaria", "Galanthus", "Geum", "Helichrysum", "Hypochaeris", "Hypochaeris uniflora", "Leontodon", "Scorzoneroides", "Malva", "Medicago", "Myosotis", "Narcissus", "Oenothera", "Oxalis", "Plant sp.", "Prunus", "Prunus", "Quercus sp", "Rosaceae tree/shrub seedling", "Rosette sp.", "Arabis"
#       ) ~ NA_character_,
#     TRUE ~ "NotThreatened"
#   ))


## Merge together, keeping all of x and adding y if it fits --------------------
# NewVascular <- merge(x = Vascular, y = ForVas, by = "Taxon", all.x = TRUE)
# colnames(NewVascular)
# 
# 
# ## Reorder Columns and Rows ----------------------------------------------------
# NewVascular <- NewVascular[, c("Number", "PlotNo", "Subplot", "Location", "Date", "MonthLabel", "Taxon", "1.00", "3.16", "phen", "Seedling", "Juvenile", "FlowerBud", "Flowering", "Fruiting", "PostFruiting", "height", "HeightMean", "MowFreq", "SR10m2", "SR1m2", "Biomass10m2", "Biomass1m2", "count", "NeoOrThreatened", "corrected", "CFs", "HeightSplit", "Layer", "0.01", "0.03", "0.10", "0.32")]
# NewVascular <- NewVascular %>% arrange(Number)


# Subsets ---------------------------------------------------------------------
VasBC01 <- subset(NewVascular, PlotNo == "BC01")
VasBC02 <- subset(NewVascular, PlotNo == "BC02")
VasBC03 <- subset(NewVascular, PlotNo == "BC03")
VasBC04 <- subset(NewVascular, PlotNo == "BC04")
VasBC05 <- subset(NewVascular, PlotNo == "BC05")
VasBC06 <- subset(NewVascular, PlotNo == "BC06")
VasBC07 <- subset(NewVascular, PlotNo == "BC07")
VasBC08 <- subset(NewVascular, PlotNo == "BC08")
VasBC09 <- subset(NewVascular, PlotNo == "BC09")
VasBC10 <- subset(NewVascular, PlotNo == "BC10")
VasBC11 <- subset(NewVascular, PlotNo == "BC11")
VasBC12 <- subset(NewVascular, PlotNo == "BC12")
VasBC13 <- subset(NewVascular, PlotNo == "BC13")
VasBC14 <- subset(NewVascular, PlotNo == "BC14")
VasBC15 <- subset(NewVascular, PlotNo == "BC15")
VasBC16 <- subset(NewVascular, PlotNo == "BC16")
VasBC17 <- subset(NewVascular, PlotNo == "BC17")
VasBC18 <- subset(NewVascular, PlotNo == "BC18")
VasBC19 <- subset(NewVascular, PlotNo == "BC19")
VasBC20 <- subset(NewVascular, PlotNo == "BC20")



# max Biomass is 100%, sum of all Biomass per Plot?
# max(Vascular$Biomass10m2, na.rm=F) O_O

# checking every plot for complete covers --------------------------------------

ggplot(VasBC01, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC02, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC03, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC04, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC05, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC06, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC07, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC08, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC09, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC10, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()
# Artemisia absinthum and campestris, are they seperate?
# Verbascum and Oenothera

ggplot(VasBC11, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC12, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC13, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC14, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC15, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC16, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC17, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()
# Festuca arund. and pratensis seperate?

ggplot(VasBC18, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()

ggplot(VasBC19, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()
# 2. May sampling big plot covers missing

ggplot(VasBC20, aes(x=MonthLabel, y=`3.16`)) +
  geom_bar(stat="identity", color="black", aes(fill=Layer)) +
  facet_wrap(~Taxon) +
  scale_y_sqrt()



# messy plots ------------------------------------------------------------------
ggplot(VasBC01, aes(x=Taxon, y=HeightMean, color=MonthLabel)) +
  geom_bar(stat="identity", aes(fill=MonthLabel), position=position_dodge(width=1, preserve='single')) +
  theme(axis.text.x=element_text(angle = 90)) +
  scale_y_sqrt()

ggplot(data=NewVascular, mapping=aes(x=PlotNo, y=HeightMean)) +
  geom_boxplot(width = 0.7, aes(fill=Subplot)) +
  theme(axis.text.x=element_text(angle = 90)) +
  scale_y_sqrt() +
  scale_fill_manual(values=c("skyblue3", "firebrick", "firebrick","goldenrod1", "goldenrod1"))

ggplot(data=NewVascular, mapping=aes(x=PlotNo, y=HeightMean)) +
  geom_boxplot(width = 0.7, aes(fill=MonthLabel)) +
  theme(axis.text.x=element_text(angle = 90)) +
  scale_y_sqrt() +
  scale_fill_manual(values=c("skyblue3", "firebrick", "goldenrod1"))

ggplot(data=NewVascular, mapping=aes(x=PlotNo, y=HeightMean)) +
  geom_boxplot(width = 0.7, aes(fill=MowFreq)) +
  facet_wrap(~MonthLabel) +
  theme(axis.text.x=element_text(angle = 90)) +
  scale_y_sqrt()

ggplot(data=NewVascular, mapping=aes(x=MonthLabel, y=HeightMean)) +
  geom_violin() +
  geom_jitter(width = 0.3, aes(color=PlotNo)) +
  facet_wrap(~MowFreq) +
  theme(axis.text.x=element_text(angle = 90)) +
  scale_y_sqrt()

ggplot(NewVascular, aes(x=SR10m2, y=HeightMean, color=MonthLabel)) + 
  geom_jitter(alpha=0.5, size=2, height=0, width=0.4)+
  geom_smooth(method=lm) +
  scale_y_sqrt()

ggplot(NewVascular, aes(x=SR10m2, y=`3.16`, color=MonthLabel)) + 
  geom_jitter(alpha=0.5, size=2, height=0.1, width=0.4)+
  geom_smooth(method=lm) +
  scale_y_sqrt()

ggplot(NewVascular, aes(x=SR10m2, y=`1.00`, color=MonthLabel)) + 
  geom_jitter(alpha=0.5, size=2, height=0.1, width=0.4)+
  geom_smooth(method=lm) +
  scale_y_sqrt()

ggplot(NewVascular, aes(x=SR10m2, y=`1.00`)) + 
  geom_cat(size=1.5) +
  scale_y_sqrt()


# Models -----------------------------------------------------------------------





