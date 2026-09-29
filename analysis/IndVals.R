# Indicator species analysis ---------------------------------------------------

library(tidyverse)
library(indicspecies)
library(permute)

citation()

# Data ------------------------------------------------------------------------

# species cover data
Dat <- read_csv("data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv") %>%
  transmute(PlotNo,
            Subplot,
            Month,
            Management=factor(MowFreq,
                              levels=c("regular", "reduced", "reduced_sown")),
            Species=EuroMed,
            Cover=cover)

str(Dat)


# Prepare species matrix ------------------------------------------------------

# sum cover in case a species occurs more than once within the same sampling unit
species_sample <- Dat %>%
  group_by(PlotNo, Subplot, Month, Management, Species) %>%
  summarise(Cover=sum(Cover, na.rm=TRUE),
            .groups="drop")

# transform species data to wide format: samples in rows and species in columns
species_sample_wide <- species_sample %>%
  pivot_wider(id_cols=c(PlotNo, Subplot, Month, Management),
              names_from=Species,
              values_from=Cover,
              values_fill=0)

# average species cover across sampling months for each permanent 1 m2 subplot
species_subplot <- species_sample_wide %>%
  group_by(PlotNo, Subplot, Management) %>%
  summarise(across(-Month, ~mean(.x, na.rm=TRUE)),
            .groups="drop")

# number of subplots in each management regime
species_subplot %>%
  count(Management)

# management grouping vector used in multipatt()
management <- species_subplot$Management

# species matrix used in indicator species analysis
species_matrix <- species_subplot %>%
  select(-PlotNo, -Subplot, -Management) %>%
  as.data.frame()

rownames(species_matrix) <- paste(species_subplot$PlotNo,
                                  species_subplot$Subplot,
                                  sep="_")

# basic consistency checks
stopifnot(nrow(species_matrix) == length(management),
          !anyNA(management),
          all(vapply(species_matrix, is.numeric, logical(1))))

table(management)
dim(species_matrix)

# Indicator species analysis -------------------------------------------------

# IndVal.g is the group-equalised indicator value
# max.order=2 allows associations with one management regime or combinations
# of two management regimes
set.seed(123)
indval_subplot <- multipatt(species_matrix,
                            management,
                            func="IndVal.g",
                            duleg=FALSE,
                            max.order=2,
                            control=how(nperm=999))


# show species with P < 0.10 and indicator-value components A and B
summary(indval_subplot,
        alpha=0.10,
        indvalcomp=TRUE)

# Extract results -------------------------------------------------------------

# raw significance table returned by multipatt()
indval_raw <- indval_subplot$sign %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  as_tibble()

# identify columns describing membership in management groups
group_names <- levels(management)
group_columns <- intersect(group_names, names(indval_raw))

# fallback for package outputs in which group columns are not named directly
if(length(group_columns) == 0){
  possible_columns <- names(indval_raw)[
    vapply(indval_raw,
           function(x){
             is.numeric(x) && all(na.omit(unique(x)) %in% c(0, 1))
           },
           logical(1))
  ]

  group_columns <- setdiff(possible_columns, "index")[seq_along(group_names)]
}

if(length(group_columns) != length(group_names)){
  stop(paste("Could not identify management columns. Available columns:",
             paste(names(indval_raw), collapse=", ")))
}

# matrix indicating which management regime(s) are associated with each species
membership_matrix <- as.matrix(indval_raw[, group_columns, drop=FALSE])
colnames(membership_matrix) <- group_names

associated_regime <- apply(membership_matrix, 1,
                           function(x){
                             paste(group_names[x == 1], collapse=" + ")
                           })

# match species order in the significance table to specificity/fidelity matrices
species_rows <- match(indval_raw$Species, rownames(indval_subplot$A))

if(anyNA(species_rows)){
  stop("Species names do not match between $sign and the A/B matrices.")
}

# index identifies the selected single group or group combination for each species
combination_index <- indval_raw$index

# A = specificity: probability that a sampled site belongs to the target group
# when the species is found there
specificity <- indval_subplot$A[cbind(species_rows, combination_index)]

# B = fidelity: probability of finding the species within the target group
fidelity <- indval_subplot$B[cbind(species_rows, combination_index)]

# combine all indicator-value components in one table
indval_results_subplot <- tibble(
  Species=indval_raw$Species,
  `Associated regime`=associated_regime,
  `Specificity (A)`=specificity,
  `Fidelity (B)`=fidelity,
  `IndVal.g`=indval_raw$stat,
  `P-value`=indval_raw$p.value) %>%
  mutate(`Associated regime`=recode(
    `Associated regime`,
    "regular"="Regular mowing",
    "reduced"="Reduced mowing",
    "reduced_sown"="Reduced mowing and sowing",
    "regular + reduced"="Regular mowing + reduced mowing",
    "regular + reduced_sown"="Regular mowing + reduced mowing and sowing",
    "reduced + reduced_sown"="Reduced mowing + reduced mowing and sowing"))

# Supplementary table ---------------------------------------------------------

# retain all associations with P < 0.10 for the supplementary material
supplementary_subplot <- indval_results_subplot %>%
  filter(`P-value` < 0.10) %>%
  arrange(`Associated regime`, desc(`IndVal.g`)) %>%
  mutate(`Specificity (A)`=round(`Specificity (A)`, 3),
         `Fidelity (B)`=round(`Fidelity (B)`, 3),
         `IndVal.g`=round(`IndVal.g`, 3),
         `P-value`=round(`P-value`, 4))

supplementary_subplot

# species associated with individual management regimes
species_regular_subplot <- supplementary_subplot %>%
  filter(`Associated regime` == "Regular mowing")

species_reduced_subplot <- supplementary_subplot %>%
  filter(`Associated regime` == "Reduced mowing")

species_reduced_sown_subplot <- supplementary_subplot %>%
  filter(`Associated regime` == "Reduced mowing and sowing")

# species associated with both reduced-mowing treatments
species_both_reduced_subplot <- supplementary_subplot %>%
  filter(`Associated regime` ==
           "Reduced mowing + reduced mowing and sowing")


# inspect results by management regime
species_regular_subplot
species_reduced_subplot
species_reduced_sown_subplot
species_both_reduced_subplot

# Save results ----------------------------------------------------------------

dir.create("results", showWarnings=FALSE, recursive=TRUE)

write_csv(supplementary_subplot,
          "results/Table_S_IndVal_species_subplot_level.csv")


write_csv(indval_results_subplot,
          "results/IndVal_all_species_subplot_level.csv")
