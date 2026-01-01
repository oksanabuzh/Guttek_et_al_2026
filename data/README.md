# Readme for folder 'data'
This folder contains datasets used for analysis

# Files in folder 'data'

`BC_2025_Cover_Data.csv`- contains cover data for all plots in 2025 
(herb cover, cryptogam cover, bare ground cover, 
old litter cover left from 2024, new litter cover from mowing 2025, litter cover (total, both old and new), 
dead wood cover,	tree leaves cover, stones & rocks cover, gravel cover, fine soil cover,
maximum herb height, max cryptogam height)

`Sampling5.0Data.csv` - contains species data for each plot and sampling campaign in 2025
(raw vegetation survey data)


# Files in folder 'processed_data'

`Commun_Spec&Phenolog_Composition_1m2.csv` - contains community species and phenological composition for each subplot and sampling campaign;
produced by [Calculate_Diversity&Phenology_1m2](../prepare_data/Calculate_Diversity&Phenology_1m2.R)

`Diversity_phenology_1m2.csv` - contains plant diversity indices for each subplot and sampling campaign;
produced by [calculate_diversity](../prepare_data/Calculate_Diversity&Phenology_1m2.R)



# Metadata for files

## [Commun_Spec&Phenolog_Composition_1m2.csv](../data/processed_data/Commun_Spec&Phenolog_Composition_1m2.R)

| Variable      | Units                       | Meaning |
|---------------|-----------------------------|---------|
| PlotNo        | character                   | Plot identifier. |
| Month         | ordered factor              | Month of observation. Ideally stored as full month name (e.g., "April") or month number; ordering should follow the calendar. |
| Taxon         | character                   | Species name (or taxon code) as used in the source data. |
| cover         | numeric (%)                 | Group-level total cover: sum(cover) across raw rows for the PlotNo × Month × Taxon group. Units are as in `vegetation2025_1m2.csv` (commonly % cover or proportion per 1 m²). |
| Seedling      | numeric (% cover)           | "Seedling" phenophase for each species weighted on it's cover. `Seedling == 0` when this phenophase was not present for this species, month and plot. |
| Juvenile      | numeric (% cover)           | "Juvenile" phenophase for each species weighted on it's cover. `Juvenile == 0` when this phenophase was not present for this species, month and plot. |
| FlowerBud     | numeric (% cover)           | "FlowerBud" phenophase for each species weighted on it's cover. `FlowerBud == 0` when this phenophase was not present for this species, month and plot. |
| Flowering     | numeric (% cover)           | "Flowering" phenophase for each species weighted on it's cover. `Flowering == 0` when this phenophase was not present for this species, month and plot. |
| Fruiting      | numeric (% cover)           | "Fruiting" phenophase for each species weighted on it's cover. `Fruiting == 0` when this phenophase was not present for this species, month and plot. |
| PostFruiting  | numeric (% cover)           | "PostFruiting" phenophase for each species weighted on it's cover. `PostFruiting == 0` when this phenophase was not present for this species, month and plot. |
