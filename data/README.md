# Readme for folder 'data'
This folder contains datasets used for analysis
## Data files and metadata 


# Files in folder 'processed_data'

### [vegetation2025_1m2.csv](../data/processed_data/vegetation2025_1m2.csv)
Contains species cover, height, biomass and phenological stage for each species, subplot and sampling campaign;
produced by [data_wrangling](../prepare_data/data_wrangling.R)

| Variable      | Units / 1m²-plot each sampling month | Meaning |
|---------------|--------------------------------------|---------|
| PlotNo        | character                            | Plot identifier. |
| Subplot       | character                            | Subplot identifier: NW (N), SE (S). |
| Month         | ordered factor                       | Month of observation. |
| Taxon         | character                            | Species name  |
| EuroMed       | character                            | Species name based on the Euro+Med Taxon |
| cover         | numeric (%)                          | Cover % for each species. |
| height        | numeric (cm)                         | Mean height (cm) of each species (average of n individuals measured on each 10m2-plot measured in field for living aboveground plant). |
| Biomass       | numeric (index)                      | Above-ground biomass (index by [Fischer et al.](https://doi.org/10.1111/jvs.13182)) estimated for each species from the cover and height of living plant parts. | 
| Seedling      | binary (0/1)                         | "Seedling" phenophase for each species, coded as 0-absent, 1-present.|
| Juvenile      | binary (0/1)                         | "Juvenile" phenophase for each species, coded as 0-absent, 1-present.| 
| FlowerBud     | binary (0/1)                         | "FlowerBud" phenophase for each species, coded as 0-absent, 1-present.|
| Flowering     | binary (0/1)                         | "Flowering" phenophase for each species, coded as 0-absent, 1-present.|
| Fruiting      | binary (0/1)                         | "Fruiting" phenophase for each species, coded as 0-absent, 1-present.|
| PostFruiting  | binary (0/1)                         | "PostFruiting" phenophase for each species, coded as 0-absent, 1-present.|


### [Commun_Spec&Phenolog_Composition_1m2.csv](../data/processed_data/Commun_Spec&Phenolog_Composition_1m2.csv)
Contains community species and phenological composition for each subplot and sampling campaign;
produced by [Calculate_Diversity&Phenology_1m2](../prepare_data/Calculate_Diversity&Phenology_1m2.R)

| Variable      | Units / 1m²-plot each sampling month | Meaning |
|---------------|--------------------------------------|---------|
| PlotNo        | character                            | Plot identifier. |
| Subplot       | character                            | Subplot identifier: NW (N), SE (S). |
| Month         | ordered factor                       | Month of observation. |
| Taxon         | character                            | Species name  |
| EuroMed       | character                            | Species name based on the Euro+Med Taxon |
| cover         | numeric (%)                          | Cover % for each species. |
| height        | numeric (cm)                         | Mean height (cm) of each species (average of n individuals measured on each 10m2-plot measured in field for living aboveground plant). |
| Biomass       | numeric (index)                      | Above-ground biomass (index by [Fischer et al.](https://doi.org/10.1111/jvs.13182)) estimated for each species from the cover and height of living plant parts. | 
| Seedling      | numeric (% cover)                    | "Seedling" phenophase for each species weighted on it's cover. `Seedling == 0` when this phenophase was not present for this species, month and plot. |
| Juvenile      | numeric (% cover)                    | "Juvenile" phenophase for each species weighted on it's cover. `Juvenile == 0` when this phenophase was not present for this species, month and plot. |
| FlowerBud     | numeric (% cover)                    | "FlowerBud" phenophase for each species weighted on it's cover. `FlowerBud == 0` when this phenophase was not present for this species, month and plot. |
| Flowering     | numeric (% cover)                    | "Flowering" phenophase for each species weighted on it's cover. `Flowering == 0` when this phenophase was not present for this species, month and plot. |
| Fruiting      | numeric (% cover)                    | "Fruiting" phenophase for each species weighted on it's cover. `Fruiting == 0` when this phenophase was not present for this species, month and plot. |
| PostFruiting  | numeric (% cover)                    | "PostFruiting" phenophase for each species weighted on it's cover. `PostFruiting == 0` when this phenophase was not present for this species, month and plot. |


### [Diversity_phenology_1m2.csv](../data/processed_data/Diversity_phenology_1m2.csv)
contains plant diversity indices and community phenology weighted on species cover for each subplot and sampling campaign;
produced by [Calculate_Diversity&Phenology_1m2](../prepare_data/Calculate_Diversity&Phenology_1m2.R)

| Variable      | Units / 1m²-plot each sampling month | Meaning |
|---------------|--------------------------------------|---------|
| PlotNo        | character                            | Plot identifier. |
| Subplot       | character                            | Subplot identifier: NW (N), SE (S). |
| Month         | ordered factor                       | Month of observation. |
| cover         | numeric (%)                          | Cover % of entire community. |
| biomass       | numeric (index)                      | Community above-ground biomass (index by [Fischer et al.](https://doi.org/10.1111/jvs.13182)) sum of biomass of all species (estimated from the cover and height of living plant parts for each species). | 
| height_mean   | numeric (cm)                         | Mean height (cm) of all species in plant community. |
| height_max    | numeric (cm)                         | Maximum height (cm) of all species in plant community. |
| Seedling      | numeric (% cover)                    | Sum of cover of all species with "Seedling" phenophase.|
| Juvenile      | numeric (% cover)                    | Sum of cover of all species with "Juvenile" phenophase.|
| FlowerBud     | numeric (% cover)                    | Sum of cover of all species with "FlowerBud" phenophase for each species weighted on it's cover. |
| Flowering     | numeric (% cover)                    | Sum of cover of all species with "Flowering" phenophase for each species weighted on it's cover. |
| Fruiting      | numeric (% cover)                    | Sum of cover of all species with "Fruiting" phenophase for each species weighted on it's cover. |
| PostFruiting  | numeric (% cover)                    | Sum of cover of all species with "PostFruiting" phenophase for each species weighted on it's cover. |





# Files in folder 'raw_data'

`BC_2025_Cover_Data.csv`- contains cover data for all plots in 2025 
(herb cover, cryptogam cover, bare ground cover, 
old litter cover left from 2024, new litter cover from mowing 2025, litter cover (total, both old and new), 
dead wood cover,	tree leaves cover, stones & rocks cover, gravel cover, fine soil cover,
maximum herb height, max cryptogam height)

`Sampling5.0Data.csv` - contains species data for each plot and sampling campaign in 2025
(raw vegetation survey data)