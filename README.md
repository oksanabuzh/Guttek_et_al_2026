# Code abd data for Guttek et al (2027) [Applied Vegetation Science](https://onlinelibrary.wiley.com/journal/1654109x)

## Guttek C.I., Buzhdygan O., Tietjen B., Vynokurov D., Borovyk D. (2027) 
## Management effects on plant community in urban grasslands depend on biodiversity facet and growing season
### (title is subject to change) 


# Project Structure

This project is structured as follows:

```md
.gitignore
.Rproj.user/
data/
    Variables_selected.csv
    Community_composition_DungExperiment.csv
    Community_composition_VegetationPlots.csv
    Diversity_&_NMDS_data.csv
    Panoara_Dat.csv
analysis/
    01_Calculate_diversity_&_composition.R
    02_PCA_variables.R
    03_Prepare_data.R
    04_Summary_Statistics.R
    05_GLMMs.R
    06_SEM.R
    old/
results/
Panoara.Rproj
```

## R Files

### Folder `analysis`

In this folder, the raw data is prepared for analysis and the analysis is performed

- `01_Calculate_diversity_&_composition.R`: calculates evenness and community composition (NMDS) for each plot, performs the PERMANOVA analysis and ordination plots.
- `02_PCA_variables.R`: performs PCA for soil data
- `03_Prepare_data.R`: prepares the data for analysis. It reads in the raw data, renames and filters the variables, join the data sets.
- `04_Summary_Statistics.R`: summarises data and calculates summary statistics and correlations among the measures of plant community (using GLMMs)
- `05_GLMMs.R`: performs the GLMM analysis for plant species richness (field data).
- `06.1_SEM_SpRich.R`: performs the SEM analysis, calculates direct & indirect effects 
for species richness 
- `06.2_SEM_NMDS.R`: performs the SEM analysis, calculates direct & indirect effects
for community composition (NMDS)


## Data files

### Folder `data`

This folder contains the raw data files and the the processed data files.

Raw data files:
- `Variables_selected`: This is the raw data file used in the `03_Prepare_data.R` script. It contains environmental variables and data from vegetation surveys and seed experiment.
- `Community_composition_VegetationPlots.csv`: contains a list of plant species and their cover for each field plot. 
- `Community_composition_DungExperiment.csv`: contains a list of plant species and their abundances for the seed experiment. 
- `headers.csv`: This file contains environmental variables for all plots.

Processed data files:
- `Panoara_Dat.csv`: This file is created by the `03_Prepare_data.R` script. It contains clean and joined dataset of the environmental variables, biodiversity measures, the NMDS and PCA scores.


## Other Files

- `.gitignore`: This file specifies intentionally untracked files that Git should ignore.
- `Panoara.Rproj`: This is the R project file for this project.



