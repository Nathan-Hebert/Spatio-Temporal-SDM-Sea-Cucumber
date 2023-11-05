# Spatio-Temporal-SDM-Sea-Cucumber
Code used in Hebert et al., submitted to Ecosphere

[![DOI](https://zenodo.org/badge/678377678.svg)](https://zenodo.org/badge/latestdoi/678377678)

Recommended order to run the `.R` files:
- `Survey Data Preprocessing.R`
- `Covariate VIFs and Correlations.R`
- `Environmental Raster Processing for Prediction.R`
- `mgcv Model Fitting + Output (Includes Figure 2).R`
- `starve Model Fitting + Output (Includes Figure 5).R`
- `sdmTMB Model Fitting + Output (Includes Figure 6).R`
- `Figures 1 and 7.R`

The three model fitting scripts listed above each depend on the processed survey data from `Survey Data Preprocessing.R` and the processed environmental raster data from `Environmental Raster Processing for Prediction.R`. `Figures 1 and 7.R` depend on the three fitted models from the three model fitting scripts.
