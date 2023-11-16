# Spatio-Temporal-SDM-Sea-Cucumber
Code used in Hebert et al., submitted to Ecosphere

[![DOI](https://zenodo.org/badge/678377678.svg)](https://zenodo.org/badge/latestdoi/678377678)

Recommended order to run the `.R` files:
- `Survey Data Preprocessing.R`
- `Covariate VIFs and Correlations.R`
- `Environmental Raster Processing for Prediction.R`
- `mgcv Model Fitting + Output.R`
- `starve Model Fitting + Output.R`
- `sdmTMB Model Fitting + Output.R`
- `Figures 1, 2, and 3.R`

The three model fitting scripts listed above each depend on the processed survey data from `Survey Data Preprocessing.R` and the processed environmental raster data from `Environmental Raster Processing for Prediction.R`. The Figures 2 and 3 portion of `Figures 1, 2, and 3.R` depends on the three fitted models.
