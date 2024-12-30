# The Relative Productivity Index: Mapping human impacts on rangeland vegetation productivity with quantile regression forests

This repository contains R and Google Earth Engine code for the analysis in Lomax, G. A., Powell, T. W. R., Lenton, T. M., and Cunliffe, A. M. (submitted), The relative productivity index: mapping human impacts on rangeland vegetation productivity with quantile regression forests monitor land degradation in East Africa using remote sensing.

Publication doi: 

Repository doi: [![DOI](https://zenodo.org/badge/724703367.svg)](https://doi.org/10.5281/zenodo.13987665)

Contact: G.Lomax@exeter.ac.uk

Google Earth Engine script repository: [https://code.earthengine.google.com/?accept_repo=users/guylomax01/Lomax_relative_rangeland_productivity_mapping](https://code.earthengine.google.com/?accept_repo=users/guylomax01/Lomax_relative_rangeland_productivity_mapping)

Google Earth Engine assets: https://code.earthengine.google.com/?asset=projects/ee-guylomax01/assets/relative_rangeland_productivity_mapping

## Introduction
This repository provides code to calculate and visualise the Relative Productivity Index (RPI), a unitless index of vegetation condition. RPI is the ratio between observed annual vegetation gross primary productivity (GPP) in a grid cell (derived from remote sensing) and the modelled potential GPP in that grid cell given its biophysical conditions. RPI is designed to control for the effect of climatic variability and spatial differences in topography, hydrology and soil properties in order to better estimate the impact of land management and other local impacts on vegetation and soil condition. RPI outperforms rain use efficiency (RUE, Le Houerou (1984)) and local net/gross primary production scaling (LNS/LGS, Prince (2004)) in explaining spatial variability in primary productivity, and performs only slightly less well than residual trend analysis (RESTREND, Evans & Geerken (2004); Burrell et al. (2019)) at controlling for temporal variability. See the manuscript linked above for more details.

## Repository structure

### Google Earth Engine scripts:
**1_generateSamplePoints** - Defines the study region and generates a random sample of points within rangeland grid cells for model training.

**2_extractSampleCovariates** - Calculates and exports the required covariates (except for topographic wetness index, calculated separately) for the sample points in GeoJSON format.

**3_extractCovariateLayers** - Calculates and exports the required covariates as continuous raster layers covering the study region. Annual/dynamic covariates are exported as separate export tasks for each year.

**helpers** - Helper function imported by scripts 2 and 3 to calculate the Unranked Gini index (UGi).

### R scripts:
Main R scripts are in the /analysis/ subdirectory and should be run in numerical order.

**1_twi.R** - R script that calculates topographic wetness index from a (hydrologically conditioned) digital elevation model.

**2_quantile_regression_rf.Rmd** - R Notebook that prepares sample data and fits quantile regression forest models, including forward feature selection and hyperparameter tuning. This script also tests random forest model performance and generates diagnostic plots (Figure 1 in manuscript).

**3_rpi_map_rf.Rmd** - R Notebook that predicts the quantile regression forest model to the full study area to derive annual predicted gross primary productivity (GPP), potential GPP and RPI values.

**4_local_gpp_scaling_df.Rmd** R Notebook that runs k-means clustering with user-defined k and multi-annual mean covariate values to implement Local GPP Scaling across the study area and calculate LGS ratio values for each grid cell.

**5_restrend.Rmd** - R Notebook that fits Theil-Sen linear regression models of GPP as a function of annual precipitation and mean temperature, then calculates annual RESTREND residuals for each grid cell.

**6_results_spatial.Rmd** - R Notebook that quantifies and visualises performance of the three spatial methods (RPI, RUE and LGS) in modelling spatial variability in GPP across three metrics - mean absolute error, R-squared and residual correlation with mean annual precipitation.

**7_results_temporal.Rmd** - R Notebook that quantifies and visualises performance of the three time series methods (RPI, RUE and RESTREND) in modelling pixel-wise temporal variability in GPP across three metrics - mean absolute error, R-squared and residual correlation with annual precipitation.

**8_rpi_pattern_maps.Rmd** - R Notebook that integrates spatial and temporal dimensions of RPI into combined visualisations (Figure 6 in the manuscript)

