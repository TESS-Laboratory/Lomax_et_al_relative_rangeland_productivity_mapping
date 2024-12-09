# The Relative Productivity Index: Mapping human impacts on rangeland vegetation productivity with quantile regression forests

This repository contains R and Google Earth Engine code for the analysis in Lomax, G. A., Powell, T. W. R., Lenton, T. M., and Cunliffe, A. M. (submitted), The relative productivity index: mapping human impacts on rangeland vegetation productivity with quantile regression forests monitor land degradation in East Africa using remote sensing.

Publication doi: 

Repository doi: [![DOI](https://zenodo.org/badge/724703367.svg)](https://doi.org/10.5281/zenodo.13987665)

Contact: G.Lomax@exeter.ac.uk

Google Earth Engine script repository: https://code.earthengine.google.com/?accept_repo=users/guylomax01/Lomax_relative_rangeland_productivity_mapping

Google Earth Engine assets: https://code.earthengine.google.com/?asset=projects/ee-guylomax01/assets/relative_rangeland_productivity_mapping


## Repository structure

Main scripts are in the /analysis/ subdirectory and should be run in numerical order.

**1_twi.R** - R script that calculates topographic wetness index from a (hydrologically conditioned) digital elevation model
**2_quantile_regression_rf.Rmd** - R Notebook that prepares sample data and fits quantile regression forest models, including forward feature selection and hyperparameter tuning. This script also tests random forest model performance and generates diagnostic plots (Figure 1 in manuscript).
