#' ---
#' title: "Mapping relative productivity index (RPI)"
#' author: Guy Lomax
#' date: 2023-08-18
#' ---
#' 
#' This script uses quantile regression models fitted on sample data to map
#' the annual relative productivity index across Kenyan and Tanzanian rangelands
#' from 2000 to 2018.
#' 
#' 
## ----setup, include = FALSE-------------------------------------------------------------------------------------------------------------------------------------------

# Data handling
library(tidyverse)
library(sf)
library(terra)
library(here)

# Analysis
library(mlr3)
library(ranger)
library(furrr)

#' 
#' 
## ----load, include = FALSE--------------------------------------------------------------------------------------------------------------------------------------------

# Country boundaries
ke_tz <- st_read(here("data", "raw", "vector", "kenya_tanzania.geojson"))

# Covariate layers
# Static
static_covariates <- rast(here("data", "raw", "raster", "covariateMaps",
                                 "staticVars.tif"))

twi <- rast(here("data", "processed", "raster", "hydroSHEDS",
                 "hydrosheds_twi_fd8.tif"))

names(twi) <- "twi"

# Load dynamic covariates as a named list of rasters (one per year)
dynamic_covariates_paths <- Sys.glob(here("data", "raw", "raster", "covariateMaps",
                                          "dynamicVars*.tif"))

years <- str_extract(dynamic_covariates_paths, "\\d\\d\\d\\d")
dynamic_covariates <- map(dynamic_covariates_paths, rast)
names(dynamic_covariates) <- years

# Fitted models

spt_model_tuned <- read_rds(here("results", "rds", "rf_tuned_spt.rds"))
model <- spt_model_tuned$learner$model

#'
#' Data preparation
#'
## ----data_prep--------------------------------------------------------------------------------------------------------------------------------------------------------

# Reproject TWI to same resolution and CRS as other layers

twi_reproj <- project(twi, static_covariates)

# Crop covariates to common extent
# (GEE doesn't match extents properly due to resampling)

twi_crop <- crop(twi_reproj, dynamic_covariates[[1]])
static_crop <- crop(static_covariates, dynamic_covariates[[1]])

# Calculate mean ppt and mean of mean ppt day layers
mean_ppt <- map(dynamic_covariates, subset, subset = "precipitation") %>%
  rast() %>%
  mean()

mean_meanPptDay <- map(dynamic_covariates, subset, subset = "pptMeanDay") %>%
  rast() %>%
  mean()

names(mean_ppt) <- "mean_ppt"
names(mean_meanPptDay) <- "mean_pptMeanDay"

#'
#'
#' Prepare functions for data preparation, model fitting and prediction to raster
#'
## ----functions--------------------------------------------------------------------------------------------------------------------------------------------------------

# Combine raster layers and convert to data frame

raster_to_df <- function(year) {

  # Extract yearly covariates
  dynamic_rast <- dynamic_covariates[[year]]

  # Extract and join rasters
  all_covariates <- c(dynamic_rast,
                      mean_ppt, mean_meanPptDay,
                      static_crop, twi_crop)

  # Mask to Kenya/Tanzania
  all_covariates_masked <- mask(all_covariates, ke_tz)

  # Convert raster to data frame
  covariate_df <- as.data.frame(all_covariates_masked, xy = TRUE, cells = TRUE)

  # Add year information
  covariate_df$year <- as.numeric(year)

  covariate_df
}

# Prepare data frame variables for prediction

clean_df <- function(df) {
  df_clean <- df %>%
    filter(!is.na(GPP)) %>%
    remove_missing()

  df_clean
}

# Add derived variables: ppt_mean, ppt_anomaly and meanPptDay_anomaly

add_derived_variables <- function(df) {
  df_derived <- df %>%
    mutate(ppt_anomaly = (precipitation - mean_ppt) / mean_ppt * 100,
           pptMeanDay_anomaly = pptMeanDay - mean_pptMeanDay)

  df_derived
}

# Predict potential GPP values
predict_gpp <- function(df) {

  predictions_q <- predict(model, data = df, type = "quantiles", quantiles = 0.9)
  predictions_m <- predict(model, data = df)

  df_fitted <- df %>%
    mutate(quantile_pred = predictions_q$predictions,
           mean_pred = predictions_m$predictions)

  df_fitted
}

# Calculate RPI

calc_rpi <- function(df) {
  df_rpi <- mutate(df, rpi = GPP / quantile_pred)
  df_rpi
}

# Extract fitted GPP and RPI values as raster layers

df_to_raster <- function(df, target) {
  df_to_convert <- df %>%
    select(x, y, GPP, quantile_pred, mean_pred, rpi)

  rpi_raster <- rast(df_to_convert,
                     crs = crs(target),
                     extent = ext(target))

  rpi_raster
}

#'
#'
#' Loop through annual raster layers to create maps of predicted potential GPP and
#' RPI. Using a loop to reduce memory usage and ensure years are written to disk
#' as processing is completed.
#'
## ----calc_rpi-------------------------------------------------------------------------------------------------------------------------------------------

for (i in years) {

  # Progress message
  message(paste0("Processing year ", i))

  # Convert to df
  df <- raster_to_df(i)

  message("Conversion to data frame complete")

  # Prepare for prediction
  df_for_prediction <- df %>%
    clean_df() %>%
    add_derived_variables()

  message("Data frame preparation complete")
  message("Predicting potential GPP from model...")

  # Predict potential GPP values from model
  df_fitted <- predict_gpp(df_for_prediction)

  message("Model prediction complete")

  # Calculate RPI

  df_rpi <- calc_rpi(df_fitted)

  # Convert back to raster
  target_rast <- dynamic_covariates[[i]]
  rpi_rast <- df_to_raster(df_rpi, target_rast)

  message("Conversion to raster complete")

  # Write to raster file
  filename <- paste0("rf_rpi_", i, ".tif")
  writeRaster(rpi_rast,
              here("data", "processed", "raster", "rpi", filename),
              overwrite = TRUE)

  message("Raster written to file: ", filename)

  # Clear memory for next loop
  rm(df, df_for_prediction, df_fitted, df_rpi, target_rast, rpi_rast)
  gc()
}

# pushoverr::pushover("Raster prediction complete")

# Calculate temporal trend and save as raster

## ----rpi_trends-------------------------------------------------------------------------------------------------------------------------------------------------------

rpi_paths <- Sys.glob(here("data", "processed", "raster", "rpi", "rf_rpi_*.tif"))
rpi_rast <- map(rpi_paths, rast) %>%
  map(subset, subset = "rpi") %>%
  rast()

names(rpi_rast) <- paste0("rpi_", years)

# Trend (using Theil-Sen)

fit_theil_sen <- function(vector) {
  
  if (any(is.na(vector))) {
    rep(NA, 4)
  }
  
  else {
    regression_df <- tibble(years = seq_len(length(vector)),
                            rpi = vector)
    
    ts_regression <- RobustLinearReg::theil_sen_regression(rpi ~ years, regression_df)
    
    slope_params <- broom::tidy(ts_regression)[2,2:5] %>%
      as.vector() %>%
      unlist()
    
    slope_params
  }
}

rpi_ts_trend <- app(rpi_rast, fit_theil_sen)
names(rpi_ts_trend) <- c("slope", "std_error", "f_statistic", "p_value")

writeRaster(rpi_ts_trend,
            here("results", "figures", "rf_rpi_ts_trend_rast.tif"),
            overwrite = TRUE)

# pushoverr::pushover("Theil Sen slope fitted and raster saved")


