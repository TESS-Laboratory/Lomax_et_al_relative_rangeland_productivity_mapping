#' ---
#' title: "RESTREND and RUE"
#' date: 2023-11-27
#' author: "Guy Lomax"
#' ---
#' 
#' This script calculates Rain Use Efficiency (RUE) (Le Houerou 1984) and
#' Residual Trend analysis (RESTREND) (Evans & Geerken 2004; Burrell et al. 2019)
#' residuals for all grid cells in the study area, as well as temporal trends in
#' these values, and saves them to disk.
#' 
## ----setup------------------------------------------------------------------------------------------------------------------------------------------------------------

# Data handling
library(tidyverse)
library(sf)
library(terra)
library(data.table)
library(here)

# Analysis
library(future)
library(furrr)

#' 
#' 
#' Read in data layers
#' 
## ----load, include = FALSE--------------------------------------------------------------------------------------------------------------------------------------------

# Country boundaries
ke_tz <- st_read(here("data", "raw", "vector", "kenya_tanzania.geojson"))

# Covariate layers

# Load dynamic covariates as a named list of rasters (one per year)
dynamic_covariates_paths <- Sys.glob(here("data", "raw", "raster", "covariateMaps",
                                          "dynamicVars*.tif"))

years <- str_extract(dynamic_covariates_paths, "\\d\\d\\d\\d")
dynamic_covariates <- map(dynamic_covariates_paths, rast)
names(dynamic_covariates) <- years


#' 
#' 
#' Data preparation
#' 
## ----data_prep--------------------------------------------------------------------------------------------------------------------------------------------------------

# Simplify dynamic covariates raster to only GPP and precipitation
precipitation <- map(dynamic_covariates, function(r) r$precipitation) %>%
  rast()
tMean <- map(dynamic_covariates, function(r) r$tMean) %>%
  rast()
gpp <- map(dynamic_covariates, function(r) r$GPP) %>%
  rast()

names(precipitation) <- paste0("precipitation_", years)
names(tMean) <- paste0("tMean_", years)
names(gpp) <- paste0("gpp_", years)


#' 
#' 
#' 
#' Calculate RUE values and trend for all pixels
#' 
## ----rue--------------------------------------------------------------------------------------------------------------------------------------------------------------

# Calculate simple rain use efficiency

rue <- gpp / precipitation

writeRaster(rue, here("data", "processed", "raster", "rue", "rue_all.tif"),
            overwrite = TRUE)

# Mean RUE
rue_mean = mean(rue)
names(rue_mean) <- "mean_rue"

# RUE Trend using Theil Sen regression slope
rue_trend <- app(rue, function(ts) {
  if (any(is.na(ts))) {
    rep(NA, 3)
  } else {
    year <- 1:19
    df <- data.frame(year = year, ts = ts)

    mod <- RobustLinearReg::theil_sen_regression(ts ~ year, data = df)

    slope <- coef(mod)[2]
    r_sq <- broom::glance(mod)$r.squared[1]
    p_value <- broom::tidy(mod)$p.value[2]

    c(slope, r_sq, p_value)
  }
})

names(rue_trend) <- c("rue_slope", "r_sq", "p_value")

# Write rasters to disk
writeRaster(rue_mean, here("data", "processed", "raster", "rue", "rue_mean.tif"), overwrite = TRUE)
writeRaster(rue_trend, here("data", "processed", "raster", "rue", "rue_trend.tif"), overwrite = TRUE)


#' 
#' Calculate RESTREND residuals and residual trend for all pixels
#' 
## ----restrend---------------------------------------------------------------------------------------------------------------------------------------------------------

# Custom function to extract model params and residuals from linear models
# fitted to the time series in each pixel

fit_restrend <- function(gpp_rast, ppt_rast, tMean_rast) {
  
  # Convert rasters to a single data frame
  message("Converting to data frames\n")

  gpp_values <- as.data.frame(gpp_rast, cells = TRUE, xy = TRUE, na.rm = T)
  ppt_values <- as.data.frame(ppt_rast, cells = TRUE, xy = FALSE, na.rm = T)
  tMean_values <- as.data.frame(tMean_rast, cells = TRUE, xy = FALSE, na.rm = T)
  
  data_df <- gpp_values %>%
    left_join(ppt_values, by = "cell") %>%
    left_join(tMean_values, by = "cell") %>%
    pivot_longer(cols = starts_with(c("gpp", "precipitation", "tMean"))) %>%
    separate_wider_delim(cols = "name", delim = "_", names = c("var", "year")) %>%
    pivot_wider(id_cols = c(cell, x, y, year), names_from = var, values_from = value) %>%
    drop_na()
  
  # Apply RESTREND function to each cell in data frame
  message("Building models and calculating residuals\n")
  
  apply_restrend <- function(df) {
    
    model <- lm(gpp ~ precipitation + tMean, data = df)

    yint <- coef(model)[1]
    ppt_slope <- coef(model)[2]
    t_slope <- coef(model)[3]

    r_sq <- broom::glance(model)$r.squared
    ppt_p_value <- broom::tidy(model)$p.value[2]
    t_p_value <- broom::tidy(model)$p.value[3]

    resids <- residuals(model)
    
    if (length(resids) == nlyr(gpp_rast)) {
      
      output <- c(yint, ppt_slope, t_slope, r_sq, ppt_p_value, t_p_value, resids)
      
    } else {
      
      output <- rep(NA, nlyr(gpp_rast) + 6)
      
    }
    
    names(output) <- c("yint", "ppt_slope", "t_slope",
                         "rsq", "ppt_p_value", "t_p_value", 
                         paste0("resid_", years))
    
    output
  }
  
  restrend_df <- data_df %>%
    group_by(cell, x, y) %>%
    nest() %>%
    ungroup() %>%
    mutate(restrend_list = future_map(data, apply_restrend)) %>%
    select(x,y,restrend_list) %>%
    unnest_wider(restrend_list)
  
  # Convert to raster
  message("Converting back to raster")
  restrend_rast <- rast(restrend_df, crs = crs(gpp_rast), extent = ext(gpp_rast))
}

restrend_rast <- fit_restrend(gpp, precipitation, tMean)
writeRaster(restrend_rast, here("data", "processed", "raster", "restrend", "restrend_resids.tif"),
            overwrite = TRUE)

# pushoverr::pushover("RESTREND residuals calculated")

# Calculate residual trend and save

restrend_rast <- rast(here("data", "processed", "raster", "restrend", "restrend_resids.tif"))

calculate_trend <- function(x) {
  
  if (any(is.na(x))) {
    
    rep(NA, 3)
  
  } else {
    
    year <- 1:19
    mod <- RobustLinearReg::theil_sen_regression(x ~ year)
    
    resid_slope <- coef(mod)[2]
    resid_r_sq <- broom::glance(mod)$r.squared[1]
    resid_p_value <- broom::tidy(mod)$p.value[2]
    
    c(resid_slope, resid_r_sq, resid_p_value)
  }
  
}

restrend_results <- app(restrend_rast[[7:25]], calculate_trend)

names(restrend_results) <- c("resid_slope", "resid_rsq", "resid_p_value")

writeRaster(restrend_results, here("data", "processed", "raster", "restrend_results.tif"),
            overwrite = TRUE)


# pushoverr::pushover("RESTREND analysis complete")

#' 
