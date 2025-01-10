#' ---
#' title: "Local GPP Scaling"
#' output: html_notebook
#' date: 2024-01-09
#' author: "Guy Lomax"
#' editor_options: 
#'   chunk_output_type: console
#' ---
#' 
#' This notebook implements a clustering algorithm using input geospatial layers
#' to define similar "land capability classes" as defined by Prince et al. (2009)
#' and developed by Noojipady et al. (2015) and Li et al. (2020). These references
#' typically use remotely sensed NPP products or use NDVI as a proxy, and hence
#' refer to the method as "local NPP scaling" (LNS). Using GPP data, we thus
#' apply "local GPP scaling" (LGS) as an alternative.
#' 
#' 
## ----setup------------------------------------------------------------------------------------------------------------------------------------------------------------

# Data handling
library(tidyverse)
library(sf)
library(terra)
library(here)

# Analysis
library(mlr3verse)
library(tictoc)

# Visualisation
library(tmap)


#' 
#' 
#' We define land capability classes (LCCs) based on the multi-annual mean values of
#' precipitation, air temperature, photosynthetically active radiation and potential ET,
#' as well as % tree cover, soil sand fraction and slope. This is to prevent the instability
#' in LCCs that results if new classes are derived for each year in the dataset.
#' 
#' We conduct two analyses, using the following variables:
#' - Mean annual precipitation, temperature, sand fraction, tree cover, slope, PAR and PET
#' - The above variables plus mean values of precipitation intensity and timing variables
#' 
## ----load-------------------------------------------------------------------------------------------------------------------------------------------------------------

# Country boundaries
ke_tz <- st_read(here("data", "raw", "vector", "kenya_tanzania.geojson"))

# Training data used for RF
sample_points <- read_csv(here("data", "processed", "csv", "rangelandSample_processed.csv"))

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

# Reproject TWI to resolution and extent of other variables
twi_reproj <- project(twi, dynamic_covariates[[1]])


#' 
#' Prepare data for analysis:
#' 
## ----process----------------------------------------------------------------------------------------------------------------------------------------------------------

dynamic_vars <- c("GPP", "precipitation", "tMean", "parMean", "potentialET",
                  "pptIntensity", "ugi")
static_vars <- c("sand", "slope", "wriTCFrac")

# Filter to desired variables
sample_points_subset <- select(sample_points, all_of(c(dynamic_vars, static_vars)), index)

# Map over yearly vars to generate mean raster layers

dynamic_vars_mean <- map(dynamic_vars, function(name) {
  message("Layer: ", name)
  var_all_years <- map(dynamic_covariates, function(r) {
    r[[name]]
  })
  
  var_mean <- var_all_years %>%
    rast() %>%
    mean()
  
  names(var_mean) <- name
  
  var_mean
}) %>% rast()

# Combine with selected static covariates
static_vars_rast <- static_covariates[[static_vars]]

static_vars_crop <- crop(static_vars_rast, dynamic_vars_mean)

combined_rast <- c(dynamic_vars_mean, static_vars_crop, twi_reproj)

# Convert to data.frame
combined_vars_df <- as.data.frame(combined_rast, xy = TRUE, na.rm = TRUE)


#' 
#' 
#' We apply k-means clustering algorithm using only environmental covariates (i.e.,
#' not including spatial coordinates) to derive LCCs. For k-means clustering, the
#' value of k (the desired number of clusters) must be manually specified by the
#' user. We test a range of k values from 100 to 1000. A higher value of k increases
#' the granularity of the analysis, allowing smaller subsets of covariate space
#' to be mapped and reducing heterogeneity in clusters. However, higher values of
#' k also reduce the number of pixels within each cluster, reducing the precision
#' of the quantile method and making it less likely that undegraed
#' 
#' 
## ----clustering, eval = FALSE-----------------------------------------------------------------------------------------------------------------------------------------
## 
## # Set up mlr3 tasks
## 
## tsk_lgs <- sample_points_subset %>%
##   select(-GPP) %>%
##   st_drop_geometry() %>%
##   group_by(index) %>%
##   summarise(across(.cols = everything(), .fns = mean)) %>%
##   ungroup() %>%
##   select(-index) %>%
##   as_task_clust()
## 
## # Function to implement k-means with specified k
## fit_k_means <- function(k, task) {
## 
##   message("Finding clusters for k = ", k)
## 
##   # Set up learners
##   lrn_km <- lrn("clust.kmeans",
##                 predict_type = "partition",
##                 centers = k,
##                 algorithm = "Lloyd",
##                 nstart = 5,
##                 iter.max = 1000
##   )
## 
##   # Preprocessing pipelines
## 
##   po_scale <- po("scale")
## 
##   ppl_km <- as_learner(po_scale %>>% lrn_km)
## 
##   # Save clustering learner
##   message("Clustering complete")
## 
##   cluster_model <- ppl_km$train(task)
## 
## }
## 
## k_values <- c(100, 200, 1000)
## 
## tic()
## set.seed(999)
## cluster_learners <- map(k_values, fit_k_means, task = tsk_lgs)
## toc()
## 
## pushoverr::pushover("Clustering done")
## 

#' 
#' 
#' Once points are assigned to clusters, we can extract the 90th percentile of GPP
#' for each cluster as a proxy for a potential or reference GPP. The LGS scaled
#' GPP value can then be calculated as either the ratio between estimated and
#' potential GPP (possibly the difference between estimated and minimum GPP) or as
#' the simple difference between estimated and potential. Here, we use the ratio
#' in order to avoid biasing the result to areas with a larger range of GPP values
#' and for consistency with the RPI method.
#' 
## ----lgs_predict------------------------------------------------------------------------------------------------------------------------------------------------------

# Assign all pixels to cluster using nearest neighbour

all_pixels_clustered <- combined_vars_df

for (i in seq_along(k_values)) {
  col_name <- paste0("cluster_k_", k_values[i])
  
  model <- cluster_learners[[i]]
  
  message("Predicting clusters: k = ", k_values[i])
  predictions <- model$predict_newdata(combined_vars_df)$partition
  
  all_pixels_clustered[, col_name] <- predictions
}

write_rds(all_pixels_clustered, here("data", "processed", "rds", "clusters_subset.rds"))

###

all_pixels_clustered <- read_rds(here("data", "processed", "rds", "clusters_subset.rds"))

# Calculate 90th percentile GPP for each cluster
quantiles <- all_pixels_clustered %>%
  pivot_longer(starts_with("cluster_"), names_to = "k", values_to = "cluster") %>%
  mutate(k = substring(k, 11) %>% as.numeric()) %>%
  group_by(cluster, k) %>%
  mutate(potential = quantile(GPP, 0.9),
         mean = mean(GPP)) %>%
  ungroup()

# Calculate LGS ratio

lgs <- quantiles %>%
  group_by(cluster, k) %>%
  mutate(lgs_ratio = GPP / potential) %>%
  ungroup()

hist(lgs$lgs_ratio, breaks = 100)

# Convert back to raster

k_list <- c(100, 200, 1000)

for (i in k_list) {
  message("Processing for k = ", i)
  lgs_k <- filter(lgs, k == i)
  
  lgs_rast <- lgs_k %>%
    select(x, y, cluster, GPP, mean, potential, lgs_ratio) %>%
    rast(crs = crs(dynamic_covariates[[1]]), extent = ext(dynamic_covariates[[1]]))
  
  filename <- paste0("lgs_subset_", i, ".tif")
  
  writeRaster(lgs_rast,
              here("data", "processed", "raster", "lgs", filename),
              overwrite = TRUE)
}




