#' ---
#' title: "Quantile Regression of Rangeland Productivity - Quantile Regression Forests"
#' author: Guy Lomax
#' date: 2023-06-19
#' ---
#' 
#' This script conducts quantile regression on rangeland productivity
#' data using Quantile Regression Forests (QRF), modelling annual gross
#' primary productivity (GPP) as a function of environmental variables
#' 
## ----setup, include = FALSE-------------------------------------------------------------------------------------------------------------------------------------------

# Data handling
library(tidyverse)
library(sf)
library(here)

# Analysis
library(mlr3verse)
library(mlr3spatiotempcv)
library(ranger)
library(future)

# # Parallelisation
nc <- 16
# plan(multicore, workers = nc)


#' 
## ----load--------------------------------------------------------------------------------------------------------------------------------------------

# Country boundaries
ke_tz <- st_read(here("data", "raw", "vector", "kenya_tanzania.geojson"))

# Sample point data
sample_points_raw <- st_read(here("data", "raw", "vector", "rangelandSample.geojson"))

twi <- read_csv(here("data", "processed", "csv", "twi_sample.csv"))


#' 
#' Prepare data and add derived variables:
#' 
## ----clean-------------------------------------------------------------------------------------------------------------------------------------------

message("Preparing data\n")

yearly_vars <- c("GPP", "precipitation", "pptIntensity", "pptMeanDay", "ugi",
                 "tMean", "parMean", "potentialET")

sample_points <- sample_points_raw %>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) %>%
  bind_cols(twi) %>%
  rename(twi = hydrosheds_twi_fd8) %>%
  select(-allRangelands, -igbp, -shrubland, -grassland, -trees, -ID, -id) %>%
  pivot_longer(cols = starts_with(yearly_vars),
               names_to = "var", values_to = "value") %>%
  separate_wider_delim(var, "_", names = c("data", "year")) %>%
  pivot_wider(names_from = data, values_from = value) %>%
  mutate(year = as.numeric(year)) %>%
  filter(GPP > 0)

# Calculate derived and normalised variables
# Convert landform to factor variable
sample_points_derived <- sample_points %>%
  group_by(index) %>%
  mutate(mean_ppt = mean(precipitation),
         ppt_anomaly = (precipitation - mean_ppt) / mean_ppt * 100,
         mean_pptMeanDay = mean(pptMeanDay),
         pptMeanDay_anomaly = pptMeanDay - mean_pptMeanDay,
         landform = factor(landform)) %>%
  ungroup()

# Save sample points data frame with all variables
write_csv(sample_points_derived, here("data", "processed", "csv", "rangelandSample_processed.csv"))

#' 
#' # Data preparation for RF modelling
#' 
#' We select variables to use for modelling, remove rows containing NA
#' values (the TWI raster has some pixels with NA values) and convert to a
#' spatiotemporal regression task. We assign point id ("id") as our index
#' of spatial location and "year" as our index of time.
#' 
#' We then choose our learner (random forest regression through the ranger
#' package) with quantile regression enabled, and define two alternative
#' resampling (cross- validation) strategies: random CV and
#' "Leave-location-and-time-out" CV that ensures the validation set for
#' each fold does not overlap in either time or space with the training
#' set.
#' 
#' When defining the learner, we initially choose crude, plausible
#' hyperparameters of mtry.ratio = 0.5 and sample.fraction = 0.5. However,
#' we set the initial minimum node size to 100 (around 0.05% of the dataset)
#' to mitigate overfitting and to allow sufficient terminal node size that
#' quantiles can be calculated. For our purposes, absolute predictive
#' accuracy is less important than retaining enough variability in nodes to
#' allow assessment of quantiles.
#' 
## ----data_prep--------------------------------------------------------------------------------------------------------------------------------------------------------

# Create task

vars_to_retain <- c(
  "index",
  "x", "y",
  "year",
  "GPP",
  # "precipitation",  # Exclude raw precipitation data
  "mean_ppt",
  "ppt_anomaly",
  "pptMeanDay_anomaly",
  "pptIntensity",
  "ugi",
  "slope",
  "distToRiver",
  "landform",
  "twi",
  "tMean",
  "parMean",
  "potentialET",
  "sand",
  "wriTCFrac"
)

sample_points_subset <- sample_points_derived %>%
  drop_na() %>%
  select(all_of(vars_to_retain))

task_gpp <- as_task_regr_st(
  x = sample_points_subset,
  id = "potential_gpp",
  target = "GPP",
  coordinate_names = c("x", "y"),
  coords_as_features = FALSE,
  crs = "epsg:4326"
)

task_gpp$set_col_roles("index", roles = "space")
task_gpp$set_col_roles("year", roles = "time")

# Create learner (ranger random forest) with initial tuning values

lrn_ranger_untuned <- lrn("regr.ranger", 
                          predict_type = "response",
                          num.trees = 1001,
                          mtry.ratio = 0.5,
                          min.node.size = 100,
                          sample.fraction = 0.5
)

# Create resampling strategy

spt_cv_plan <- rsmp("sptcv_cstf", folds = 10)
# nspt_cv_plan <- rsmp("cv", folds = 10)  # Non-spatiotemporal CV (not used in main analysis)


#' 
#' # Feature selection
#' 
#' We use forward feature selection with the untuned model to identify
#' which covariates have unique information to predict GPP. For this, we
#' start by building all possible single-variable models. Then, selecting
#' the variable that gives the best performance (lowest RMSE in this case),
#' we build all possible two-variable models that include this variable. We
#' do the same with a third variable, and so on, until the point at which
#' adding an additional variable fails to improve model performance.
#' 
#' We do this separately for models using random CV and spatio-temporal CV,
#' and save the optimum feature set for each.
#' 
## ----rf_feature_select, message = FALSE, results = FALSE, include = FALSE---------------------------------------------------------------------------------------------

message("Beginning feature selection\n")

# Create forward feature selection with untuned model, minimising rmse

perf_msr <- msr("regr.mae")

fs_method <- fs("sequential", min_features = 1)

fs_term <- trm("stagnation_batch")

spt_feature_select <- fsi(
  task = task_gpp,
  learner = lrn_ranger_untuned,
  resampling = spt_cv_plan,
  measures = perf_msr,
  terminator = fs_term
)

# Identify optimal feature set for each resampling strategy and store
# Time ~ 10-12 hours

set.seed(123)

progressr::with_progress(
  spt_feature_set <- fs_method$optimize(spt_feature_select)
)

write_rds(spt_feature_select, here("results", "rds", "spt_feature_selector.rds"))
write_rds(spt_feature_set, here("results", "rds", "spt_features.rds"))
rm(spt_feature_select, spt_feature_set)
gc()
# pushoverr::pushover("Feature selection complete")


#' 
#' # Hyper-parameter optimisation (model tuning)
#' 
#' We create a new task with the variables selected from the
#' previous step. Then we define a new learner allowing for tuning of
#' hyperparameters. Finally, we pass the learner and task to an auto-tuner
#' object. We use a random search to identify hyperparameters that 
#' give decent performance.
#' 
## ----rf_tuning_prep, include = FALSE----------------------------------------------------------------------------------------------------------------------------------

message("Beginning hyperparameter tuning\n")

# Load feature set and compare best performance

spt_feature_set <- read_rds(here("results", "rds", "spt_features.rds"))

spt_feature_select <- read_rds(here("results", "rds", "spt_feature_selector.rds"))

spt_feature_set

# Create new task using features selected above
task_gpp_spt <- task_gpp$clone()
task_gpp_spt$select(unlist(spt_feature_set$features))

# Define new learner for tuning

lrn_ranger_tuned_spt <- lrn(
  "regr.ranger", 
  predict_type = "response",
  quantreg = TRUE,
  keep.inbag = TRUE,
  importance = "permutation",
  num.trees = 1001,
  mtry.ratio = to_tune(p_dbl(0, 1)),
  min.node.size = to_tune(p_int(100, 10000, logscale = TRUE)),
  sample.fraction = to_tune(p_dbl(0.1, 0.9))
)

set_threads(lrn_ranger_tuned_spt, 2)

# Define auto-tuner objects for each model

random_tuner <- tnr("random_search")

perf_msr <- msr("regr.rmse")

at_spt <- auto_tuner(
  tuner = random_tuner,
  learner = lrn_ranger_tuned_spt,
  resampling = spt_cv_plan,
  measure = perf_msr,
  term_evals = 25
)

#' 
#' 
#' 
## ----rf_tuning, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
# Spatial model

set.seed(789)

rf_tuned_spt <- at_spt$train(task_gpp_spt)

write_rds(rf_tuned_spt, here("results", "rds", "rf_tuned_spt.rds"))

#' 
#' Performance assessment of spatial and non-spatial models, including code for
#' Predicted vs. Observed GPP (Figure 2).
#' 
## ----performance------------------------------------------------------------------------------------------------------------------------------------------------------

rf_tuned_spt <- read_rds(here("results", "rds", "rf_tuned_spt.rds"))

# Measures to evaluate model performance
measures <- msrs(c("regr.mae", "regr.rmse", "regr.rsq"))

rf_predictions_spt <- rf_tuned_spt$predict(task_gpp_spt)

rf_predictions_spt$score(measures)

# Predicted quantile values
rf_predictions_spt_qu <- predict(rf_tuned_spt$learner$model, task_gpp_spt$data(), type = "quantiles", quantiles = 0.9)

data_with_quantiles_spt <- task_gpp_spt$data() %>%
  bind_cols(rf_predictions_spt_qu$predictions) %>%
  rename(quantile_pred = "quantile= 0.9") %>%
  mutate(quantile_pred = quantile_pred / 1000,
         GPP = GPP / 1000)

#' Performance plot (Figure 2(b))
density_breaks <- c(1, 10, 100, 1000, 10000)

density_plot <- ggplot(data_with_quantiles_spt, aes(x = quantile_pred, y = GPP)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(direction = 1, trans = "log", breaks = density_breaks) +
  geom_abline(slope = seq(0.2, 1, 0.2), intercept = 0, colour = "grey", lwd = 0.8, linetype = "longdash") +
  geom_abline(slope = 1, intercept = 0, colour = "grey", lwd = 1.6) +
  theme_classic() +
  theme(legend.position = c(0.18, 0.75), legend.key.height = unit(1.5, "cm"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent")) +
  xlim(-0.01, 15) +
  ylim(-0.01, 15) +
  labs(x = expression(atop("Potential GPP", (kg~C~m^-2~yr^-1))),
       y = expression(atop("Actual GPP", (kg~C~m^-2~yr^-1))),
       fill = "Number of points")

ggsave(
  here("results", "figures", "rf_density_plot.png"),
  density_plot,
  width = 24, height = 24, units = "cm", dpi = 250
)


#' 
#' 
#' Variable importance plot (Figure 2(a)):
#' 
## ----importance-------------------------------------------------------------------------------------------------------------------------------------------------------

# Variable importance

importance_spt <- rf_tuned_spt$learner$importance()

static_vars <- c("mean_ppt", "slope", "distToRiver", "landform", "twi", "sand", "wriTCFrac")

var_labels <- c(
  mean_ppt = "Mean annual precipitation",
  tMean = "Mean annual temperature",
  pptIntensity = "Precipitation intensity",
  pptMeanDay_anomaly = "Anomaly in\nmean precipitation day",
  parMean = "Mean PAR",
  ppt_anomaly = "Annual precipitation anomaly",
  potentialET = "Potential evapotranspiration",
  ugi = "Unranked Gini index",
  sand = "Soil sand fraction",
  wriTCFrac = "Tree cover fraction",
  slope = "Slope",
  distToRiver = "Distance to river",
  landform = "Landform",
  twi = "Topographic Wetness Index"
)

var_labels_df <- data.frame(variable = names(var_labels), label = unname(var_labels))

importance_spt_df <- tibble(variable = names(importance_spt),
                            importance = as.vector(importance_spt)) %>%
  mutate(type = ifelse(variable %in% static_vars, "Static", "Dynamic")) %>%
  left_join(var_labels_df)

importance_plot <- ggplot(
  importance_spt_df,
  aes(x = importance, y = reorder(label, importance), colour = type)) +
  geom_point(size = 6) +
  theme_bw() +
  scale_colour_manual(values = c("darkblue", "red")) +
  labs(x = "Variable importance", y = "Variable", colour = "Variable type") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

ggsave(here("results", "figures", "importance_plot.png"),
       importance_plot,
       width = 20, height = 20, units = "cm", dpi = 250)
  

#' 
#' 
