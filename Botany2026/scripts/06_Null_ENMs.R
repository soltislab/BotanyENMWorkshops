# Null_models.R
# ---------------------------------------------------------------
# Purpose: Assess performance of niche models using null models
# Authors: Tyler Radtke and Sebastian Fernandez
# ---------------------------------------------------------------

## ---- Load Required Packages ----
library(ENMeval)
library(ggplot2)
library(dplyr)

## ---- A) Load Cleaned Occurrence Data ----

# Load all occurrence records
alldf <- read.csv("data/02_cleaning/maxent_ready/diapensiaceae_maxentready_2025_06_27.csv")

## ---- B) Load Optimal Model Parameters ----

# Load optimal feature class (fc) and regularization multiplier (rm)
opt.seq <- read.delim("data/05_ENMs/Galax_urceolata_OptModel.txt")

# Extract parameters
fc <- opt.seq$fc
rm <- opt.seq$rm

## ---- C) Load ENMeval Object ----

# Load saved ENMeval object for Galax urceolata
load("data/05_ENMs/Galax_urceolata_ENM_eval.RData")

## ---- D) Run Null Model Simulations ----

# Run ENMnulls with optimal parameters and 100 iterations
spec.mod.null <- ENMnulls(
  eval,
  mod.settings = list(fc = fc, rm = rm),
  no.iter = 100
)

## ---- E) Save Null Model Comparison Results ----

# Extract comparison results between empirical and null models
null_comparison_results <- null.emp.results(spec.mod.null)

# Save results as CSV
write.csv(
  null_comparison_results,
  file = "data/05_ENMs/Galax_urceolata_null_comparison_results.csv",
  row.names = FALSE
)

## ---- F) Plot Null Model Results ----

# Generate histogram plots for AUC and OR.10p statistics
spec.null <- evalplot.nulls(
  spec.mod.null,
  stats = c("or.10p", "auc.val"),
  plot.type = "histogram"
)

# Display plot in R console
plot(spec.null)

# Save plot to PNG
ggsave(
  filename = "Galax_urceolata_null_histogram.png",
  plot = spec.null,
  path = "data/05_ENMs/",
  height = 12,
  width = 13
)
