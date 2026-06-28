# 00_Setup.R
# ---------------------------------------------------------------
# Purpose:
#   - Install all R packages required for ENM workflows
#   - Check that all packages load successfully
#   - Install any needed GitHub packages
# ---------------------------------------------------------------

## ---- Define CRAN Packages ----

list_of_packages <- c(
  # Core ENM tools
  "terra", "dismo", "ENMeval", "ENMTools", "biomod2",

  # Spatial / Mapping
  "sf", "rnaturalearth", "ggspatial", "leaflet", "fields", "rangeBuilder",

  # Visualization
  "ggplot2", "viridis", "gridExtra",

  # Data manipulation
  "dplyr", "tidyr", "stringr", "gtools",

  # Phylogenetics
  "ape", "phytools",

  # Biodiversity data
  "ridigbio", "gatoRs",

  # Other utilities
  "multcompView", "usdm", "predicts", "rJava"
)

## ---- Install Missing CRAN Packages ----

# Identify any packages from the list not currently installed.
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

if (length(new.packages)) {
  cat("\nInstalling missing CRAN packages:\n")
  print(new.packages)
  install.packages(new.packages)
} else {
  cat("\nAll required CRAN packages already installed.\n")
}

## ---- Load All CRAN Packages ----

# Try loading all required CRAN packages.
loaded <- sapply(list_of_packages, require, character.only = TRUE)

if (any(!loaded)) {
  cat("\nWARNING: These CRAN packages failed to load:\n")
  print(list_of_packages[!loaded])
} else {
  cat("\nAll CRAN packages loaded successfully.\n")
}

## ---- Install GitHub Packages ----

# devtools is required for installing from GitHub.
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}

# List of GitHub repositories for required packages.
github_repos <- c(
  "vqv/ggbiplot"
)

# Check whether these GitHub packages are installed, and install if not.
for (repo in github_repos) {
  pkg <- tail(strsplit(repo, "/")[[1]], 1)
  if (!pkg %in% installed.packages()[,"Package"]) {
    cat(paste("\nInstalling GitHub package:", repo, "\n"))
    devtools::install_github(repo, upgrade = "never")
  } else {
    cat(paste("\nGitHub package already installed:", pkg, "\n"))
  }
}

# Try loading GitHub packages.
github_packages <- c("ggbiplot")
github_loaded <- sapply(github_packages, require, character.only = TRUE)

if (any(!github_loaded)) {
  cat("\nWARNING: These GitHub packages failed to load:\n")
  print(github_packages[!github_loaded])
} else {
  cat("\nAll GitHub packages loaded successfully.\n")
}
