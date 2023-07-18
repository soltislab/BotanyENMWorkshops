# Set-up for Botany 2023
## 00_Setup.R
## ML Gaynor

# Install R packages
## Make a list of packages
list_of_packages <- c("ade4","ape", "biomod2",
                      "caret", "CoordinateCleaner", 
                      "devtools","dismo", 
                      "dplyr", "ecospat", 
                      "ENMeval",  "ENMTools", 
                      "fields",  "gatoRs",  
                      "ggplot2", "ggspatial",
                      "gridExtra",  "gtools", 
                      "hypervolume",  "kuenm", 
                      "lattice", "leaflet",
                      "magrittr",   "maps",   
                      "multcompView", "parsedate", 
                      "phytools", "plyr", 
                      "rangeBuilder", "raster",
                      "Rcpp",  "ridigbio",
                      "rJava","rgbif",
                      "sf",   "sp",
                      "spThin","spam",
                      "spatstat.geom","stringr",
                      "terra", "tidyr", 
                      "usdm",  "utils", 
                      "viridis", "viridisLite")

## Here we identify which packages are not installed and install these for you
### Please do install all package which need compilation
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Check version of packages
### Compare to package version used for demo creation
### This data frame has two columns: (1) package (2) version, indicating which
### version of the package was used to develop this workshop
versiondf <- read.csv("data/setup/setup.csv", stringsAsFactors = FALSE)
### Save your version under a new column named "current_version"
versiondf$current_version <- as.character(do.call(c, lapply(list_of_packages, packageVersion)))
### Compare the version the workshop was developed with to your current version
updatelist <- versiondf[which(versiondf$version != versiondf$current_version), ]
### Update packages with old versions
lapply(as.character(updatelist$packages), install.packages, ask = FALSE)

## Make sure all packages load
## Packages that are FALSE did not load
### If anything prints here, then something went wrong and a package did not install
loaded <- lapply(list_of_packages, require, character.only = TRUE)
list_of_packages[loaded == FALSE]

# Install packages not in CRAN
library(devtools)
install_github('johnbaums/rmaxent')
install_github("marlonecobos/kuenm")

## Check and make sure all github packages load
github_packages <- c("rmaxent", "kuenm", "gatoRs")
github_loaded <- lapply(github_packages, require, character.only = TRUE)
### If anything prints here, then something went wrong and a package did not install
github_packages[github_loaded == FALSE]

##################################################################################
# Debugging 

# spatstat, spThin, and feilds
## for spatstat, for OS make sure xcode is installed
## in terminal - xcode-select --install
## do not compile from source for spThin or feilds

## ENMTools and/or Java issues - please check the html file for more help!
