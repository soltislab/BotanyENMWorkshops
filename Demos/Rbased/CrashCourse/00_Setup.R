# Set-up for Botany 2021

# Install R packages prior to the Fridays workshop
## Make a list of packages
list_of_packages <- c("dplyr", 
                      "tidyr",
                      "plyr", 
                      "spocc", 
                      "ridigbio",
                      "tibble", 
                      "tidyverse",
                      "rbison",
                      "CoordinateCleaner",
                      "lubridate",
                      "ggplot2",
                      "gtools",
                      "raster", 
                      "sp", 
                      "spatstat", 
                      "spThin", 
                      "fields", 
                      "ggspatial", 
                      "rgdal", 
                      "rangeBuilder", 
                      "sf", 
                      "dismo", 
                      "devtools", 
                      "ENMeval", 
                      "caret", 
                      "usdm", 
                      "stringr", 
                      "factoextra", 
                      "FactoMineR", 
                      "multcompView", 
                      "ggsci",
                      "gridExtra", 
                      "ecospat", 
                      "rJava", 
                      "viridis", 
                      "ENMTools", 
                      "ape", 
                      "RStoolbox", 
                      "hypervolume", 
                      "phytools",
                      "picante")


## If you do not have these packahes installed, install the package
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Check version of packages
### Compare to version demo was made with
versiondf <- read.csv("data/setup/setup.csv", stringsAsFactors = FALSE)
current_version <- as.character(do.call(c, lapply(list_of_packages, packageVersion)))
versiondf$current_version <- current_version
updatelist <- versiondf[which(versiondf$version != versiondf$current_version), ]
### Update packages with old versions
lapply(as.character(updatelist$packages), install.packages, ask = FALSE)

## Make sure all packages load
lapply(list_of_packages, require, character.only = TRUE)

# Install packages not in CRAN
library(devtools)
install_github('johnbaums/rmaxent')
install_github("marlonecobos/kuenm")


##################################################################################
# Debugging 

# spatstat, spThin, and feilds
## for spatstat, for OS make sure xcode is installed
## in terminal - xcode-select --install
## do not compile from source for spThin or feilds

## ENMTools and/or Java issues - please check the html file for more help!
