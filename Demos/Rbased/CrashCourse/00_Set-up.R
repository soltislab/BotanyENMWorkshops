# Set-up for Botany 2020

# Install R packages prior to the Fridays workshop
## Make a list of packages
list_of_packages <- c("tidyverse", "ridigbio", "spocc", 
					 "scrubr", "RCurl", "rjson", 
					 "sp", "raster", "maptools", 
					 "maps",  "rgdal", "mapproj", 
					 "ggplot2", "ENMeval", "dismo", 
					 "RStoolbox", "hypervolume", "gtools",
					 "caret", "phytools", "factoextra",
					 "FactoMiner", "alphanull", "picante", 
					 "devtools")

## If you do not have these packahes installed, install the package
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Make sure all packages load
lapply(list_of_packages, require, character.only = TRUE)

# Install packages directly from github
library(devtools)

## ggbiplot
install_github("vqv/ggbiplot")
### If you have R v.3.6
#install_github("vqv/ggbiplot", ref = "experimental")
### If you are having issues installing ggbiplot, you may have to uninstall the package 'backports' and reset R. 

## ENMTools
library(devtools)
install_github("danlwarren/ENMTools")
library(ENMTools)

# If ENMTools does not install 
# library(devtools)
# install_github(repo = "ecospat/ecospat/ecospat", ref="3.1")
# library(ecospat)
# install_github("danlwarren/ENMTools")
## If you are still having issues, you may need to reset R and try again. 

## ENMeval
library(devtools)
install_github("bobmuscarella/ENMeval@master", force = TRUE)
## rSDM
devtools::install_github("Pakillo/rSDM", force = TRUE)
## spThin
library(devtools)
install_github('johnbaums/rmaxent',"spThin")