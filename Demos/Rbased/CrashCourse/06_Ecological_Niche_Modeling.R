# Ecological_Niche_Modeling.R
## Ecological Niche Modeling 
## Script by Anthony Melton and ML Gaynor.  

# This script is for generating and testing ENMs using ENMEval. Please see 
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13628
# for the paper describing ENMEval and 
# https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0.0-vignette.html
# for the vignette.

# Set up java memory 
options(java.parameters = "-Xmx16g") # increase memory that can be used

# Load Packages
library(raster)
library(gtools)
library(dplyr)
library(dismo)
library(ENMeval)
library(ggplot2)
library(viridis)
library(kuenm)

# Load Function
source("functions/ENMevaluation.R")

# Load data file
alldf <- read.csv("data/cleaning_demo/maxent_ready/diapensiaceae_maxentready_20230605.csv")

## Subset for each species
Galax_urceolata <- dplyr::filter(alldf, species == "Galax urceolata")
Pyxidanthera_barbulata <- dplyr::filter(alldf, species == "Pyxidanthera barbulata")
Pyxidanthera_brevifolia <- dplyr::filter(alldf, species == "Pyxidanthera brevifolia")
Shortia_galacifolia <- dplyr::filter(alldf, species == "Shortia galacifolia")


# Raster layers
list <- list.files("data/climate_processing/PresentLayers/all", full.names = TRUE, recursive = FALSE) 
list <- mixedsort(sort(list))
allstack <- stack(list)

## Read in species training layers
gstack <- stack(mixedsort(sort(list.files("data/climate_processing/PresentLayers/Galax_urceolata/", full.names = TRUE))))
pbastack <- stack(mixedsort(sort(list.files("data/climate_processing/PresentLayers/Pyxidanthera_barbulata/", full.names = TRUE))))
pbrstack <- stack(mixedsort(sort(list.files("data/climate_processing/PresentLayers/Pyxidanthera_brevifolia/", full.names = TRUE))))
sstack <- stack(mixedsort(sort(list.files("data/climate_processing/PresentLayers/Shortia_galacifolia/", full.names = TRUE))))

## Fix projection
projection(allstack) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
projection(gstack) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
projection(pbastack) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
projection(pbrstack) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
projection(sstack) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

# Model Generation
## dismo Model Generations
## Warning: Time intensive on Apple M1 chip!
### Match MaxEnt GUI settings (see word document and powerpoint)
evaldis <- dismo::maxent(x = gstack, p = Galax_urceolata[, c("longitude", "latitude")], nbg = 10000,
                         args = c("projectionlayers=data/climate_processing/PresentLayers/all",
                                  "responsecurves", "jackknife",  "outputformat=logistic",
                                  "randomseed", "randomtestpoints=25",  "replicates=5", 
                                  "replicatetype=subsample",  "maximumiterations=5000", "writebackgroundpredictions",
                                  "responsecurvesexponent", "writeplotdata"), 
                         removeDuplicates = TRUE
                         ,path = "data/Ecological_Niche_Modeling/enm_output/Galax_urceolata/"
              )


## ENMeval Model generation
#### ENMeval will generate multiple models and test them per the specified test method.
#### Two important variables here are the regularization multiplier value and the feature 
#### class (fm). FC will allow for different shapes in response curves (linear, hinge, 
#### quadratic, product, and threshold) can be used in the model and RM will influence how 
#### many parameters are included in the model.
eval1 <- ENMeval::ENMevaluate(occ = Galax_urceolata[, c("longitude", "latitude")], 
                              env = gstack,
                              tune.args = list(fc = c("L","Q"), rm = 1:2), # test the feature classes L = linear and Q = quadratic
                              partitions = "block",
                              n.bg = 10000,
                              parallel = FALSE,
                              algorithm = 'maxent.jar', 
                              user.eval = proc)

##########################################################################################
# Model Statistics
## dismo
### Inspect the dismo model based on the html
browseURL(evaldis@html)

## ENMeval 
### Inspect the many models
#### Visualize
maps <- eval1@predictions
plot(maps)

#### Look at model overlap
mod_overlap <- calc.niche.overlap(eval1@predictions, overlapStat = "D")
mod_overlap


### Inspect the results
### Identify the best model
#### selecting models with the lowest average test omission rate and 
#### the highest average validation AUC 
results <- eval.results(eval1)
opt.seq <- results %>% 
            dplyr::filter(or.10p.avg == min(or.10p.avg)) %>% 
            dplyr::filter(auc.val.avg == max(auc.val.avg))
opt.seq 

#### Subset model
mod.seq <- eval.models(eval1)[[opt.seq$tune.args]]

#### Inspect
mod.seq@results

### Look at variable contribution
plot(mod.seq)
### Look at the response curves
dismo::response(mod.seq)

## Project model to allstack
p <- predict(mod.seq, allstack) 

### Visualize
### Make p plottable 
p_df <-  as.data.frame(p, xy = TRUE)
### Plot
ggplot() +
  geom_raster(data = p_df, aes(x = x, y = y, fill = layer)) +
  geom_point(data= Galax_urceolata, 
             mapping = aes(x = longitude, y = latitude), 
             col='red', cex=0.05) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours = viridis::viridis(99),
                       na.value = "black")

##########################################################################################
# Save outputs
## R saved dataset
save(mod.seq, file = "data/Ecological_Niche_Modeling/enm_output/ENMeval/GalaxENM.rda")

## Save Raster 
writeRaster(x = p, filename = "data/Ecological_Niche_Modeling/enm_output/ENMeval/GalaxENM.asc",
            format = "ascii", NAFlag = "-9999", overwrite = T)

