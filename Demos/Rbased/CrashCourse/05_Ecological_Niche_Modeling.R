# Ecological_Niche_Modeling.R
## Ecological Niche Modeling 
### We will not demo these scripts during the workshop due to system requirements of rJava.
## Script by Anthony Melton and ML Gaynor.  

# This script is for generating and testing ENMs using ENMEval. Please see 
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12261
# for the paper describing ENMEval and 
# https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html#cb2
# for the vignette.

# Set up java memory 
options(java.parameters = "- Xmx16g") # increase memory that can be used

# Load Packages
library(dplyr)
library(ENMeval)
library(rSDM)
library(dismo)
library(spThin)
library(readr)
library(rmaxent)
library(gtools)


# Load Datafiles 
florida_points <- read.csv("data/cleaning_demo/MaxEntPointsInput_cleaned.csv")

## Subset occurence records for each species
Asclepias_curtissii <- florida_points %>%
                       filter(species == "Asclepias_curtissii")

Asimina_obovata <- florida_points %>%
                   filter(species == "Asimina_obovata")

Pinus_palustris <- florida_points %>%
                   filter(species == "Pinus_palustris")

# Raster layers
list <- list.files("data/climate_processing/PresentLayers/", full.names = T, recursive = FALSE) 
list <- mixedsort(sort(list))
envtStack <- stack(list)

## Removing highly correlated layers
library(caret)
c <- data.matrix(read.csv("data/climate_processing/correlationBioclim.csv", header = TRUE, row.names = 1,
                          sep = ","))
### Take absolute value
c <- abs(c)
# Find rows that should be removed
envtCor <- findCorrelation(c, cutoff = 0.80, names = TRUE, exact = TRUE)
sort(envtCor)
# keep: bio2, bio3, bio5, bio8, bio9, bio12,  bio14, and bio18

### Subset layers
envt.subset<-subset(envtStack, c(3, 4, 6, 9, 10, 13, 15, 19)) 
envt.subset

## Thin occurrence records
### If there are still a lot of points after cleaning, you can further reduce to minimize
### biases in g-space. The following code requires the thin_max.R script by Dan Warren available at
### http://enmtools.blogspot.com/

nrow(Pinus_palustris)

if (nrow(Pinus_palustris) > 50) {
source("functions/thin.max.R")
Pinus_palustris <- thin.max(Pinus_palustris, c("long", "lat"), 50) # trim down to N; check column names
}

nrow(Pinus_palustris)
plot(x = Pinus_palustris$lat, y = Pinus_palustris$long, col = "black")

## Designate background data
### Randomly sample 10,000 background points from one background extent raster (only one per cell 
### without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all 
### pixels will be used for background. We will be sampling from the biome variable because it is 
### missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envt.subset[[1]], n=10000)
bg <- as.data.frame(bg)

### Check the bg point distribution in g-space
plot(envt.subset[[1]], legend=FALSE)
plot(x = bg$x, y = bg$y, pch = 16, col='red')

### Model generation
#### ENMeval will generate multiple models and test them per the specified test method.
#### Two important variables here are the regularization multiplier value and the feature 
#### class (fm). FC will allow for different shapes in response curves (linear, hinge, 
#### quadratic, product, and threshold) can be used in the model and RM will influence how 
#### many parameters are included in the model.
modeval <- ENMevaluate(occ = Pinus_palustris[, c("long", "lat")], 
                       env = envt.subset,
                       bg.coords = bg,
                       #categoricals = "Taxonomy",
                       algorithm = 'maxent.jar',
                       #RMvalues = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                       #fc = c("L", "H", "LQ", "LQH", "LQP", "LQPH", "LQPHT"),
                       RMvalues = c(1, 2.5, 5),
                       fc = c("L", "H", "LQH"), 
                       method = "block",
                       aggregation.factor = c(2,2),
                       #overlap = TRUE,
                       clamp = TRUE, 
                       rasterPreds = TRUE,
                       parallel = TRUE,
                       numCores = 4,
                       bin.output = TRUE,
                       progbar = TRUE)

# Evaluating the Models 
## AICc score
### Get the "best" model based on AICc score
aic.opt <- modeval@models[[which(modeval@results$delta.AICc==0)]]
aic.opt
aic.opt@lambdas

## Lambda
### The Parse_lambdas.R script by John Baums can be found at 
### https://github.com/johnbaums/rmaxent/blob/master/R/parse_lambdas.R
source("functions/parse_lambdas.R")
parse_lambdas(aic.opt) 
### features that have a lambda of > 0. These are
### the features that are included in the model.

## Create a table  with the stats for all models
results <- modeval@results
kable(results)

## Visualize model statistics 
### Relative occurrence rate
plot(modeval@predictions[[which(modeval@results$delta.AICc==0)]], main="Relative occurrence rate")

### Variable importance 
(df <- var.importance(aic.opt))
barplot(df$permutation.importance, names.arg=df$variable, las=2, ylab="Permutation Importance")

### delta.AICc
eval.plot(results)

## Visualize Predicted Models
### Plots all of the predictions by the different models
maps <- modeval@predictions
plot(maps)

### Plots the predictions of the best model
plot(maps[[which(results$delta.AICc == 0)]], main = "Models with lowest AICc")

## Response curves
### Plots all of the response curves for the best model
for (i in which(results$delta.AICc == 0)) {
  response(modeval@models[[i]])
}

# Pick "best" model based on AICc
which(results$delta.AICc == 0)
n <- which(results$delta.AICc == 0)

## Isolate the "best" model
modeval@predictions[[n]]  # raw output; Check model name!
p <- predict(modeval@models[[n]], envt.subset) 
plot(p, zlim = c(0,1))  # plots the predictions in a 0-1 scale; logistic output

##########################################################################################
# Save outputs
## R saved dataset
save(modeval, file = "PineENM.rda")
# load(file = "FILENAME.rda")

## Raster 
### Save "best" model 
writeRaster(x = p, filename = "PineENM.asc", format = "ascii", NAFlag = "-9999", overwrite = T)

## Clear memory and get session summary
# rm(list = ls())
# devtools::session_info()
