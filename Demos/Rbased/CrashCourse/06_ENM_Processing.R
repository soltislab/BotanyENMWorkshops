# ENM_Processing.R
# Processing ENMs
## Script by Anthony Melton and ML Gaynor.

# Load Packages
library(raster)
library(ENMTools)
library(ENMeval)
library(gtools)
library(dplyr)
library(RStoolbox)
library(hypervolume)
library(phytools)
library(alphanull)

# See Ecological_Niche_Modeling.R to generate models in R
# You can also read the models you generated with the MaxEnt GUI into R.
Asclepias_curtissii_enm.mx.b <- raster("data/enm_output/Asclepias_curtissii/Asclepias_curtissii_avg.asc")
Asimina_obovata_enm.mx.b <- raster("data/enm_output/Asclepias_obovata/Asimina_obovata_avg.asc")
Pinus_palustris_enm.mx.b <- raster("data/enm_output/Pinus_palustris/Pinus_palustris_avg.asc")

# You can also read the models you generated with ENMeval or ENMTools in R.
load(file = "data/enm_output/ENMTools/Asclepias_curtissii.rda")
load(file = "data/enm_output/ENMTools/Asimina_obovata.rda")
load(file = "data/enm_output/ENMTools/Pinus_palustris.rda")

# Read in Bioclim layers
## These steps are described in the previoius script
list <- list.files("data/climate_processing/PresentLayers/", full.names = T, recursive = FALSE) 
list <- mixedsort(sort(list))
envtStack <- stack(list)
envt.subset<-subset(envtStack, c(3, 4, 6, 9, 10, 13, 15, 19)) 
envt.subset <- setMinMax(envt.subset)
proj4string(envt.subset) <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# ENM Breadth 
## Niche breadth is the breadth of environmental factors for a species' niche, 
## ranging from 0 to 1. 
## When breadth is closer to 1 the more generalist species with wider tolerances.
## Values closer to 0 indicate a more specialized species
## The raster.breadth command in ENMTools measures the smoothness of suitability
## scores across a projected landscape. The higher the score, the more of the 
## available niche space a species occupies.

## Rbased
Ac_breadth <- ENMTools::raster.breadth(x = Asclepias_curtissii_enm.mx$suitability)
Ac_breadth$B2
Ao_breadth <- ENMTools::raster.breadth(x = Asimina_obovata_enm.mx$suitability)
Ao_breadth$B2
Pp_breadth <- ENMTools::raster.breadth(x= Pinus_palustris_enm.mx$suitability)
Pp_breadth$B2

## MaxEnt GUI
Ac_breadth.b <- ENMTools::raster.breadth(x = Asclepias_curtissii_enm.mx.b)
Ac_breadth.b$B2
Ao_breadth.b <- ENMTools::raster.breadth(x = Asimina_obovata_enm.mx.b)
Ao_breadth.b$B2
Pp_breadth.b <- ENMTools::raster.breadth(x = Pinus_palustris_enm.mx.b)
Pp_breadth.b$B2

# ENM Overlap
## Calculating niche overlap, Schoener’s D, with ENMEval - 
## Schoener’s D ranges from 0 to 1,
## where zero represents no similarity between the projections 
## and one represents completely identical projections.

## ## Rbased
### Stack the projections and make sure the stack is named.
enm_stack <- stack(Asclepias_curtissii_enm.mx$suitability, Asimina_obovata_enm.mx$suitability,Pinus_palustris_enm.mx$suitability)
names(enm_stack) <- c("Asclepias_curtissii", "Asimina_obovata", "Pinus_palustris")

### Calculate 
calc.niche.overlap(enm_stack, stat = "D")

## MaxEnt GUI
### Stack the projections and make sure the stack is named.
enm_stack.b <- stack(Asclepias_curtissii_enm.mx.b, Asimina_obovata_enm.mx.b,Pinus_palustris_enm.mx.b)
names(enm_stack.b) <- c("Asclepias_curtissii", "Asimina_obovata", "Pinus_palustris")

### Calculate 
calc.niche.overlap(enm_stack.b, stat = "D")

####################################################################################################################################################################################
# Hypervolume

# The following command generates a binary predicted occurrence map. To run this command,
# a function must be generated. There are two functions in the Functions.R script.
# By sourcing the script, the functions will be generated and you can run them here!

source("functions/Functions_AEM.R") # Only needs to be ran once per session to generate the functions.

# this will generate a binary map with any cell that contains a suitability score greater than
# or equal to the lowest score with an occurrence point as a presence. There are other ways
# to set the threshold, including percentiles or training statistics.

## Load occurrence data 
florida_points <- read.csv("data/cleaning_demo/MaxEntPointsInput_cleaned.csv")

Asclepias_curtissii.pts <- florida_points %>%
  filter(species == "Asclepias_curtissii")

Asimina_obovata.pts <- florida_points %>%
  filter(species == "Asimina_obovata")

Pinus_palustris.pts <- florida_points %>%
  filter(species == "Pinus_palustris")

## Create the binary maps 

Asimina.dist <- make.binary.map(model = Asimina_obovata_enm.mx.b, occ.dat = Asimina_obovata.pts)
Pinus.dist <- make.binary.map(model = Pinus_palustris_enm.mx.b, occ.dat = Pinus_palustris.pts)
Asclepias.dist <- make.binary.map(model = Asclepias_curtissii_enm.mx.b, occ.dat = Asclepias_curtissii.pts)

### Plot
par(mfrow = c(1,3))
plot(Asimina.dist)
plot(Pinus.dist)
plot(Asclepias.dist)


## Next, let's work on getting some data from the predicted distributions!
## Niche space can be thought of as a multi-dimensional hypervolume. We're 
## using climatic data in this case, so we're measuring the hypervolume
## of climatic niche space occupied by these species in FL.

Asimina.hv <- get_hypervolume(binary_projection = Asimina.dist, envt = envt.subset)
Pinus.hv <- get_hypervolume(binary_projection = Pinus.dist, envt = envt.subset)
Asclepias.hv <- get_hypervolume(binary_projection = Asclepias.dist, envt = envt.subset)

## Compare the hypervolumes
hv_set <- hypervolume_set(hv1 = Asimina.hv, hv2 = Pinus.hv, check.memory = F)
hypervolume_overlap_statistics(hv_set)
plot(hv_set)
get_volume(hv_set)

####################################################################################################################################################################################
# Phylogenetic correlations
## Let's look at niche overlap over a tree!
## Load the tree file
tree <- read.newick(file = "data/AOC_Test_Demo/tree.tre")
tree <- ape::compute.brlen(phy = tree, method = "grafen")
plot(tree)

## Need to add a species to match the tree
Liatris_chapmanii.pts <- read.csv("data/AOC_Test_Demo/clean_L_chapmaniiPointsSPOCC_EDITS.csv")[,2:3]

## Generate species objects; should have first three from previous exercise
Asclepias_curtissii <- enmtools.species(species.name = "Asclepias_curtissii",
                                            presence.points = Asclepias_curtissii.pts[,3:2])
Asclepias_curtissii$range <- background.raster.buffer(Asclepias_curtissii$presence.points, 25000, mask = envt.subset)
Asclepias_curtissii$background.points = background.points.buffer(points = Asclepias_curtissii$presence.points, radius = 5000, n = 10000, mask =  envt.subset[[1]])
##############################
Asimina_obovata <- enmtools.species(species.name = "Asimina_obovata",
                                        presence.points = Asimina_obovata.pts[,3:2])
Asimina_obovata$range <- background.raster.buffer(Asimina_obovata$presence.points, 25000, mask = envt.subset)
Asimina_obovata$background.points = background.points.buffer(points = Asimina_obovata$presence.points, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#############################
Pinus_palustris <- enmtools.species(species.name = "Pinus_palustris",
                                        presence.points = Pinus_palustris.pts[,3:2])
Pinus_palustris$range <- background.raster.buffer(Pinus_palustris$presence.points, 25000, mask = envt.subset)
Pinus_palustris$background.points = background.points.buffer(points = Pinus_palustris$presence.points, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#############################
Liatris_chapmanii <- enmtools.species(species.name = "Liatris_chapmanii",
                                          presence.points = Liatris_chapmanii.pts)
Liatris_chapmanii$range <- background.raster.buffer(Liatris_chapmanii$presence.points, 25000, mask = envt.subset)
Liatris_chapmanii$background.points = background.points.buffer(points = Liatris_chapmanii$presence.points, radius = 5000, n = 10000, mask =  envt.subset[[1]])
  

## Create "clade" object with all the species in the tree
clade=enmtools.clade(species = list(Asclepias_curtissii, Pinus_palustris, Asimina_obovata, Liatris_chapmanii), tree = tree)
check.clade(clade)


## Age-Range Correlation Test
range.aoc <- enmtools.aoc(clade = clade,  nreps = 10, overlap.source = "range")
summary(range.aoc)

## Age-Overlap Correlation Test
glm.aoc <- enmtools.aoc(clade = clade,  env = envt.subset, nreps = 10, overlap.source = "glm")
summary(glm.aoc)

# The results here are pretty meaningless since we're looking at very few, distantly related species, but it serves
# the purpose of the demo. For range AOCs, an intercept >0.5 and negative slope are indicative of sympatric species,
# while an intercept of <0.5 and a positive slope are indicative on non-sympatric speciation. A low intercept and 
# positive slope for niche overlap would indicate niche divergence.


