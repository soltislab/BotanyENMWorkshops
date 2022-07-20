# ENM_Processing.R
# Processing ENMs
## Script by Anthony Melton and ML Gaynor.

# Load Packages
library(raster)
library(gtools)
library(dplyr)
library(ENMTools)
library(ENMeval)
library(ape)
library(RStoolbox)
library(hypervolume)
library(phytools)

# Load functions
## The following command generates a binary predicted occurrence map. This was written by Anthony Melton
source("functions/Functions_AEM.R")

# Read in models you generated with the MaxEnt GUI into R.
sp1_enm.mx.b <- raster("data/Ecological_Niche_Modeling/enm_output/Galax_urceolata_avg.asc")
sp2_enm.mx.b <- raster("data/Ecological_Niche_Modeling/enm_output/Pyxidanthera_barbulata_avg.asc")
sp3_enm.mx.b <- raster("data/Ecological_Niche_Modeling/enm_output/Pyxidanthera_brevifolia_avg.asc")
sp4_enm.mx.b <- raster("data/Ecological_Niche_Modeling/enm_output/Shortia_galacifolia_avg.asc")

# Read in Bioclim layers
## These steps are described in the previous script
system("rm -rf data/climate_processing/PresentLayers/all/maxent.cache/")
list <- list.files("data/climate_processing/PresentLayers/all/", 
                   full.names = TRUE, recursive = FALSE) 
list <- mixedsort(sort(list))
allstack <- stack(list)
projection(allstack) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"


# ENM Breadth 
## Niche breadth is the breadth of environmental factors for a species' niche, 
## ranging from 0 to 1. 
## When breadth is closer to 1 the more generalist species with wider tolerances.
## Values closer to 0 indicate a more specialized species
## The raster.breadth command in ENMTools measures the smoothness of suitability
## scores across a projected landscape. The higher the score, the more of the 
## available niche space a species occupies.

## MaxEnt GUI
sp1_breadth <- ENMTools::raster.breadth(x = sp1_enm.mx.b)
sp1_breadth$B2
sp2_breadth <- ENMTools::raster.breadth(x = sp2_enm.mx.b)
sp2_breadth$B2
sp3_breadth <- ENMTools::raster.breadth(x = sp3_enm.mx.b)
sp3_breadth$B2
sp4_breadth <- ENMTools::raster.breadth(x = sp4_enm.mx.b)
sp4_breadth$B2

# ENM Overlap
## Calculating niche overlap, Schoener’s D, with ENMEval - 
## Schoener’s D ranges from 0 to 1,
## where zero represents no similarity between the projections 
## and one represents completely identical projections.

## MaxEnt GUI
### Stack the projections and make sure the stack is named.
enm_stack.b <- stack(sp1_enm.mx.b, sp2_enm.mx.b, sp3_enm.mx.b, sp4_enm.mx.b)
names(enm_stack.b) <- c("Galax urceolata", "Pyxidanthera barbulata",  "Pyxidanthera brevifolia","Shortia galacifolia" )

### Calculate 
calc.niche.overlap(enm_stack.b, overlapStat = "D")

####################################################################################################################################################################################
# Phylogenetic correlations
## Let's look at niche overlap over a tree!
## Load the tree file
tree <- ape::read.tree(file = "data/Ecological_Niche_Modeling/AOC_Test_Demo/diapensiaceae_subset.tre")
tree <- ape::compute.brlen(phy = tree, method = "grafen")
plot(tree)

## Drop the outgroup
tree <- drop.tip(tree, "Cyrilla_racemiflora")

# Load data file
alldf <- read.csv("data/cleaning_demo/maxent_ready/diapensiaceae_maxentready_20220712.csv")

## Subset for each species
Galax_urceolata <- dplyr::filter(alldf, name == "Galax urceolata")
Pyxidanthera_barbulata <- dplyr::filter(alldf, name == "Pyxidanthera barbulata")
Pyxidanthera_brevifolia <- dplyr::filter(alldf, name == "Pyxidanthera brevifolia")
Shortia_galacifolia <- dplyr::filter(alldf, name == "Shortia galacifolia")

## Generate species objects for each tree member!
sp1 <- enmtools.species(species.name = "Galax_urceolata",
                        presence.points = Galax_urceolata[,3:2])
sp1$range <- background.raster.buffer(sp1$presence.points, 25000, mask = allstack)
sp1$background.points = background.points.buffer(points = sp1$presence.points, radius = 5000, n = 10000, mask =  allstack[[1]])
##############################
sp2 <- enmtools.species(species.name = "Pyxidanthera_barbulata",
                                        presence.points = Pyxidanthera_barbulata[,3:2])
sp2$range <- background.raster.buffer(sp2$presence.points, 25000, mask = allstack)
sp2$background.points = background.points.buffer(points = sp2$presence.points, radius = 5000, n = 10000, mask =  allstack[[1]])
#############################
sp3 <- enmtools.species(species.name = "Pyxidanthera_brevifolia",
                                        presence.points = Pyxidanthera_brevifolia[,3:2])
sp3$range <- background.raster.buffer(sp3$presence.points, 25000, mask = allstack)
sp3$background.points = background.points.buffer(points = sp3$presence.points, radius = 5000, n = 10000, mask =  allstack[[1]])
#############################
sp4 <- enmtools.species(species.name = "Shortia_galacifolia",
                                          presence.points = Shortia_galacifolia[,3:2])
sp4$range <- background.raster.buffer(sp4$presence.points, 25000, mask = allstack)
sp4$background.points = background.points.buffer(points = sp4$presence.points, radius = 5000, n = 10000, mask =  allstack[[1]])
  

## Create "clade" object with all the species in the tree
clade=enmtools.clade(species = list(sp1, sp2, sp3, sp4), tree = tree)
check.clade(clade)


## Age-Range Correlation Test
range.aoc <- enmtools.aoc(clade = clade,  nreps = 10, overlap.source = "range")
summary(range.aoc)
plot(range.aoc)

## Age-Overlap Correlation Test
glm.aoc <- enmtools.aoc(clade = clade,  env = allstack, nreps = 10, overlap.source = "glm")
summary(glm.aoc)
plot(glm.aoc)

# The results here are pretty meaningless since we're looking at very few, distantly related species, but it serves
# the purpose of the demo. For range AOCs, an intercept >0.5 and negative slope are indicative of sympatric species,
# while an intercept of <0.5 and a positive slope are indicative on non-sympatric speciation. A low intercept and 
# positive slope for niche overlap would indicate niche divergence.
####################################################################################################################################################################################

# Hypervolume
# this will generate a binary map with any cell that contains a suitability score greater than
# or equal to the lowest score with an occurrence point as a presence. There are other ways
# to set the threshold, including percentiles or training statistics.


## Create the binary maps 
sp1.dist <- make.binary.map(model = sp1_enm.mx.b, occ.dat = Galax_urceolata)
sp2.dist <- make.binary.map(model = sp2_enm.mx.b, occ.dat = Pyxidanthera_barbulata)
sp3.dist <- make.binary.map(model = sp3_enm.mx.b, occ.dat = Pyxidanthera_brevifolia)
sp4.dist <- make.binary.map(model = sp4_enm.mx.b, occ.dat = Shortia_galacifolia)

### Plot
par(mfrow = c(2,2),  mar = c(1,1,1,1))
plot(sp1.dist)
plot(sp2.dist)
plot(sp3.dist)
plot(sp4.dist)

## Next, let's work on getting some data from the predicted distributions!
## Niche space can be thought of as a multi-dimensional hypervolume. We're 
## using climatic data in this case, so we're measuring the hypervolume
## of climatic niche space occupied by these species. 
## Warning: This takes forever!
#sp1.hv <- get_hypervolume(binary_projection = sp1.dist, envt = allstack)
#sp2.hv <- get_hypervolume(binary_projection = sp2.dist, envt = allstack)
#sp3.hv <- get_hypervolume(binary_projection = sp3.dist, envt = allstack)
#sp4.hv <- get_hypervolume(binary_projection = sp4.dist, envt = allstack)

## Compare the hypervolumes
#hv_set <- hypervolume_set(hv1 = sp1.hv, hv2 = sp2.hv, check.memory = F)
#hypervolume_overlap_statistics(hv_set)
#plot(hv_set)
#get_volume(hv_set)

