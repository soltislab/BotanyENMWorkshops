# Phylogenetic Diversity
## Phylogenetic_Diversity.R
## Original script by J. Janzten. 
## Modified by ML Gaynor.

# Load Packages
library(raster)
library(ape)
library(phytools)
library(picante)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(rgdal)
library(rgeos)
library(sp)
library(spatialEco)


# Load the tree file
tree <- read.tree(file = "data/phylogenetic_diversity/tree_withBL.tre")
tree <- ape::compute.brlen(phy = tree, method = "grafen")
plot(tree)

# Read ENM models
Asclepias_curtissii_enm.mx.b <- raster("data/enm_output/Asclepias_curtissii/Asclepias_curtissii_avg.asc")
Asimina_obovata_enm.mx.b <- raster("data/enm_output/Asclepias_obovata/Asimina_obovata_avg.asc")
Pinus_palustris_enm.mx.b <- raster("data/enm_output/Pinus_palustris/Pinus_palustris_avg.asc")
Liatris_chapmanii_enm.mx.b <- raster("data/phylogenetic_diversity/Liatris_chapmanii/Liatris_chapmanii_avg.asc")

# Reclassify rasters
reclassify_raster <- function(OGraster){
  ## Reclassify the raster by the suitability score
  OGraster[OGraster >= 0.25] <- 1
  OGraster[OGraster < 0.25] <- 0
  OGraster[is.na(OGraster)] <- 0
  crs(OGraster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  return(OGraster)
}

Asclepias_curtissii_enm.mx.b  <- reclassify_raster(Asclepias_curtissii_enm.mx.b)
Asimina_obovata_enm.mx.b <- reclassify_raster(Asimina_obovata_enm.mx.b)
Pinus_palustris_enm.mx.b <- reclassify_raster(Pinus_palustris_enm.mx.b)
Liatris_chapmanii_enm.mx.b <- reclassify_raster(Liatris_chapmanii_enm.mx.b)

## Inspect the reclassified layers
par(mfrow=c(2,4))
plot(Asclepias_curtissii_enm.mx.b)
plot(Asimina_obovata_enm.mx.b)
plot(Pinus_palustris_enm.mx.b)
plot(Liatris_chapmanii_enm.mx.b)
par(mfrow=c(1,1))

# Stack rasters and create dataframe
rasterstack <- stack(Asclepias_curtissii_enm.mx.b, Asimina_obovata_enm.mx.b,
                   Pinus_palustris_enm.mx.b, Liatris_chapmanii_enm.mx.b)
## Rename
names(rasterstack) <- c("Asclepias_curtissii1", "Asimina_obovata1", "Pinus_palustris1", "Liatris_chapmanii1")

## Convert to dataframe
rasterstack_df <- as.data.frame(rasterstack, xy = TRUE)

# Extract Ecoregions for points
## Load Florida shapefile
### Obtained at https://www.epa.gov/eco-research/ecoregion-download-files-state-region-4
shp <- ("data/phylogenetic_diversity/Florida_Raster/reg4_eco_l4/reg4_eco_l4.shp")
floridashape <- shapefile(shp)
## Correct CRS
newcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
floridashape2 <- spTransform(floridashape, newcrs)
## Crop to florida
floridashape3 <- raster::crop(floridashape2, extent(-87.625, -80.04167, 24.54167, 31))
plot(floridashape3)


## Convert dataframe to points
df.sp <- rasterstack_df[,1:2]
coordinates(df.sp)<-~x+y
crs(df.sp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Extract points from shapefile
df.sp.info <- sp::over(df.sp, floridashape3)

## Bind dataframe with points and with ecosystems
rasterstack_df.info <- cbind(df.sp.info, rasterstack_df)
rasterstack_df.info2 <- rasterstack_df.info %>%
                        dplyr::group_by(L4_KEY) %>% 
                        dplyr::summarize(Asclepias_curtissii = max(Asclepias_curtissii1),
                               Asimina_obovata = max(Asimina_obovata1),
                               Pinus_palustris = max(Pinus_palustris1), 
                               Liatris_chapmanii = max(Liatris_chapmanii1))
rasterstack_df.info2df <- as.data.frame(rasterstack_df.info2, row.names = FALSE)
rasterstack_df.info4 <- rasterstack_df.info2df[1:20,]
rasterstack_df.info3 <- rasterstack_df.info4[,-1]
row.names(rasterstack_df.info3) <- NULL
row.names(rasterstack_df.info3) <- rasterstack_df.info4[,1]
rasterstack_df.info3

                              
                               
####################################################################################################################################################################################
# Phylogenetic Diversity Calculations
# Match tree and dataframe
matched <- match.phylo.comm(tree, rasterstack_df.info3)

#Make object for tree including matching taxa only (non-matching taxa will be pruned out)
matchtree <- matched$phy
#Make object for dataset including matching taxa only (non-matching taxa removed)
matchcomm <- matched$comm 

# Calculate phylogenetic distance matrix for use in MPD and MNTD calculations
## MPD is mean pairwise distance, or mean phylogenetic distance among all pairs of 
## species within a community. MNTD is mean nearest taxon distance, 
## or the mean distance between each species within a community and its closest relative.
phydist <- cophenetic(matched$phy)

# Calculate indices using picante
## Null model options include taxa.labels (shuffle tips of phylogeny) among 
## others Use ses.pd, ses.mpd and ses.mntd for more info on alternative null
## models and other parameters and output format. 
## Number of runs - only using 99 runs (comparing to null model) due to time constraints.

## PD calculates Faith’s PD and standard effect size of PD
### PD is phylogenetic distance
pd_result <-ses.pd(matchcomm, matchtree, null.model="taxa.labels", runs=99)

## Calculates mpd and ses mpd - equivalent to -NRI
### NRI is net-relatedness index, standardizes MPD or mean pairwise distance. 
### ses stands for standardized effect size.
### If we had abundance data, we could include that for mpd and mntd. 
### The default is abundance.weighted = FALSE.
mpd_result <- ses.mpd(matchcomm, phydist, null.model="taxa.labels", abundance.weighted = FALSE, runs=99)

## Calculates mntd and ses mntd - equivalent to -NTI
## NTI or nearest taxon index, standardizes MNTD mean nearest taxon distance. 
## ses stands for standardized effect size.
mntd_result <- ses.mntd(matchcomm, phydist, null.model="taxa.labels", runs=99)

# Look at results
## Positive ses values (obs.z values) and p values > 0.95 indicate overdispersion (greater than expected).
## Negative ses values (obs.z values) and p values < 0.05 indicate clustering (less than expected).
## Values not significantly different from zero indicate taxa in community randomly distributed across tree.

pd_result
mpd_result
mntd_result

## Plot PD results by community
ggplot()+
  geom_point(data = pd_result, aes(x = rownames(pd_result), y = pd.obs))+
  ggtitle("PD values for FL by ecoregion")+
  xlab("Community")+
  ylab("PD")+
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1))

## Plot MPD results by community
ggplot()+
  geom_point(data = mpd_result, aes(x = rownames(mpd_result), y = mpd.obs))+
  ggtitle("MPD values for FL by ecoregion")+
  xlab("Community")+
  ylab("MPD")+
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1))

## Plot MNTD results by community
ggplot()+
  geom_point(data = mntd_result, aes(x = rownames(mntd_result), y = mntd.obs))+
  ggtitle("MNTD values for FL by ecoregion")+
  xlab("Community")+
  ylab("MNTD")+
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1))
