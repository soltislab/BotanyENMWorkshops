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


# Load the tree file
tree <- ape::read.tree(file = "data/Ecological_Niche_Modeling/AOC_Test_Demo/diapensiaceae_subset.tre")
tree <- ape::compute.brlen(phy = tree, method = "grafen")
plot(tree)

## Drop the outgroup
tree <- drop.tip(tree, "Cyrilla_racemiflora")


# Read ENM models
sp1_enm.mx.b <- raster("data/Ecological_Niche_Modeling/enm_output/Galax_urceolata_avg.asc")
sp2_enm.mx.b <- raster("data/Ecological_Niche_Modeling/enm_output/Pyxidanthera_barbulata_avg.asc")
sp3_enm.mx.b <- raster("data/Ecological_Niche_Modeling/enm_output/Pyxidanthera_brevifolia_avg.asc")
sp4_enm.mx.b <- raster("data/Ecological_Niche_Modeling/enm_output/Shortia_galacifolia_avg.asc")

# Reclassify rasters
reclassify_raster <- function(OGraster){
  ## Reclassify the raster by the suitability score
  OGraster[OGraster >= 0.25] <- 1
  OGraster[OGraster < 0.25] <- 0
  OGraster[is.na(OGraster)] <- 0
  crs(OGraster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  return(OGraster)
}

sp1_enm.mx.b  <- reclassify_raster(sp1_enm.mx.b)
sp2_enm.mx.b  <- reclassify_raster(sp2_enm.mx.b)
sp3_enm.mx.b  <- reclassify_raster(sp3_enm.mx.b)
sp4_enm.mx.b  <- reclassify_raster(sp4_enm.mx.b)

## Inspect the reclassified layers
par(mfrow=c(2,2))
plot(sp1_enm.mx.b)
plot(sp2_enm.mx.b)
plot(sp3_enm.mx.b)
plot(sp4_enm.mx.b)


# Stack rasters and create dataframe
enm_stack <- stack(sp1_enm.mx.b, sp2_enm.mx.b, sp3_enm.mx.b, sp4_enm.mx.b)
names(enm_stack) <- c("Galax_urceolata", "Pyxidanthera_barbulata",  "Pyxidanthera_brevifolia","Shortia_galacifolia" )

## Convert to dataframe
rasterstack_df <- as.data.frame(enm_stack, xy = TRUE)

# Extract Ecoregions for points
## Load USA shapefile
### Obtained at https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states
shp <- ("data/phylogenetic_diversity/USraster/us_eco_l4/us_eco_l4_no_st.shp")
USAshape <- readOGR(shp, layer="us_eco_l4_no_st")

## Correct CRS
newcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
USAshape2 <- spTransform(USAshape, newcrs)


## Convert dataframe to points
df.sp <- rasterstack_df[,1:2]
coordinates(df.sp)<-~x+y
crs(df.sp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Extract points from shapefile
### Warning: this will take a while!! 
df.sp.info <- sp::over(df.sp, USAshape2)

## Bind dataframe with points and with ecosystems
rasterstack_df.info <- cbind(df.sp.info, rasterstack_df)
rasterstack_df.info2 <- rasterstack_df.info %>%
                        dplyr::group_by(L4_KEY) %>% 
                        dplyr::summarize(Galax_urceolata = max(Galax_urceolata),
                                         Pyxidanthera_barbulata = max(Pyxidanthera_barbulata),
                                         Pyxidanthera_brevifolia = max(Pyxidanthera_brevifolia), 
                                         Shortia_galacifolia = max(Shortia_galacifolia))
rasterstack_df.info2df <- as.data.frame(rasterstack_df.info2, row.names = FALSE)
## Removing NA
rasterstack_df.info4 <- rasterstack_df.info2df[1:137,]
## Removing the first column
rasterstack_df.info3 <- rasterstack_df.info4[,-1]
## Renaming the rows
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
  geom_point(data = pd_result, aes(x = rownames(pd_result), y = pd.obs, col = pd.obs))+
  ggtitle("PD values by ecoregion")+
  xlab("Community")+
  ylab("PD")+
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1))

## Plot MPD results by community
ggplot()+
  geom_point(data = mpd_result, aes(x = rownames(mpd_result), y = mpd.obs, col = mpd.obs))+
  ggtitle("MPD values by ecoregion")+
  xlab("Community")+
  ylab("MPD")+
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1))

## Plot MNTD results by community
ggplot()+
  geom_point(data = mntd_result, aes(x = rownames(mntd_result), y = mntd.obs, col = mntd.obs))+
  ggtitle("MNTD values by ecoregion")+
  xlab("Community")+
  ylab("MNTD")+
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1))
