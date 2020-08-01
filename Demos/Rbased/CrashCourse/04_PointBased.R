# PointBased.R
# Ecological Analysis using Points
# This is based off scripts created by Hannah Owens and James Watling, and Anthony Melton.  
# Modified and created by ML Gaynor. 

# Load Packages
library(raster)
library(rgdal)
library(ENMTools)
library(dismo)
library(RStoolbox)
library(hypervolume)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)
library(ggbiplot)

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

####################################################################################################################################################################################
# ENM - Realized Niche

## For each occurence record, extract the value for each bioclim variables using the package raster.
ptExtract_Asclepias_curtissii <- raster::extract(envt.subset, Asclepias_curtissii[3:2])
ptExtract_Asimina_obovata <- raster::extract(envt.subset, Asimina_obovata[3:2])
ptExtract_Pinus_palustris <- raster::extract(envt.subset, Pinus_palustris[3:2])

## Convert the extracted information for all three species and combined. I wrote a quick function to convert the files to dataframes.
convert_ptExtract <- function(value, name){
                    # convert ptExtract to a data frame
                    value_df <- as.data.frame(value)
                    # add a column with species name
                    value_df_DONE <- value_df %>% 
                      mutate(species = name)
                    # return data frame
                    return(value_df_DONE)
}

ptExtract_Asclepias_curtissii_df <- convert_ptExtract(ptExtract_Asclepias_curtissii, "Asclepias_curtissii")
ptExtract_Asimina_obovata_df <- convert_ptExtract(ptExtract_Asimina_obovata, "Asimina_obovata")
ptExtract_Pinus_palustris_df <- convert_ptExtract(ptExtract_Pinus_palustris, "Pinus_palustris")

### Combined the three converted files and remove any values where NA appears 
pointsamples_combined <- rbind(ptExtract_Asclepias_curtissii_df,ptExtract_Asimina_obovata_df, ptExtract_Pinus_palustris_df)
pointsamples_combined <- pointsamples_combined %>% 
                         drop_na(bio2, bio3, bio5, bio8, bio9, bio12, bio14, bio18)


## PCA 
### Create two dataframes.
data.bioclim <- pointsamples_combined[, 1:8]
data.species <- pointsamples_combined[, 9]

### Using only the bioclim columns to run the principal components analysis.
data.pca <- prcomp(data.bioclim, scale. = TRUE) 

#### Understanding the PCA - Optional 
library(factoextra)
library(FactoMineR)

##### When you use the command prcomp your loading variables show up as rotational variables. 
##### Thanks to a really great answer on stack overflow: https://stackoverflow.com/questions/43407859/how-do-i-find-the-link-between-principal-components-and-raw-datas-variables
##### you can even convert the rotational 
##### variable to show the relative contribution.

loadings <- data.pca$rotation
summary(loadings)

##### There are two options to convert the loading to show the relative contribution, 
##### they both give the same answer so either can be used.
loadings_relative_A <- t(t(abs(loadings))/rowSums(t(abs(loadings))))*100
summary(loadings_relative_A)

loadings_relative_B <- sweep(x = abs(loadings), MARGIN = 2, STATS = colSums(abs(loadings)), FUN = "/")*100
summary(loadings_relative_B)

#### Plotting the PCA
##### First, I made a theme to change the background of the plot. Next, I changed the plot margins and the text size.
theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"), axis.text = element_text(size = 12), 
               legend.text = element_text(size = 12), legend.title = element_text(size = 12),
               text = element_text(size = 12))
##### Next, ggbiplot where obs.scale indicates the scale factor to apply to observation, 
##### var.scale indicates the scale factor to apply to variables, 
##### ellipse as TRUE draws a normal data ellipse for each group, 
##### and circle as TRUE draws a correlation circle.
g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = data.species, ellipse = TRUE, circle = TRUE)
g <- g + theme
g <- g + scale_colour_manual(values = c("#73dfff", "#ff8700", "#7a8ef5"))
g <- g + theme(legend.direction = 'horizontal', legend.position = 'bottom')
g

####################################################################################################################################################################################
# Multivariate environment similarity surface (MESS) analysis and mahalanobis distances
## Load layers 
### Training region layers
env.files.c <- list.files(path = "data/climate_processing/PresentLayers/", pattern = ".asc", full.names = TRUE) #Always check for tiff or asc
current <- stack(env.files.c)
names(current) <- c("alt", "bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")
current <- setMinMax(current)
projection(current) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
plot(current[[1]])

### Projection layers
#### These can be a different geographic region or different time period
env.files.p <- list.files(path = "data/climate_processing/Coastal_Plains_Layers/", pattern = ".asc", full.names = TRUE) #Always check for tiff or asc
proj <- stack(env.files.p)
names(proj) <- c("alt", "bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")
proj <- setMinMax(proj)
projection(proj) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
plot(proj[[1]])

## MESS
### For MESS, the more negative the score, the more likely you are to project 
### into novel environments. 

### Generate a data frame from the raster data
ref <- as.data.frame(current) 
head(ref)

### Use dismo to run the MESS. Negative results means dissimilarity
mess <- mess(x = proj, v = ref, full = FALSE)  
plot(mess)

## Mahalanobis distances.
### For mahalanobis distances, zero equals completely overlapping. 
### Comparatively, the higher the distance, 
### the more likely you are to project into an environment that falls 
### outside the range of the training region, leading to unreliable projections.

## Greater numbers means the variables are more dissimilar.
 
### Generate data frames containing the environmental data
refdat <- as.data.frame(ref) 
head(refdat)

projdat <- as.data.frame(proj) 
head(projdat)

### Calculate the average and covariance matrix of the variables 
### in the reference set
ref.av <- colMeans(refdat, na.rm=TRUE)
ref.cov <- var(refdat, na.rm=TRUE)

### Calculate the mahalanobis distance of each raster cell to the environmental center of the reference 
### set for both the reference and the projection data set and calculate the ratio between the two.
mah.ref <- mahalanobis(x=refdat, center=ref.av, cov=ref.cov, tol=1e-20)
mah.pro <- mahalanobis(x=projdat, center=ref.av, cov=ref.cov, tol=1e-20)
mah.max <- max(mah.ref[is.finite(mah.ref)])
nt2 <- as.data.frame(mah.pro / mah.max)

### Create and plot the raster layer
NT2 <- read.asciigrid(env.files.p[[1]])
NT2@data <- nt2
NT2rast <- raster(NT2)
plot(NT2rast, col=rainbow(100))
