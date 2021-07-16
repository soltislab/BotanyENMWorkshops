# Climate Processing
## Climate Layer Processing 
## Modified and created by ML Gaynor. 
## Based on code by Mike Belitz (mbelitz/Odo_SDM_Rproj)

# Load Packages 
library(raster)
library(gtools)
library(dplyr)
library(rgdal)
library(sp)
library(rangeBuilder)
library(sf)
library(caret)
library(usdm)
library(dismo)
library(stringr)

# Load functions
source("functions/VIFLayerSelect.R")

# Load bioclim layers
biolist <- list.files("data/climate_processing/bioclim/", pattern = "*.tif", full.names = TRUE)

## Order list using gtools
biolist <- mixedsort(sort(biolist))

### Load rasters
biostack <- raster::stack(biolist)

# Load occurrence records
alldf <- read.csv("data/cleaning_demo/maxent_ready/diapensiaceae_maxentready_20210625.csv")
alldf$name <- as.character(alldf$name)

# Present Layers - all
## Here we will make the projection layers, or the layers which include shared space. 
## First we have to define the accessible space.
alldfsp <- alldf

## Make into a spatial point data frame
coordinates(alldfsp) <- ~ long + lat
proj4string(alldfsp) <- CRS("+proj=longlat +datum=WGS84")

## Create alpha hull
hull <- getDynamicAlphaHull(x = alldfsp@coords, 
                             fraction = 1, # min. fraction of records we want included
                             partCount = 1, # number of polygons
                             initialAlpha = 20, # initial alpha size, 20m
                             clipToCoast = "terrestrial",
                             proj = "+proj=longlat +datum=WGS84")
### Visualize
plot(hull[[1]], col=transparentColor('gray50', 0.5), border = NA)
points(x = alldf$long, y = alldf$lat, cex = 0.5, pch = 3)

## Add buffer to hull
### Transform into CRS related to meters
hullTrans <- spTransform(hull[[1]], "+proj=cea +lat_ts=0 +lon_0=0")
alldfspTrans <- spTransform(alldfsp, "+proj=cea +lat_ts=0 +lon_0")

### Calculate buffer size
#### Here we take the 80th quantile of the max distance between points
buffDist <- quantile(x = (apply(spDists(alldfspTrans), 2, FUN = function(x) sort(x)[2])), 
                     probs = 0.80, na.rm = TRUE) 

### Buffer the hull
buffer_m <- buffer(x = hullTrans, width = buffDist, dissolve = TRUE)
buffer <- spTransform(buffer_m, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

### Visualize
plot(buffer, col=transparentColor('gray50', 0.5), border = NA)
points(x = alldf$long, y = alldf$lat, cex = 0.5, pch = 3)

## Mask and crop bioclim layers
path <- "data/climate_processing/all/"
end <- ".asc"
for(i in 1:length(biolist)){
    # Subset raster layer
    rast <- biostack[[i]]
    # Setup file names
    name <- names(rast)
    out <- paste0(path, name)
    outfile <- paste0(out, end)
    # Crop and mask
    c <- crop(rast, extent(buffer))
    c <- mask(c, buffer)
    # Write raster
    writeRaster(c, outfile, format = "ascii", overwrite = TRUE)
}


# Select layers for MaxEnt
## We only want to include layers that are not highly correlated.
## To assess which layers we will include, we want to look at the pearson correlation coefficient among layers.
### Stack all layers
clippedlist <- list.files("data/climate_processing/all/", pattern = "*.asc", full.names = TRUE)
clippedlist <- mixedsort(sort(clippedlist))
clippedstack <- raster::stack(clippedlist)

### Then calculate the correlation coefficient
corr <- layerStats(clippedstack, 'pearson', na.rm=TRUE)

### Isolate only the pearson correlation coefficient and take absolute value
c <- abs(corr$`pearson correlation coefficient`)

### Write file and view in excel
write.csv(c, "data/climate_processing/correlationBioclim.csv", row.names = FALSE)
## Highly correlated layers (> |0.80|) can impact the statistical significance 
## of the niche models and therefore must be removed. 

# Randomly select variables to remove
envtCor <- mixedsort(sort(findCorrelation(c, cutoff = 0.80, names = TRUE, exact = TRUE)))
envtCor

## Variable inflation factor (VIF)
### VIF can detect for multicollinearity in a set of multiple regression variables. 
### Run a simple maxent model for every species and calculate the average permutation contribution
#### Loop through each species and save permutation importance in list
set.seed(195)
m <- c()
for(i in  1:length(unique(alldf$name))){
    species <- unique(alldf$name)[i]
    spp_df <-  alldf %>%
               dplyr::filter(name == species)
    coordinates(spp_df) <- ~ long + lat
    model <- maxent(x = clippedstack, p = coordinates(spp_df), progress = "text", silent = FALSE) 
    m[[i]] <- vimportance(model)
}

#### Bind the dataframes
mc <- do.call(rbind, m)

#### Calculate the mean and rename columns
mc_average <- aggregate(mc[, 2], list(mc$Variables), mean)
mc_average <- mc_average %>%
              dplyr::select(Variables = Group.1, permutation.importance = x)
mc1 <- mc_average

# Use VIF and the MaxEnt permutation importance to select the best variables for your model.
## Note, this leads to different layers when the models are rerun 
## without setting seed due to permutations having randomness to them. 
selectedlayers <- VIF_layerselect(clippedstack, mc_average)
mixedsort(sort(names(selectedlayers)))

## Since this can vary per system (despite setting seed), we added this line to keep our files consistent for the workshop
sl <- c("bio_3", "bio_7", "bio_8", "bio_9", "bio_14", "bio_15", "bio_18", "elev")
selectedlayers <- raster::subset(clippedstack, sl)

## Copy selected layers to Present Layer folder
for(i in 1:length(names(selectedlayers))){
    name <- names(selectedlayers)[i]
    from <- paste0("data/climate_processing/all/", name, ".asc")
    to <- paste0("data/climate_processing/PresentLayers/all/", name, ".asc")
    file.copy(from, to,
              overwrite = TRUE, recursive = FALSE, 
              copy.mode = TRUE)
}

# Create Species Training Layers

for(i in 1:length(unique(alldf$name))){
    species <- unique(alldf$name)[i]
    # Subset species from data frame
    spp_df <-  alldf %>%
               dplyr::filter(name == species)
    # Make spatial
    coordinates(spp_df) <- ~ long + lat
    proj4string(spp_df) <- CRS("+proj=longlat +datum=WGS84")
    
    ## Create alpha hull
    sphull <- getDynamicAlphaHull(x = spp_df@coords, 
                                fraction = 1, # min. fraction of records we want included
                                partCount = 1, # number of polygons
                                initialAlpha = 20, # initial alpha size, 20m
                                clipToCoast = "terrestrial",
                                proj = "+proj=longlat +datum=WGS84")
    
    ### Transform into CRS related to meters
    sphullTrans <- spTransform(sphull[[1]], "+proj=cea +lat_ts=0 +lon_0=0")
    spp_dfTrans <- spTransform(spp_df, "+proj=cea +lat_ts=0 +lon_0")
    
    ### Calculate buffer size
    #### Here we take the 80th quantile of the max distance between points
    spbuffDist <- quantile(x = (apply(spDists(spp_dfTrans), 2, FUN = function(x) sort(x)[2])), 
                         probs = 0.80, na.rm = TRUE) 
    
    ### Buffer the hull
    spbuffer_m <- buffer(x = sphullTrans, width = spbuffDist, dissolve = TRUE)
    spbuffer <- spTransform(spbuffer_m, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    ### Crop and Mask
    spec <- gsub(" ", "_", species)
    path <- paste0("data/climate_processing/PresentLayers/", spec,"/")
    end <- ".asc"
    for(j in 1:length(names(selectedlayers))){
        # Subset raster layer
        rast <- selectedlayers[[j]]
        # Setup file names
        name <- names(rast)
        out <- paste0(path, name)
        outfile <- paste0(out, end)
        # Crop and mask
        c <- crop(rast, extent(spbuffer))
        c <- mask(c, spbuffer)
        # Write raster
        writeRaster(c, outfile, format = "ascii", overwrite = TRUE)
    }
}

