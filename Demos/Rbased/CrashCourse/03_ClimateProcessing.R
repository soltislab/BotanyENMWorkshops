# ClimateProcesssing.R
## Climate Layer Processing
## Original script by Charlotte Germain-Aubrey.
## Modified and created by ML Gaynor. 

# Load Packages 
library(maptools)
library(raster)
library(rgdal)
library(sp)
library(maps)
library(mapproj)

# Load bioclim layers
alt_l <- raster("data/climate_processing/Bioclim/alt.bil")
bio1_l <- raster("data/climate_processing/Bioclim/bio1.bil")
bio2_l <- raster("data/climate_processing/Bioclim/bio2.bil")
bio3_l <- raster("data/climate_processing/Bioclim/bio3.bil")
bio4_l <- raster("data/climate_processing/Bioclim/bio4.bil")
bio5_l <- raster("data/climate_processing/Bioclim/bio5.bil")
bio6_l <- raster("data/climate_processing/Bioclim/bio6.bil")
bio7_l <- raster("data/climate_processing/Bioclim/bio7.bil")
bio8_l <- raster("data/climate_processing/Bioclim/bio8.bil")
bio9_l <- raster("data/climate_processing/Bioclim/bio9.bil")
bio10_l <- raster("data/climate_processing/Bioclim/bio10.bil")
bio11_l <- raster("data/climate_processing/Bioclim/bio11.bil")
bio12_l <- raster("data/climate_processing/Bioclim/bio12.bil")
bio13_l <- raster("data/climate_processing/Bioclim/bio13.bil")
bio14_l <- raster("data/climate_processing/Bioclim/bio14.bil")
bio15_l <- raster("data/climate_processing/Bioclim/bio15.bil")
bio16_l <- raster("data/climate_processing/Bioclim/bio16.bil")
bio17_l <- raster("data/climate_processing/Bioclim/bio17.bil")
bio18_l <- raster("data/climate_processing/Bioclim/bio18.bil")
bio19_l <- raster("data/climate_processing/Bioclim/bio19.bil")

# Define desired extent
## We are going to use Florida as our desired extent
FL <- rgdal::readOGR("data/climate_processing/FL/FLstate2.shp")

# Crop bioclim layers to the desired extent
# Visualize the first layer
plot(alt_l)
## First, mask the bioclim layer with Florida.
alt <- mask(alt_l, FL)
## Visualize the masked layers 
plot(alt)
## Next, crop the bioclim layers to Florida's extent.
alt <- crop(alt, extent(FL))
## Visualize the final layer
plot(alt)
## Now save the layer - note we already saved these layers for you, so you have to overwrite the file to save yours
writeRaster(alt, "data/climate_processing/PresentLayers/alt.asc", format="ascii", overwrite = TRUE)

# Repeat with all additional files 
bio1 <- mask(bio1_l, FL)
bio1 <- crop(bio1, extent(FL))
#writeRaster(bio1, "data/climate_processing/PresentLayers/bio1.asc", format="ascii")

bio2 <- mask(bio2_l, FL)
bio2 <- crop(bio2, extent(FL))
#writeRaster(bio2, "data/climate_processing/PresentLayers/bio2.asc", format="ascii")

bio3 <- mask(bio3_l, FL)
bio3 <- crop(bio3, extent(FL))
#writeRaster(bio3, "data/climate_processing/PresentLayers/bio3.asc", format="ascii")

bio4 <- mask(bio4_l, FL)
bio4 <- crop(bio4, extent(FL))
#writeRaster(bio4, "data/climate_processing/PresentLayers/bio4.asc", format="ascii")

bio5 <- mask(bio5_l, FL)
bio5 <- crop(bio5, extent(FL))
#writeRaster(bio5, "data/climate_processing/PresentLayers/bio5.asc", format="ascii")

bio6 <- mask(bio6_l, FL)
bio6 <- crop(bio6, extent(FL))
#writeRaster(bio6, "data/climate_processing/PresentLayers/bio6.asc", format="ascii")

bio7 <- mask(bio7_l, FL)
bio7 <- crop(bio7, extent(FL))
#writeRaster(bio7, "data/climate_processing/PresentLayers/bio7.asc", format="ascii")

bio8 <- mask(bio8_l, FL)
bio8 <- crop(bio8, extent(FL))
#writeRaster(bio8, "data/climate_processing/PresentLayers/bio8.asc", format="ascii")

bio9 <- mask(bio9_l, FL)
bio9 <- crop(bio9, extent(FL))
#writeRaster(bio9, "data/climate_processing/PresentLayers/bio9.asc", format="ascii")

bio10 <- mask(bio10_l, FL)
bio10 <- crop(bio10, extent(FL))
#writeRaster(bio10, "data/climate_processing/PresentLayers/bio10.asc", format="ascii")

bio11 <- mask(bio11_l, FL)
bio11 <- crop(bio11, extent(FL))
#writeRaster(bio11, "data/climate_processing/PresentLayers/bio11.asc", format="ascii")

bio12 <- mask(bio12_l, FL)
bio12 <- crop(bio12, extent(FL))
#writeRaster(bio12, "data/climate_processing/PresentLayers/bio12.asc", format="ascii")

bio13 <- mask(bio13_l, FL)
bio13 <- crop(bio13, extent(FL))
#writeRaster(bio13, "data/climate_processing/PresentLayers/bio13.asc", format="ascii")

bio14 <- mask(bio14_l, FL)
bio14 <- crop(bio14, extent(FL))
#writeRaster(bio14, "data/climate_processing/PresentLayers/bio14.asc", format="ascii")

bio15 <- mask(bio15_l, FL)
bio15 <- crop(bio15, extent(FL))
#writeRaster(bio15, "data/climate_processing/PresentLayers/bio15.asc", format="ascii")

bio16 <- mask(bio16_l, FL)
bio16 <- crop(bio16, extent(FL))
#writeRaster(bio16, "data/climate_processing/PresentLayers/bio16.asc", format="ascii")

bio17 <- mask(bio17_l, FL)
bio17 <- crop(bio17, extent(FL))
#writeRaster(bio17, "data/climate_processing/PresentLayers/bio17.asc", format="ascii")

bio18 <- mask(bio18_l, FL)
bio18 <- crop(bio18, extent(FL))
#writeRaster(bio18, "data/climate_processing/PresentLayers/bio18.asc", format="ascii")

bio19 <- mask(bio19_l, FL)
bio19 <- crop(bio19, extent(FL))
#writeRaster(bio19, "data/climate_processing/PresentLayers/bio19.asc", format="ascii")

# Select layers for MaxEnt
## We only want to include layers that are not highly correlated.
## To assess which layers we will include, we want to look at the pearson correlation coefficient among layers.
### Stack all layers
stack <- stack(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)
### Then calculate the correlation coefficient
corr <- layerStats(stack, 'pearson', na.rm=TRUE)
### Isolate only the pearson correlation coefficient 
c <- corr$`pearson correlation coefficient`

### Write file and view in excel
#write.csv(c, "data/climate_processing/correlationBioclim.csv")

## Highly correlated layers (> |0.80|) can impact the statistical significance 
## of the niche models and therefore must be removed. 



