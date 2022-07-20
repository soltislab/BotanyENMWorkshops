# Occurrence_Data_Cleaning.R
## Occurrence Data Cleaning
## Modified and created by ML Gaynor

### Load Packages
library(tidyverse)
library(raster)
library(sp)
library(sf)
library(spatstat)
library(spThin)
library(fields)
library(lubridate)
library(CoordinateCleaner)
library(ggspatial)
library(leaflet)

### Load functions 
#### This is a function I created with Natalie Patten.
source("functions/gators.R")

## Read in downloaded data frame
rawdf <- read.csv("data/download/raw/Shortia_galacifolia_raw_20220712.csv")

# Cleaning
## Inspect the data frame

### How many observations do we start with?
nrow(rawdf)

# 1. Resolve taxon names
## Inspect scientific names included in the raw df. 
unique(rawdf$name)

## Create a list of accepted names based on the name column in your data frame
search <-  c("Shortia galacifolia", "Sherwoodia galacifolia")


## Filter to only include accepted names:
df <-  filter_fix_names(rawdf, listofsynonyms = search, acceptedname = "Shortia galacifolia")

## How many observations do we have now?
nrow(df)

# 2. Decrease number of columns
## Merge the two locality columns
df$Latitude <- dplyr::coalesce(df$Latitude, df$spocc.latitude)
df$Longitude <- dplyr::coalesce(df$Longitude, df$spocc.longitude)

# Merge the two date columns
df$date <- dplyr::coalesce(df$date, df$spocc.date)

## Subset the columns
df <- df %>%
      dplyr::select(ID = ID, 
                    name = new_name, 
                    basis = basis, 
                    coordinateUncertaintyInMeters = coordinateUncertaintyInMeters, 
                    informationWithheld = informationWithheld, 
                    lat = Latitude, 
                    long = Longitude, 
                    date = date)

# 3. Clean localities 
## Filtering out NA's
df <- df %>%
      filter(!is.na(long)) %>%
      filter(!is.na(lat))

## How many observations do we have now?
nrow(df)

## Precision 
### Round to two decimal places
df$lat <- round(df$lat, digits = 2)
df$long <- round(df$long, digits = 2)

## Remove unlikely points
### Remove points at 0.00, 0.00
df <- df %>%
      filter(long != 0.00) %>%
      filter(lat != 0.00)

### Remove coordinates in cultivated zones, botanicals gardens, and outside our desired range
df <- CoordinateCleaner::cc_inst(df, 
              lon = "long", 
              lat = "lat", 
              species = "name")
### Next, we look for geographic outliers and remove outliers. 
df <- CoordinateCleaner::cc_outl(df, 
              lon = "long", 
              lat = "lat", 
              species = "name")

## How many observations do we have now?
nrow(df)

# 4. Remove Duplicates
## Fix dates
## Parse dates into the same format
df$date <- lubridate::ymd(df$date)

### Separate date into year, month, day
df <- df %>%
      mutate(year = lubridate::year(date), 
             month = lubridate::month(date), 
             day = lubridate::day(date))

## Remove rows with identical lat, long, year, month, and day
df <- distinct(df, lat, long, year, month, day, .keep_all = TRUE)

## How many observations do we have now?
nrow(df)

# 5. Spatial Correction
## Maxent will only retain one point per pixel. 
## To make the ecological niche analysis comparable, we will retain only one pt per pixel. 
## Read in raster file
bio1 <- raster("data/climate_processing/bioclim/bio_1.tif")
# Set resolution
rasterResolution <- max(res(bio1))

# Remove a point which nearest neighbor distance is smaller than the resolution size
## aka remove one point in a pair that occurs within one pixel
while(min(nndist(df[,6:7])) < rasterResolution){
  nnD <- nndist(df[,6:7])
  df <- df[-(which(min(nnD) == nnD) [1]), ]
}

## How many observations do we have now?
nrow(df)

## Spatial thinning
### Reduce the effects of sampling bias using randomization approach
## Calculate minimum nearest neighbor distance in km
nnDm <- rdist.earth(as.matrix(data.frame(lon = df$long, lat = df$lat)), miles = FALSE, R = NULL)
nnDmin <- do.call(rbind, lapply(1:5, function(i) sort(nnDm[,i])[2]))
min(nnDmin)

# Identify points to keep using spThin
keep <- spThin::thin(loc.data =  df, 
        verbose = FALSE, 
        long.col = "long", 
        lat.col = "lat",
        spec.col = "name",
        thin.par = 0.002, # Studies found 2m distance was enough to collect unique genets
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE)[[1]]

## Filter df to only include those lat/long
df <- df %>%
       filter((lat %in% keep$Latitude +
                long %in% keep$Longitude) == 2)
nrow(df)

# 6. Plot Cleaned Records
## Make points spatial 
df_fixed <- st_as_sf(df, coords = c("long", "lat"), crs = 4326)

## Set basemap
USA <- borders(database = "usa", colour = "gray50", fill = "gray50")
state <- borders(database = "state", colour = "black", fill = NA)

## Plot 
simple_map <- ggplot() +
              USA +
              state +
              geom_sf(df_fixed, 
                       mapping = aes(col = name), 
                       col = "blue") +
              coord_sf(xlim = c(min(df$long) - 3, max(df$long) + 3),
                       ylim = c(min(df$lat) - 3, max(df$lat) + 3)) +
              xlab("Longitude") +
              ylab("Latitude") +
              annotation_scale() +
              annotation_north_arrow(height = unit(1, "cm"), 
                                     width = unit(1, "cm"), 
                                     location = "tl")
simple_map

## Extra - Another fun way to view these points is with leaflet
leaflet(df_fixed) %>% 
  addMarkers(label = paste0(df$long, ", ", df$lat)) %>% 
  addTiles()

# 7. Save Cleaned.csv
write.csv(df, "data/cleaning_demo/Shortia_galacifolia_20220712-cleaned.csv", row.names = FALSE)

# 8. Make maxent ready
## Read in all cleaned files
alldf <- list.files("data/cleaning_demo/", full.names = TRUE, 
                    recursive = FALSE, include.dirs = FALSE, pattern = "*.csv")
alldf <- lapply(alldf, read.csv)
alldf <- do.call(rbind, alldf)


## Plot all records
### Make points spatial 
alldf_fixed <- st_as_sf(alldf, coords = c("long", "lat"), crs = 4326)

### Set basemap
USA <- borders(database = "usa", colour = "gray50", fill = "gray50")
state <- borders(database = "state", colour = "black", fill = NA)

### Plot 
all_map <- ggplot() +
            USA +
            state +
            geom_sf(alldf_fixed, 
                    mapping = aes(col = factor(name))) +
            coord_sf(xlim = c(min(alldf$long) - 3, max(alldf$long) + 3),
                     ylim = c(min(alldf$lat) - 3, max(alldf$lat) + 3)) +
            xlab("Longitude") +
            ylab("Latitude") +
            annotation_scale() +
            annotation_north_arrow(height = unit(1, "cm"), 
                                   width = unit(1, "cm"), 
                                   location = "tl")
all_map

## Select needed columns
alldf <- alldf %>%
         dplyr::select(name, lat, long)

## Save Maxent.csv
write.csv(alldf, "data/cleaning_demo/maxent_ready/diapensiaceae_maxentready_20220712.csv", row.names = FALSE)

    