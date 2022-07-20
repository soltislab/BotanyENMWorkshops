# 03_Georeferencing.R
## Create a batch file for GeoLocate. 
### Created by ME Mabry 

### Load Packages
library(dplyr) 
library(tidyr)
library(plyr) 
library(tibble) 

### Read in downloaded data frame
rawdf <- read.csv("data/download/raw/Shortia_galacifolia_raw_20220712.csv")

### Identify records for geoferencing
#### Merge the two locality columns
rawdf$Latitude <- dplyr::coalesce(rawdf$Latitude, rawdf$spocc.latitude)
rawdf$Longitude <- dplyr::coalesce(rawdf$Longitude, rawdf$spocc.longitude)

#### Subset for columns of interest
rawdf <- rawdf %>%
         dplyr::select(ID = ID, 
                      name = name, 
                      basis = basis, 
                      coordinateUncertaintyInMeters = coordinateUncertaintyInMeters, 
                      informationWithheld = informationWithheld,
                      country = country,
                      locality = locality,
                      lat = Latitude, 
                      long = Longitude, 
                      date = date,
                      state = state,
                      county = county)

#### Filter for NA locality
rawdf_GeoRef <- filter(rawdf, is.na(lat))

### Format for GeoLocate
#### Subset needed columns 
rawdf_GeoRef <- rawdf_GeoRef %>%
                dplyr::select("locality string" = locality,
                              country = country,
                              state = state,
                              county = county,
                              latitude = lat, 
                              longitude = long,
                              ID = ID, 
                              name = name, 
                              basis = basis)

#### Add required columns
rawdf_GeoRef$'correction status' <- ""
rawdf_GeoRef$precision <- ""
rawdf_GeoRef$'error polygon' <- ""
rawdf_GeoRef$'multiple results' <- ""

#### Reorder
rawdf_GeoRef2<- rawdf_GeoRef[,c("locality string", "country", "state", "county", "latitude", "longitude", "correction status", "precision", "error polygon", "multiple results", "ID", "name", "basis")]

### Save file for georeferencing
write.csv(rawdf_GeoRef2, file = "data/georeferencing/Shortia_galacifolia_Needing_GeoRef_20220712.csv",row.names = FALSE)


