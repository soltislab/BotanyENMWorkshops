# 03_Georeferencing.R
## Create a batch file for GeoLocate. 
### Created by ME Mabry 
### Modified by ML Gaynor.

### Load Packages
library(dplyr) 
library(gatoRs)

### Read in downloaded data frame
rawdf <- read.csv("data/download/raw/Shortia_galacifolia_raw_20230605.csv")

### Identify records for geoferencing
#### Filter for NA locality
rawdf_GeoRef <- need_to_georeference(rawdf)

### Format for GeoLocate
#### Subset needed columns 
rawdf_GeoRef <- rawdf_GeoRef %>%
                dplyr::select("locality string" = locality,
                              country = country,
                              state = stateProvince,
                              county = county,
                              latitude = latitude, 
                              longitude = longitude,
                              ID = ID, 
                              name = scientificName, 
                              basis = basisOfRecord)

#### Add required columns
rawdf_GeoRef$'correction status' <- ""
rawdf_GeoRef$precision <- ""
rawdf_GeoRef$'error polygon' <- ""
rawdf_GeoRef$'multiple results' <- ""

#### Reorder
rawdf_GeoRef2<- rawdf_GeoRef[,c("locality string", "country", "state", "county", "latitude", "longitude", "correction status", "precision", "error polygon", "multiple results", "ID", "name", "basis")]

### Save file for georeferencing
write.csv(rawdf_GeoRef2, file = "data/georeferencing/Shortia_galacifolia_Needing_GeoRef_20230605.csv",row.names = FALSE)


