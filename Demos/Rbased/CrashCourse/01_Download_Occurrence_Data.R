# Download_Occurrence_Data.R
## Download Occurrence Data
### Modified and created by ML Gaynor 

### Load packages 
library(dplyr) 
library(tidyr) 
library(plyr) 
library(spocc) 
library(ridigbio) 
library(tibble) 
library(rbison)

### Load functions 
#### This is a function I created with Natalie Patten.
#### It will be part of her R package gatoRs (Geographic And Taxonomic Occurrence R-based Scrubbing).
source("functions/gators.R")

### Data download from iDigBio
#### Search for the species Galax urceolata
iDigBio_GU <- idig_search_records(rq=list(scientificname="Galax urceolata"))

#### Search for the family ADiapensiaceae
iDigBio_GU_family <- idig_search_records(rq=list(family="Diapensiaceae"), limit=1000)


### What if you want to read in all the points for a family within an extent?
#### Hint: Use the [iDigBio portal](https://www.idigbio.org/portal/search) to
#### determine the bounding box for your region of interest. 
rq_input <- list("scientificname"=list("type"="exists"),
                 "family"="Diapensiaceae", 
                 geopoint=list(
                   type="geo_bounding_box",
                   top_left=list(lon = -98.16, lat = 48.92),
                   bottom_right=list(lon = -64.02, lat = 23.06)
                   )
                 )

##### Search using the input you just made
iDigBio_GU_family_USA <- idig_search_records(rq_input, limit=1000)

#### Save as csv files 
write.csv(iDigBio_GU, "data/download/iDigBio_GU_20210614.csv", row.names = FALSE)
write.csv(iDigBio_GU_family, "data/download/iDigBio_GU_family_20210614.csv", row.names = FALSE)

### Data download using spocc_combined
#### Make synonym lists
Shortia_galacifolia <- c("Shortia galacifolia", "Sherwoodia galacifolia")
Galax_urceolata <- c("Galax urceolata", "Galax aphylla")
Pyxidanthera_barbulata <- c("Pyxidanthera barbulata","Pyxidanthera barbulata var. barbulata")
Pyxidanthera_brevifolia <- c("Pyxidanthera brevifolia", "Pyxidanthera barbulata var. brevifolia")

#### Use the spocc_combine function
spocc_combine(Shortia_galacifolia, "data/download/raw/Shortia_galacifolia_raw_20210614.csv")
spocc_combine(Galax_urceolata, "data/download/raw/Galax_urceolata_raw_20210614.csv")
spocc_combine(Pyxidanthera_barbulata, "data/download/raw/Pyxidanthera_barbulata_raw_20210614.csv")
spocc_combine(Pyxidanthera_brevifolia, "data/download/raw/Pyxidanthera_brevifolia_raw_20210614.csv")

