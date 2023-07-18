# 01_Download_Occurrence_Data.R
## Download Occurrence Data
### ML Gaynor 

### Load packages 
library(ridigbio) 
library(gatoRs)
library(leaflet)


### Data download from iDigBio
#### Search for the species Galax urceolata
iDigBio_GU <- idig_search_records(rq=list(scientificname="Galax urceolata"))

#### Search for the family Diapensiaceae
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
write.csv(iDigBio_GU, "data/download/iDigBio_GU_20230605.csv", row.names = FALSE)
write.csv(iDigBio_GU_family, "data/download/iDigBio_GU_family_20230605.csv", row.names = FALSE)

### Data download using gatoRs
#### Make synonym lists
Shortia_galacifolia <- c("Shortia galacifolia", "Sherwoodia galacifolia")
Galax_urceolata <- c("Galax urceolata", "Galax aphylla")
Pyxidanthera_barbulata <- c("Pyxidanthera barbulata","Pyxidanthera barbulata var. barbulata")
Pyxidanthera_brevifolia <- c("Pyxidanthera brevifolia", "Pyxidanthera barbulata var. brevifolia")

#### Use the gators_download function
gators_download(synonyms.list = Shortia_galacifolia,
                write.file = TRUE, 
                filename = "data/download/raw/Shortia_galacifolia_raw_20230605.csv")
gators_download(synonyms.list = Galax_urceolata, 
              write.file = TRUE, 
              filename = "data/download/raw/Galax_urceolata_raw_20230605.csv")
gators_download(synonyms.list = Pyxidanthera_barbulata, 
              write.file = TRUE, 
              filename = "data/download/raw/Pyxidanthera_barbulata_raw_20230605.csv")
gators_download(synonyms.list = Pyxidanthera_brevifolia, 
              write.file = TRUE, 
              filename = "data/download/raw/Pyxidanthera_brevifolia_raw_20230605.csv")

### Quick-look at the downloaded files
rawdf <- read.csv("data/download/raw/Shortia_galacifolia_raw_20230605.csv")

#### Inspect the data frame
#### What columns are included?
names(rawdf)

### How many observations do we start with?
nrow(rawdf)

### Where are these points?
#### The error message here indicates many points do not have long/lat values (more in 02)
leaflet(rawdf) %>% 
  addMarkers(label = paste0(rawdf$longitude, ", ", rawdf$latitude)) %>% 
  addTiles()

