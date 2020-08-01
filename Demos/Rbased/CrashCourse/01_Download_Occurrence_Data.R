# Download_Occurrence_Data.R
## Download Occurrence Data
## Original script by AE Melton. 
### Modified and created by ML Gaynor 

### Load packages 
library(ridigbio)
library(spocc)
library(scrubr)


### Data download from iDigBio
#### Search for the species Liatris cylindracea
iDigBio_LC <- idig_search_records(rq=list(scientificname="Liatris cylindracea"))

#### Search for the family Asteraceae
iDigBio_LC_family <- idig_search_records(rq=list(family="Asteraceae"), limit=1000)


### What if you want to read in all the points for a family within an extent?
#### Hint: Use the [iDigBio portal](https://www.idigbio.org/portal/search) to
#### determine the bounding box for your region of interest. 

rq_input <- list("scientificname"=list("type"="exists"),
                 "family"="asteraceae", 
                 geopoint=list(
                   type="geo_bounding_box",
                   top_left=list(lon=-87.86, lat=30.56),
                   bottom_right=list(lon=-79.21, lat= 24.78)
                 )
)
##### Search using the input you just made
iDigBio_LC_family_florida <- idig_search_records(rq_input, limit=1000)

#### Save as csv files 
write.csv(iDigBio_LC, "data/download/iDigBio_LC_072819.csv", row.names = FALSE)
write.csv(iDigBio_LC_family, "data/download/iDigBio_LC_family_072819.csv", row.names = FALSE)


### Data download using spocc
#### Using the package *spocc* we download data from iDigBio and GBIF.
#### We specified that coordinates are needed. Notice that the occurrence records for iDigBio and GBIF are stored sperately. 
(spocc_LC <- occ(query = "Liatris cylindracea", from = c('gbif','idigbio'), has_coords = T))


#### Processing spocc download
##### Synoymns
###### Sometimes the 'name' in the database does not match the 'scientific name'- this can be due to synoymns or variety names. 
cbind(spocc_LC$gbif$data$Liatris_cylindracea$acceptedScientificName[340], spocc_LC$gbif$data$Liatris_cylindracea$name[340])

###### Since this could be an issue when processing the data, we are going to fix the name using the package *scurbr*. 
spocc_LC_namefixed <- fixnames(spocc_LC, how = "query")

##### Data format
###### Since the occurrence records for iDigBio and GBIF are stored sperately,
###### next we have to convert the occurrence download into a single data.frame (df) where all query are combined.

spocc_LC_df <- occ2df(spocc_LC_namefixed)



#### Save as csv files 
write.csv(spocc_LC_df, "data/download/spocc_LC_df_072819.csv", row.names = FALSE)


