# Occurrence_Data_Cleaning.R
## Occurrence Data Cleaning
## Modified based on script written by Charlotte Germain-Aubrey.
## Modified and created by ML Gaynor

### Load Packages
library(dplyr)
library(tidyr)
library(rjson)  
library(RCurl) 
library(raster)
library(sp)

### Load data files
# The datafiles for this demo are stored in "data/cleaning_demo/". 
raw.data <- read.csv("data/cleaning_demo/SampleFL-data.csv")
head(raw.data)


#### Exploring the datafile 
##### How many species are included in this file?
(species_count <- raw.data %>%
                  group_by(species) %>%
                  tally())

### Taxonomic Cleaning
#### Why? Lots of bad names due to misspelling, synonymns, and more  
#### How? **iPlant TNRS Tool**  

## In addition to accessing the iPlant TNRS Tool in R, you can also use the 
## website to access the database (http://tnrs.iplantcollaborative.org/). 
## In this tutorial, we only focus on checking synoymns using one tool, 
## however there are many other tools to check taxonomy.
## For example, taxize(https://ropensci.github.io/taxize-book/index.html) can 
## be used to accesss taxomny databases including Encylopedia of Life (eol), 
## Taxonomic Name Resolution Service (tnrs), Integrated Taxonomic Information Service (itis), 
## and [many more](https://ropensci.github.io/taxize-book/data-sources.html). 

#### Using the iPlant TNRS tool

##### Transform the species names list
###### First, we need to subset the names from the list we made above.

# Subset from the list we made above
(names <- species_count$species)

# Next, remove the underscores between the genus and species name, replace with spaces.
(names<-gsub('_',' ',names))

# Turn the list into a string by adding commas.
(names<-paste(names,collapse=','))

# Using the package **RCurl** makes the string URL-encoded.
names<-curlEscape(names)

##### Query the API
# Using the packages **rjson** and **RCurl** we can query the Application programming interface (API) for the iPlant TNRS tool.  
# We have to set an object equal to the API url address.
tnrs.api<-'http://tnrs.iplantc.org/tnrsm-svc'

# Next, we send a request to the address specified above. 
url<-paste(tnrs.api,'/matchNames?retrieve=best&names=',names,sep='')
tnrs.json<-getURL(url) 

##### Process the results
# The response is saved as tnrs.json above. Now, the response needs to be converted to a readable format from JSON. 
tnrs.results<-fromJSON(tnrs.json)

# Next, we need to summarize the names and change the format to match our raw data.
# First we isolate the names submitted and the found accepted name.
# Next, we replace the space in between the genus and species names with an underscore. 
# Then, we convert the data file to a dataframe. Finally, we rename the columns to match the names in raw.data.

corrected_names<-sapply(tnrs.results[[1]], function(x) c(x$nameSubmitted,x$acceptedName))
corrected_names <- gsub(" ", "_", corrected_names)
corrected_names<-as.data.frame(t(corrected_names),stringsAsFactors=FALSE)
names(corrected_names) <- c("species", "new")
corrected_names

##### Correcting names
# We now need to correct the names in the raw data file. To do this, we are going to merge the datasets, 
# so we can save the new name with all the information in our raw data file.
merged_datasets <- merge(raw.data, corrected_names, by = "species")
head(merged_datasets)

### Date Cleaning
# In addition to taxonmy fixes, we need to check the format of other columns, like dates.
# If you noticed above, when the year is missing there is an odd symbol. Lets replace that symbol with "NA"
merged_datasets$year <- gsub("\\N", NA, merged_datasets$year)

## Is there any other issues? 
# Visualize the year and data associated. 

(year_count <- merged_datasets %>%
               group_by(year) %>%
               tally())

### Location Cleaning
#### Precision
# How precise is the locality information? 

# Precision is influenced by the number of decimal places associated with latitude and longitude. 
# See below how degree and distance are related. 

location_table <- data.frame(places = c(0,1,2,3,4), 
degree = c(1, 0.1, 0.01, 0.001, 0.0001), 
distance = c("111 km", "11.1 km", "1.11 km", "111 m", "11.1 m"))

location_table 

# Here we are going to round to two decimal points. 

(merged_datasets$lat <- round(merged_datasets$lat, digits = 2))
merged_datasets$long <- round(merged_datasets$long, digits = 2)

#### Removing impossible points
# Prior to the next step we have 142 points, after we have 137.

merged_datasets <- merged_datasets %>%
                   filter(lat != 0, long != 0)


# Remove points that are botanical gardens.  
# First, load the data file of the botanical gardens in florida. 
bg.points <- read.csv("data/cleaning_demo/BotanicalGardensFloridaCoordinates.csv", header=TRUE)

# Next, we are going to filter the merged_dataset to exclude any coordinates that are in the bg.points dataset.
merged_datasets <- merged_datasets %>%
                   filter(!lat %in% bg.points$Lat & !long %in% bg.points$Long)

#### Removing duplicates
# We currently have 137 observations. Once only unique records are obtained, we have 121 observation. 

merged_datasets.unique <- merged_datasets %>%
                          distinct
head(merged_datasets.unique) 

#### Trimming to the desired area 
# Our desired area is Florida.
# Using the package **raster**  we are able to create a raster layer of the USA.
#If you want to download other countries using getData, check out the [avaliable maps](https://gadm.org/maps.html).
# Next, we can subset for Florida. State names are saved as "NAME_1" and county names are saved as "NAME_2" in the map file.   

# name = 'GADM' which is Database of Global Administrative Areas
# level = 2 gives the  level of subdivision (States)
usa <- getData('GADM', country='USA', level=2)

# subset
florida <- subset(usa, NAME_1=="Florida")


# Next, convert merged dataset unique into a spatial file. 

xy_data <- data.frame(x = merged_datasets.unique$long, y = merged_datasets.unique$lat)
coordinates(xy_data) <- ~ x + y
proj4string(xy_data) <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


# Extract the values of the points for the Florida polygon.
# Then, recombine everything to gather all data. 
# Next, filter for only the points in Florida, select the columns needed for MaxEnt input, and rename the new names as simply 'species'.
over <- over(xy_data, florida)
total <- cbind(merged_datasets.unique, over)
florida_points <- total %>%
                  filter(NAME_1 == "Florida") %>%
                  dplyr::select(new, lat, long) %>%
                  rename(species = new)
                  #rename(replace = c("species" = "new")) # AEM EDIT
head(florida_points)


### Finally, write a csv for the maxent steps.

# We need a csv with only species, lat, long for MaxEnt input
write.csv(florida_points, "data/cleaning_demo/MaxEntPointsInput_cleaned.csv", row.names = FALSE)

