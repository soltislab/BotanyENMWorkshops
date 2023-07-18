# Occurrence_Data_Cleaning.R
## Occurrence Data Cleaning
## Modified and created by ML Gaynor

### Load Packages
library(gatoRs)
library(fields)
library(sf)
library(ggplot2)
library(ggspatial)
library(leaflet)


## Read in downloaded data frame
rawdf <- read.csv("data/download/raw/Shortia_galacifolia_raw_20230605.csv")

# Cleaning
## Inspect the data frame
### How many observations do we start with?
nrow(rawdf)

# 1. Resolve taxon names
## Here we are going to harmonize taxonomy using taxa_clean()
## This function has three filter options: exact, fuzzy, or interactive.
## Inspect scientific names included in the raw df. 
unique(rawdf$scientificName)

## Create a list of accepted names based on the name column in your data frame
search <-  c("Shortia galacifolia", "Sherwoodia galacifolia")

## Filter to only include accepted names:
df <-  taxa_clean(df = rawdf, 
                  synonyms.list = search,
                  taxa.filter = "fuzzy",
                  accepted.name = "Shortia galacifolia")

## How many observations do we have now?
nrow(df)

# 2. Clean localities
## Here we remove any records with missing coordinates, 
## impossible coordinates, coordinates at (0,0), 
## and any that are flagged as skewed.
## We also round the provided latitude and longitude values
## to a specified number of decimal places.   
df <- basic_locality_clean(df = df,  
                           remove.zero = TRUE, # Records at (0,0) are removed
                           precision = TRUE, # latitude and longitude are rounded 
                           digits = 2, # round to 2 decimal places
                           remove.skewed = TRUE)

## How many observations do we have now?
nrow(df)

### Remove coordinates in cultivated zones, 
### botanicals gardens, and outside our desired range
df <- process_flagged(df, interactive = FALSE)

## How many observations do we have now?
nrow(df)

# 3. Remove Duplicates
## Here we identify and remove both (1) specimen duplicates and 
## (2) aggregator duplicates based on each specimens coordinates,
## occurrenceID, and eventDate. To leverage all date information 
## available, set `remove.unparseable = FALSE` to manually populate 
## the year, month, and day columns.Here, we also confirm all ID 
## (UUID and key) are unique to remove any  within-aggregator 
## duplicates that may accumulate due to processing errors.
df <- remove_duplicates(df, remove.unparseable = TRUE)

## How many observations do we have now?
nrow(df)

## Bonus question: How many more records do you retain when `remove.unparseable = FALSE`?


# 4. Spatial Correction
## One point per pixel
## Maxent will only retain one point per pixel. 
## To make the ecological niche analysis comparable, 
## we will retain only one pt per pixel. 
### Note: Default is a raster with 30 arc sec resolution.
df <- one_point_per_pixel(df)

## How many observations do we have now?
nrow(df)

## Spatial thinning
### Reduce the effects of sampling bias using randomization approach

### Step 1: What should your minimum distance be?
## Calculate minimum nearest neighbor distance in km
nnDm <- rdist.earth(as.matrix(data.frame(lon = df$longitude, lat = df$latitude)), miles = FALSE, R = NULL)
nnDmin <- do.call(rbind, lapply(1:5, function(i) sort(nnDm[,i])[2]))
min(nnDmin)

### Here the current minimum distance is 2.22 km. 
### Based on literature, we find a 2 meters (or 0.002 km)
### distance was enough to collect unique genets, 
### so we do not need to thin our points. 

###  Step 2: Thin occurrence records using spThin through gatoRs.   
## When you do need to thin your records, 
## here is a great function to do so!
df <- thin_points(df, 
                    distance = 0.002, # in km 
                    reps = 100)


# 5. Plot Cleaned Records
## Make points spatial 
df_fixed <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)

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
              coord_sf(xlim = c(min(df$longitude) - 3, max(df$longitude) + 3),
                       ylim = c(min(df$latitude) - 3, max(df$latitude) + 3)) +
              xlab("Longitude") +
              ylab("Latitude") +
              annotation_scale(plot_unit = "km") +
              annotation_north_arrow(height = unit(1, "cm"), 
                                     width = unit(1, "cm"), 
                                     location = "tl")
simple_map

## Extra - Another fun way to view these points is with leaflet
leaflet(df_fixed) %>% 
  addMarkers(label = paste0(df$longitude, ", ", df$latitude)) %>% 
  addTiles()

# 6. Save Cleaned.csv
write.csv(df, "data/cleaning_demo/Shortia_galacifolia_20230605-cleaned.csv", 
          row.names = FALSE)

# 7. Repeat for all species
## Set up for loop
files <- list.files("data/download/raw", full.names = TRUE)[1:3]
synonymns <- list(Galax_urceolata = c("Galax urceolata",  "Galax urceolata (Poir.) Brummitt",  "Galax urceolata (Poiret) Brummitt", "Galax urceolaa", "Galax aphylla L.", "Galax aphylla"),  
               Pyxidanthera_barbulata =c("Pyxidanthera barbulata", "Pyxidanthera barbulata Michx.", "Pyxidanthera barbulata var. barbulata" , "Pyxidenthera barbulata"), 
               Pyxidanthera_brevifolia = c("Pyxidanthera brevifolia", "Pyxidanthera brevifolia Wells", "Pyxidanthera barbulata var. brevifolia (Wells) H.E.Ahles", "Pyxidanthera barbulata var. brevifolia"  ))

## Repeat cleaning steps for the remaining taxa
for(i in 1:3){
    df <- read.csv(files[i])
    # Taxa clean
    search <-  synonymns[[i]]
    df <-  taxa_clean(df = df, synonyms.list = search, 
                      taxa.filter = "exact", 
                      accepted.name = search[1])
    # Locality clean 
    df <- basic_locality_clean(df = df,  remove.zero = TRUE, 
                               precision = TRUE, digits = 2, 
                               remove.skewed = TRUE)
    df <- process_flagged(df, interactive = FALSE)
    # Remove duplicates
    df <- remove_duplicates(df, remove.unparseable = TRUE)
    # Spatial correct
    df <- one_point_per_pixel(df)
    df <- thin_points(df,  distance = 0.002, reps = 100)
    # Save file
    outfile <- paste0("data/cleaning_demo/",
                      gsub(" ", "_", search[1]), 
                      "_20230605-cleaned.csv")
    write.csv(df, outfile, 
              row.names = FALSE)
    rm(df, search, outfile)
}

# 8. Make maxent ready
## Read in all cleaned files
alldf <- list.files("data/cleaning_demo/", full.names = TRUE, 
                    recursive = FALSE, include.dirs = FALSE, 
                    pattern = "*.csv")
alldf <- lapply(alldf, read.csv)
alldf <- do.call(rbind, alldf)


## Plot all records
### Make points spatial 
alldf_fixed <- st_as_sf(alldf, coords = c("longitude", "latitude"), 
                        crs = 4326)

### Set basemap
USA <- borders(database = "usa", colour = "gray50", fill = "gray50")
state <- borders(database = "state", colour = "black", fill = NA)

### Plot 
all_map <- ggplot() +
            USA +
            state +
            geom_sf(alldf_fixed, 
                    mapping = aes(col = factor(accepted_name))) +
            coord_sf(xlim = c(min(alldf$longitude) - 3, max(alldf$longitude) + 3),
                     ylim = c(min(alldf$latitude) - 3, max(alldf$latitude) + 3)) +
            xlab("Longitude") +
            ylab("Latitude") +
            labs(color = "Scientific name") +
            annotation_scale(plot_unit = "km") +
            annotation_north_arrow(height = unit(1, "cm"), 
                                   width = unit(1, "cm"), 
                                   location = "tl")
all_map

## Select needed columns
alldfchomp <- c()
species <- unique(alldf$accepted_name)
for(i in 1:4){
  alldfchomp[[i]] <- data_chomp(alldf[alldf$accepted_name ==species[i], ], 
                                accepted.name =  species[i])
}
alldf <- do.call(rbind, alldfchomp)

## Save Maxent.csv
write.csv(alldf, "data/cleaning_demo/maxent_ready/diapensiaceae_maxentready_20230605.csv", row.names = FALSE)

    