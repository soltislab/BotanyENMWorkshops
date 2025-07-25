# Download_Occurrence_Data_Updated.R
# ---------------------------------------------------------------
# Purpose: Download and preview species occurrence records using iDigBio and gatoRs.
# Created by ML Gaynor
# ---------------------------------------------------------------

## ---- Load Required Packages ----
library(ridigbio)     # Interface to iDigBio API
library(gatoRs)       # Unified taxonomic/occurrence tools
library(leaflet)      # Interactive mapping

## ---- A) Download from iDigBio ----

# Search for specific species (Galax urceolata)
iDigBio_GU <- idig_search_records(rq = list(scientificname = "Galax urceolata"))

# Search for all Diapensiaceae occurrences (limit to 1000 for example)
iDigBio_GU_family <- idig_search_records(rq = list(family = "Diapensiaceae"), limit = 1000)

# Search with a geographic bounding box (e.g., eastern USA extent)
rq_input <- list(
  scientificname = list(type = "exists"),
  family = "Diapensiaceae",
  geopoint = list(
    type = "geo_bounding_box",
    top_left = list(lon = -98.16, lat = 48.92),
    bottom_right = list(lon = -64.02, lat = 23.06)
  )
)
iDigBio_GU_family_USA <- idig_search_records(rq_input, limit = 1000)

# Save iDigBio results as CSV
write.csv(iDigBio_GU, "data/01_download/iDigBio_GU_2025_06_27.csv", row.names = FALSE)
write.csv(iDigBio_GU_family, "data/01_download/iDigBio_GU_family_2025_06_27.csv", row.names = FALSE)

## ---- B) Download Using gatoRs ----

# Define synonym lists for each focal species
Shortia_galacifolia <- c("Shortia galacifolia", "Sherwoodia galacifolia")
Galax_urceolata <- c("Galax urceolata", "Galax aphylla")
Pyxidanthera_barbulata <- c("Pyxidanthera barbulata", "Pyxidanthera barbulata var. barbulata")
Pyxidanthera_brevifolia <- c("Pyxidanthera brevifolia", "Pyxidanthera barbulata var. brevifolia")

# Use gatoRs to download records with resolved synonyms
gators_download(synonyms.list = Shortia_galacifolia,
                write.file = TRUE,
                filename = "data/01_download/raw/Shortia_galacifolia_raw_2025_06_27.csv")

gators_download(synonyms.list = Galax_urceolata,
                write.file = TRUE,
                filename = "data/01_download/raw/Galax_urceolata_raw_2025_06_27.csv")

gators_download(synonyms.list = Pyxidanthera_barbulata,
                write.file = TRUE,
                filename = "data/01_download/raw/Pyxidanthera_barbulata_raw_2025_06_27.csv")

gators_download(synonyms.list = Pyxidanthera_brevifolia,
                write.file = TRUE,
                filename = "data/01_download/raw/Pyxidanthera_brevifolia_raw_2025_06_27.csv")

## ---- C) Preview Downloaded Files ----

# Read one downloaded file
rawdf <- read.csv("data/01_download/raw/Shortia_galacifolia_raw_2025_06_27.csv")

# Preview columns and dimensions
names(rawdf)
nrow(rawdf)

# Visualize records interactively - The error message here indicates many points do not have long/lat values 
leaflet(rawdf) %>% 
  addTiles() %>% 
  addMarkers(label = paste0(rawdf$longitude, ", ", rawdf$latitude))
