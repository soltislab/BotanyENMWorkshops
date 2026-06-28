# Occurrence_Data_Cleaning_Updated.R
# ---------------------------------------------------------------
# Purpose: Clean and harmonize occurrence records for ecological niche modeling.
# Created and modified by ML Gaynor
# ---------------------------------------------------------------

## ---- Load Required Packages ----
library(gatoRs)         # Taxonomic and geographic cleaning
library(fields)         # For calculating spatial distances
library(sf)             # For spatial operations and visualization
library(ggplot2)        # For plotting
library(ggspatial)      # Scale bars and north arrows
library(leaflet)        # For interactive mapping

## ---- A) Read Raw Occurrence Data ----

# Read species-specific raw data
rawdf <- read.csv("data/01_download/raw/Shortia_galacifolia_raw_2025_06_27.csv")
nrow(rawdf)  # Starting number of records

## ---- B) Taxonomic Harmonization ----

# Inspect unique names
unique(rawdf$scientificName)

# Define synonym list
search <- c("Shortia galacifolia", "Sherwoodia galacifolia")

# Clean names with fuzzy matching
df <- taxa_clean(df = rawdf,
                 synonyms.list = search,
                 taxa.filter = "fuzzy",
                 accepted.name = "Shortia galacifolia")
nrow(df)  # Updated count

## ---- C) Locality Cleaning ----

# Remove invalid/skewed records and round coordinates
df <- basic_locality_clean(df = df,
                           remove.zero = TRUE,
                           precision = TRUE,
                           digits = 2,
                           remove.skewed = TRUE)
nrow(df)

## ---- D) Flag and Filter Cultivated Coordinates ----

df <- process_flagged(df, interactive = FALSE, scientific.name = "accepted_name")
nrow(df)

## ---- E) Duplicate Removal ----

# Remove duplicate specimens and aggregator artifacts
df <- remove_duplicates(df, remove.unparseable = TRUE)
nrow(df)

## ---- F) Spatial Deduplication ----

# Retain one point per raster pixel (30 arc-sec default)
df <- one_point_per_pixel(df)
nrow(df)

## ---- G) Spatial Thinning ----

# Calculate minimum nearest-neighbor distance
nnDm <- rdist.earth(as.matrix(df[, c("longitude", "latitude")]), miles = FALSE)
nnDmin <- do.call(rbind, lapply(1:5, function(i) sort(nnDm[, i])[2]))
min(nnDmin)  # e.g., 2.22 km

# If thinning is needed
df <- thin_points(df, distance = 0.002, reps = 100)
nrow(df)

## ---- H) Static Map of Cleaned Points ----

df_fixed <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
USA <- borders(database = "usa", colour = "gray80", fill = "gray80")
state <- borders(database = "state", colour = "black", fill = NA)

simple_map <- ggplot() +
  USA + state +
  geom_sf(data = df_fixed, color = "blue") +
  coord_sf(xlim = c(min(df$longitude) - 3, max(df$longitude) + 3),
           ylim = c(min(df$latitude) - 3, max(df$latitude) + 3)) +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale() +
  annotation_north_arrow(location = "tl", height = unit(1, "cm"), width = unit(1, "cm"))

simple_map

## ---- I) Interactive Map with Leaflet ----

leaflet(df_fixed) %>%
  addMarkers(label = paste0(df$longitude, ", ", df$latitude)) %>%
  addTiles()

## ---- J) Save Cleaned CSV ----

write.csv(df, "data/02_cleaning/Shortia_galacifolia_2025_06_27_cleaned.csv", row.names = FALSE)

## ---- K) Batch Clean for All Species ----

files <- list.files("data/01_download/raw", full.names = TRUE)[1:3]
synonymns <- list(
  Galax_urceolata = c("Galax urceolata", "Galax urceolata (Poir.) Brummitt", "Galax urceolata (Poiret) Brummitt", "Galax urceolaa", "Galax aphylla L.", "Galax aphylla"),
  Pyxidanthera_barbulata = c("Pyxidanthera barbulata", "Pyxidanthera barbulata Michx.", "Pyxidanthera barbulata var. barbulata", "Pyxidenthera barbulata"),
  Pyxidanthera_brevifolia = c("Pyxidanthera brevifolia", "Pyxidanthera brevifolia Wells", "Pyxidanthera barbulata var. brevifolia (Wells) H.E.Ahles", "Pyxidanthera barbulata var. brevifolia")
)

for (i in 1:3) {
  df <- read.csv(files[i])
  search <- synonymns[[i]]

  df <- taxa_clean(df, synonyms.list = search, taxa.filter = "fuzzy", accepted.name = search[1])
  df <- basic_locality_clean(df, remove.zero = TRUE, precision = TRUE, digits = 2, remove.skewed = TRUE)
  df <- process_flagged(df, interactive = FALSE, scientific.name = "accepted_name") # removing all points for brevifolia
  df <- remove_duplicates(df, remove.unparseable = TRUE)
  df <- one_point_per_pixel(df)
  df <- thin_points(df, distance = 0.002, reps = 100)

  outfile <- paste0("data/02_cleaning/", gsub(" ", "_", search[1]), "_2025_06_27_cleaned.csv")
  write.csv(df, outfile, row.names = FALSE)
  rm(df, search, outfile)
}

## ---- L) Merge Cleaned Files for Maxent ----

alldf <- list.files("data/02_cleaning", pattern = "*.csv", full.names = TRUE)
alldf <- lapply(alldf, read.csv)
alldf <- do.call(rbind, alldf)

write.csv(alldf, "data/02_cleaning/maxent_ready/diapensiaceae_maxentready_2025_06_27.csv", row.names = FALSE)

## ---- M) Map All Records ----

alldf_fixed <- st_as_sf(alldf, coords = c("longitude", "latitude"), crs = 4326)

all_map <- ggplot() +
  USA + state +
  geom_sf(data = alldf_fixed, aes(col = factor(accepted_name))) +
  coord_sf(xlim = c(min(alldf$longitude) - 3, max(alldf$longitude) + 3),
           ylim = c(min(alldf$latitude) - 3, max(alldf$latitude) + 3)) +
  xlab("Longitude") + ylab("Latitude") +
  labs(color = "Scientific name") +
  annotation_scale() +
  annotation_north_arrow(location = "tl", height = unit(1, "cm"), width = unit(1, "cm"))

all_map

## ---- N) Prepare GeoLocate Batch File ----

rawdf <- read.csv("data/01_download/raw/Shortia_galacifolia_raw_2025_06_27.csv")
rawdf_GeoRef <- need_to_georeference(rawdf)

# Format columns for GeoLocate submission
rawdf_GeoRef <- rawdf_GeoRef %>%
  dplyr::select(
    "locality string" = locality,
    country,
    state = stateProvince,
    county,
    latitude,
    longitude,
    ID,
    name = scientificName,
    basis = basisOfRecord
  )

# Add required empty columns
rawdf_GeoRef$'correction status' <- ""
rawdf_GeoRef$precision <- ""
rawdf_GeoRef$'error polygon' <- ""
rawdf_GeoRef$'multiple results' <- ""

# Reorder columns
rawdf_GeoRef2 <- rawdf_GeoRef[, c("locality string", "country", "state", "county", "latitude", "longitude",
                                  "correction status", "precision", "error polygon", "multiple results",
                                  "ID", "name", "basis")]

write.csv(rawdf_GeoRef2, "data/03_georeferencing/Shortia_galacifolia_Needing_GeoRef_2025_06_27.csv", row.names = FALSE)
