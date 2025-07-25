# ClimateProcessing_Updated.R
# ---------------------------------------------------------------
# Purpose: Process climatic raster data and occurrence records for Diapensiaceae SDMs.
# Includes alpha hull generation, raster cropping, VIF selection, and correlation analysis.
# ---------------------------------------------------------------

## ---- Load Required Packages ----
library(terra)        # For raster handling
library(sf)           # For vector spatial operations
library(ggplot2)      # For plotting
library(dplyr)        # For data manipulation
library(usdm)         # For VIF calculations
library(stringr)      # For string operations
library(gridExtra)    # For plotting

## ---- A) Load Cleaned Occurrence Data ----

# Read in cleaned occurrence data with geographic coordinates and species ID
alldf <- read.csv("data/02_cleaning/maxent_ready/diapensiaceae_maxentready_2025_06_27.csv")

# Convert occurrence data to an sf object
alldfsp <- st_as_sf(alldf,
                    coords = c("longitude", "latitude"),
                    crs = 4326)

## ---- B) Load and Stack Climatic Variables ----

# List and sort all .tif files in the BioClim directory
biolist <- list.files("data/04_climate_processing/BioClim/",
                      pattern = "\\.tif$",
                      full.names = TRUE)

# Stack the raster files into a SpatRaster
biostack <- terra::rast(biolist)

## ---- C) Crop and Mask Bioclim Layers to Accessible Area (M) for Each Species ----

# Define base output directory for cropped rasters
dir  <- "data/04_climate_processing/Cropped"

# Loop through each unique species in the occurrence data
for (species in unique(alldf$accepted_name)) {
  
  # Print current species for progress tracking
  message("Processing species: ", species)
  
  # Subset occurrence data for this species
  species_df <- alldf %>%
    filter(accepted_name == species)
  
  # Convert species occurrences to sf
  species_sf <- st_as_sf(species_df,
                         coords = c("longitude", "latitude"),
                         crs = 4326)
  
  # Run alpha hull via rangeBuilder
  hull <- rangeBuilder::getDynamicAlphaHull(
    x = species_df,
    coordHeaders = c("longitude", "latitude"),
    fraction = 1,
    partCount = 1,
    initialAlpha = 20,
    clipToCoast = "terrestrial",
    verbose = FALSE
  )
  
  # Convert hull to sf
  hull_sf <- st_as_sf(hull[[1]])
  
  # Transform hull and points to equal-area projection
  proj <- "+proj=cea +lat_ts=0 +lon_0=0"
  hull_sf_proj <- st_transform(hull_sf, crs = proj)
  species_sf_proj <- st_transform(species_sf, crs = proj)
  
  # Calculate nearest-neighbor distances
  dist_matrix <- st_distance(species_sf_proj)
  nearest_neighbor_dists <- apply(dist_matrix, 2, function(x) sort(x)[2])
  buffer_distance <- quantile(nearest_neighbor_dists, probs = 0.80, na.rm = TRUE)
  
  # Buffer the hull
  buffer_proj <- st_buffer(hull_sf_proj, dist = buffer_distance)
  
  # Union polygons in case of multiple pieces
  buffer_proj <- st_union(buffer_proj)
  
  # Transform both hull and buffered hull back to WGS84 for plotting
  hull_sf_wgs84 <- st_transform(hull_sf_proj, crs = 4326)
  buffer_wgs84 <- st_transform(buffer_proj, crs = 4326)
  
  # plot hull before buffering
  p_before <- ggplot() +
    geom_sf(data = hull_sf_wgs84, fill = "gray70", alpha = 0.5, color = NA) +
    geom_sf(data = species_sf, color = "black", size = 0.5) +
    theme_minimal() +
    ggtitle(paste("Alpha Hull (No Buffer) -", species))
  
  # plot after buffering alpha hull
  p_after <- ggplot() +
    geom_sf(data = buffer_wgs84, fill = "gray70", alpha = 0.5, color = NA) +
    geom_sf(data = species_sf, color = "red", size = 0.5) +
    theme_minimal() +
    ggtitle(paste("Buffered Hull -", species))
  
  gridExtra::grid.arrange(p_before, p_after, ncol = 2)
  
  # Convert buffered area to terra SpatVector
  buffer_vect <- terra::vect(buffer_wgs84)
  
  # Define output directory for this species
  spec <- gsub(" ", "_", species)
  species_dir <- file.path(dir, spec)
  dir.create(species_dir, showWarnings = FALSE)
  
  # Loop through each raster layer and crop/mask to buffered area
  for (i in seq_along(biolist)) {
    
    rast <- biostack[[i]]
    
    # Crop raster to extent
    rast_crop <- crop(rast, buffer_vect)
    
    # Mask raster to buffered hull shape
    rast_mask <- mask(rast_crop, buffer_vect)
    
    # Define output file name
    layer_name <- names(rast)
    outfile <- file.path(species_dir, paste0(layer_name, ".tif"))
    
    # Write cropped raster
    writeRaster(rast_mask,outfile,overwrite = TRUE)
    
    message("  Cropped and saved: ", outfile)
  }
}

## ---- D) VIF Selection for Each Species ----

# Loop through each species
for (species in unique(alldf$accepted_name)) {
  spec <- gsub(" ", "_", species)
  spec_dir <- file.path(dir, spec)
  
  # Read in species cropped rasters and stack them
  cropped_files <- list.files(spec_dir,
                              pattern = "\\.tif$",
                              full.names = TRUE)
  
  species_stack <- terra::rast(cropped_files)
  
  # Calculate VIF and exclude collinear variables
  vif_result <- usdm::vifstep(species_stack, th = 10)
  message("VIF calculated for species: ", species)
  
  reduced_stack <- exclude(species_stack, vif_result)
  
  # Create a folder for VIF-selected layers
  vif_dir <- file.path(spec_dir, "VIF")
  dir.create(vif_dir, showWarnings = FALSE)
  
  # Save each retained raster layer
  for (layer_name in names(reduced_stack)) {
    
    r_layer <- reduced_stack[[layer_name]]
    
    out_file <- file.path(vif_dir, paste0(layer_name, ".tif"))
    
    terra::writeRaster(r_layer,
                       out_file,
                       overwrite = TRUE)
    
    message("Saved layer: ", out_file)
  }
}

## ---- E) Correlation Matrix of Cropped Bioclim Layers for Each Species ----

# List species folders in the cropped directory
species_dirs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)

# Loop through each species folder
for (species_dir in species_dirs) {
  
  # Extract species name from folder path
  species_name <- basename(species_dir)
  
  # List all cropped raster files for this species
  cropped_files <- list.files(species_dir,
                              pattern = "\\.tif$",
                              full.names = TRUE)
  
  # Stack the rasters for this species
  cropped_stack <- terra::rast(cropped_files)
    
  # Calculate correlation matrix
  corr <- terra::layerCor(cropped_stack, fun = "cor", na.rm = TRUE)
    
  # Extract absolute Pearson correlation values
  cor_abs <- abs(corr$correlation)
    
  # Define output CSV path
  corr_outfile <- file.path(species_dir,paste0("correlationBioclim_", species_name, ".csv"))
    
  # Write correlation matrix to CSV
  write.csv(cor_abs,corr_outfile,row.names = TRUE)

}
