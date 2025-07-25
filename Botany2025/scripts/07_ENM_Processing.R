# ENM_Processing.R
# ---------------------------------------------------------------
# Purpose:
#   - Project ecological niche models (ENMs) for four species onto
#     current and future climate layers (Eastern Temperate Forests)
#   - Visualize habitat suitability maps
#   - Calculate niche breadth for each species
#   - Calculate pairwise niche overlap
#   - Test for phylogenetic signal in niche overlap
#   - Compare current vs future projections to assess range shifts
# ---------------------------------------------------------------

# Increase Java memory limit to avoid OutOfMemoryError
# Some ENM packages, especially dismo and older ENMeval workflows,
# may call Java routines that require more RAM.
# This sets the Java heap space to 8GB.
options(java.parameters = "-Xmx30g")

## ---- Load Required Packages ----

library(terra)            # For reading and manipulating raster data (modern replacement for raster package)
library(dismo)            # For predicting from ENMs and compatibility with Maxent models
library(ENMeval)          # For calculating niche overlap statistics
library(ENMTools)         # For calculating niche breadth (Levins' B2)
library(ggplot2)          # For creating publication-quality plots
library(rnaturalearth)    # For downloading shapefiles of countries and states
library(viridis)          # For colorblind-friendly color scales
library(ggspatial)        # For adding north arrows to plots
library(dplyr)            # For data manipulation
library(ape)              # For reading and processing phylogenetic trees
library(biomod2)          # For range size change calculations between binary rasters

## ---- A) Project Models onto Current Climate (Eastern Temperate Forests) ----

# List all climate raster files (.asc format) for the Eastern Temperate Forests (ETF) region.
# These will be used as environmental predictors for projecting the ENMs.
clim_files <- list.files("data/04_climate_processing/CurrentEasternTemperateForests",pattern = "\\.asc$", full.names = TRUE)

# Read all climate rasters into a single SpatRaster object for efficient processing.
clim_stack <- terra::rast(clim_files)

# Plot a quick preview of the first four layers to check that rasters loaded correctly.
plot(clim_stack[[1:4]])

# Define a vector of the species you want to process.
# Each species has an ENM already built and saved to disk.
species_list <- c("Galax_urceolata",
                  "Pyxidanthera_barbulata",
                  "Pyxidanthera_brevifolia",
                  "Shortia_galacifolia")

# Load the USA state boundaries as an sf object.
# This will be used as a basemap for plotting projections.
usa <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf")

# Loop through each species to project their ENM onto the current climate layers.
for (species in species_list) {

  cat("\nProjecting current climate for:", species, "\n")

  # Load the optimal ENM for this species.
  # The object `opt.mod` will be loaded into the workspace.
  mod_file <- paste0("data/05_ENMs/", species, "_opt_mod.RData")
  load(mod_file)  # loads object: opt.mod

  # Identify which layers were used to build the ENM.
  # This ensures that your projection uses the same predictor variables
  # as the model was trained on.
  vif_dir <- paste0("data/04_climate_processing/Cropped/", species, "/VIF/")
  species_layers <- list.files(vif_dir, full.names = TRUE)

  # Read in those specific rasters
  spec_stack <- terra::rast(species_layers)
  layer_names <- names(spec_stack)

  # Subset the full ETF raster stack to only the layers used in the ENM.
  ETF_Rasters <- terra::subset(clim_stack, layer_names)
  #if error, try: ETF_Rasters <- raster::stack(ETF_Rasters)

  # Project the ENM onto the current ETF climate layers.
  proj_file <- paste0("data/06_ENM_processing/", species, "_EFT_Projection.asc")

  p <- dismo::predict(opt.mod, ETF_Rasters, filename = proj_file, overwrite = TRUE, NAflag = -9999)

  # Convert the raster to a data frame for plotting in ggplot.
  p_df <- as.data.frame(p, xy = TRUE)
  colnames(p_df)[3] <- "Habitat_Suitability"

  # Get the bounding box coordinates of the raster for plotting.
  bbox_vals <- terra::ext(p)

  # Create the habitat suitability map.
  p_plot <- ggplot() +
    geom_sf(data = usa, fill = "grey70") +  # basemap
    geom_tile(data = p_df, aes(x = x, y = y, fill = Habitat_Suitability)) +
    geom_sf(data = usa, fill = NA, color = "white", linewidth = 0.5) +
    coord_sf(xlim = c(bbox_vals$xmin - 2, bbox_vals$xmax + 2),
             ylim = c(bbox_vals$ymin - 2, bbox_vals$ymax + 2),
             expand = FALSE) +
    scale_fill_viridis_c(name = "Suitability") +
    labs(title = paste0(gsub("_", " ", species), " - Eastern Temperate Forest Habitat Suitability"),
         x = "Longitude",
         y = "Latitude") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 10)) +
    annotation_north_arrow(location = "tl",
                           height = unit(1, "cm"),
                           width = unit(1, "cm"))

  # Save the plot as a PNG file.
  png_file <- paste0("data/06_ENM_processing/", species, "_ETF_Projection_Plot.png")

  ggsave(filename = png_file, plot = p_plot, width = 7, height = 5, dpi = 300)
}

## ---- B) Calculate Niche Breadth ----

# For each species, calculate niche breadth across the ETF region.
# Niche breadth (Levins' B2) indicates whether a species is a generalist (high value) or specialist (low value).
for (species in species_list) {

  cat("\nCalculating niche breadth for:", species, "\n")

  proj_file <- paste0("data/06_ENM_processing/", species, "_EFT_Projection.asc")

  proj_raster <- terra::rast(proj_file)

  # Compute Levins' B2 niche breadth
  breadth <- ENMTools::raster.breadth(proj_raster)

  cat("Niche breadth:", round(breadth$B2, 3), "\n")
}

## ---- C) Calculate Niche Overlap ----

# Stack all species projections into a single SpatRaster object.
# This allows calculation of pairwise niche overlap (Schoener's D).
enm_stack <- c(terra::rast("data/06_ENM_processing/Galax_urceolata_EFT_Projection.asc"),
               terra::rast("data/06_ENM_processing/Pyxidanthera_barbulata_EFT_Projection.asc"),
               terra::rast("data/06_ENM_processing/Pyxidanthera_brevifolia_EFT_Projection.asc"),
               terra::rast("data/06_ENM_processing/Shortia_galacifolia_EFT_Projection.asc")
)

# Assign species names as layer names in the raster stack.
names(enm_stack) <- c("Galax urceolata",
                      "Pyxidanthera barbulata",
                      "Pyxidanthera brevifolia",
                      "Shortia galacifolia")
enm_stack_raster <- stack(enm_stack)

# Calculate pairwise Schoener's D overlap between species.
# Values range from 0 (no overlap) to 1 (identical niches).
overlap_matrix <- calc.niche.overlap(enm_stack_raster, overlapStat = "D")

print(overlap_matrix)

## ---- D) Test Phylogenetic Signal in Niche Overlap ----

# Read the phylogenetic tree of the species.
tree <- ape::read.tree("data/06_ENM_processing/diapensiaceae_subset.tre")

# Drop the outgroup so only focal species remain.
tree <- ape::drop.tip(tree, "Cyrilla_racemiflora")

# Replace underscores with spaces to match overlap matrix species names.
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Calculate pairwise phylogenetic distances between species.
phylo_dist <- cophenetic(tree)

# Ensure the phylogenetic distance matrix matches the species order in the overlap matrix.
species_names <- rownames(overlap_matrix)
species_names <- gsub("\\.", " ", species_names)
phylo_dist <- phylo_dist[species_names, species_names]

# Initialize vectors to store pairwise values.
pair_names <- NULL
phylo_vals <- NULL
overlap_vals <- NULL

# Loop through the lower triangle of the matrices to extract pairwise comparisons.
for (i in 2:length(species_names)) {
  for (j in 1:(i - 1)) {
    pair_names <- c(pair_names,
                    paste(species_names[i], species_names[j], sep = " - "))
    phylo_vals <- c(phylo_vals, phylo_dist[i, j])
    overlap_vals <- c(overlap_vals, overlap_matrix[i, j])
  }
}

# Fit a linear regression between phylogenetic distance and niche overlap.
lm_res <- lm(overlap_vals ~ phylo_vals)
s <- summary(lm_res)
print(s)

# Format regression statistics for display on the plot.
label_text <- paste0("Intercept = ", format(round(coef(lm_res)[1], 3)), "\n",
                     "Slope = ", format(round(coef(lm_res)[2], 5)), "\n",
                     "P-value = ", format(signif(s$coefficients[2, 4], 3)), "\n",
                     "R² = ", format(round(s$r.squared, 3)))

# Plot the relationship between phylogenetic distance and niche overlap.
plot(phylo_vals, overlap_vals,
     xlab = "Phylogenetic Distance",
     ylab = "Niche Overlap (Schoener's D)",
     pch = 19)
abline(lm_res, col = "red", lwd = 2)
text(phylo_vals - 3,
     overlap_vals - 0.02,
     labels = 1:length(pair_names),
     pos = 3,
     cex = 0.8)
text(x = max(phylo_vals) * 0.7,
     y = max(overlap_vals) * 0.9,
     labels = label_text,
     col = "red",
     cex = 0.9, pos = 4)
legend("topleft",
       legend = paste(1:length(pair_names), pair_names, sep = ": "),
       bty = "n", cex = 0.7)

# For niche Age Overlap Correlations (AOCs), a high intercept and negative slope
# indicate niche conservatism,meaning closely related species tend to occupy
# similar ecological niches. A low intercept combined with a positive slope
# suggests niche divergence, where closely related species occupy different
# niches.The results here are pretty meaningless since we're looking at very few,
# distantly related species, but it serves the purpose of the demo. A linear
# regression indicated no significant relationship between phylogenetic distance
# and niche overlap among the four species,suggesting minimal phylogenetic signal
# in ecological niche similarity.

## ---- E) Project Models onto Future Climate ----

# List future climate layers under the SSP370 scenario for 2081-2100.
future_files <- list.files("data/04_climate_processing/ACCESS-CM2_2081-2100_ssp370/", pattern = "\\.asc$", full.names = TRUE)

# Add the current elevation layer to the future stack, since elevation remains constant.
elev_file <- list.files("data/04_climate_processing/CurrentEasternTemperateForests/", pattern = "elev\\.asc$", full.names = TRUE)

futurestack <- terra::rast(c(future_files, elev_file))

# Project each species onto future ETF climate conditions.
for (species in species_list) {

  cat("\nProjecting future climate for:", species, "\n")

  load(paste0("data/05_ENMs/", species, "_opt_mod.RData"))

  vif_dir <- paste0("data/04_climate_processing/Cropped/", species, "/VIF/")

  species_layers <- list.files(vif_dir, full.names = TRUE)
  spec_stack <- terra::rast(species_layers)
  layer_names <- names(spec_stack)

  FutureETF_Rasters <- terra::subset(futurestack, layer_names)
  #if error, try: FutureETF_Rasters <- raster::stack(FutureETF_Rasters)

  future_file <- paste0("data/06_ENM_processing/", species, "_Future_ETF_Projection.asc")

  p_future <- dismo::predict(opt.mod, FutureETF_Rasters, filename = future_file, overwrite = TRUE, NAflag = -9999)

  save(p_future, file = sub("\\.asc$", ".RData", future_file))

  # Plot the future habitat suitability map.
  p_future_df <- as.data.frame(p_future, xy = TRUE)
  colnames(p_future_df)[3] <- "Habitat_Suitability"

  bbox_vals <- terra::ext(p_future)

  p_plot <- ggplot() +
    geom_sf(data = usa, fill = "grey70") +
    geom_tile(data = p_future_df,
              aes(x = x, y = y, fill = Habitat_Suitability)) +
    geom_sf(data = usa, fill = NA, color = "white", linewidth = 0.5) +
    coord_sf(xlim = c(bbox_vals$xmin - 2, bbox_vals$xmax + 2),
             ylim = c(bbox_vals$ymin - 2, bbox_vals$ymax + 2),
             expand = FALSE) +
    scale_fill_viridis_c(name = "Future Suitability") +
    labs(title = paste0(gsub("_", " ", species), " - Future Habitat Suitability"),
         x = "Longitude",
         y = "Latitude") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 10)) +
    annotation_north_arrow(location = "tl",
                           height = unit(1, "cm"),
                           width = unit(1, "cm"))

  png_file <- paste0("data/06_ENM_processing/", species, "_Future_ETF_Projection.png")

  ggsave(filename = png_file, plot = p_plot, width = 7, height = 5, dpi = 300)
}

## ---- F) Compare Current and Future Projections ----

# Initialize a list to hold change maps for each species.
all_dfs <- list()

for (species in species_list) {

  cat("\nComparing projections for:", species, "\n")

  # Load current and future suitability rasters.
  current <- terra::rast(paste0("data/06_ENM_processing/", species, "_EFT_Projection.asc"))
  future <- terra::rast(paste0("data/06_ENM_processing/", species, "_Future_ETF_Projection.asc"))

  # Threshold each raster at 0.7 to produce binary presence/absence maps.
  current_binary <- current >= 0.7
  future_binary <- future >= 0.7

  # Compute range changes between current and future projections.
  RangeSizeDiff <- BIOMOD_RangeSize(as(current_binary, "Raster"),
                                    as(future_binary, "Raster"))

  diff_raster <- RangeSizeDiff$Diff.By.Pixel

  df <- as.data.frame(diff_raster, xy = TRUE)
  colnames(df) <- c("x", "y", "change")

  # Convert numeric change codes to categorical labels for plotting.
  df$fill <- factor(df$change,
                    levels = c(-2, -1, 0, 1),
                    labels = c("lost", "unchanged", "not occupied", "occupied in future"))

  # Add species name for faceting plots.
  df$species <- gsub("_", " ", species)

  all_dfs[[species]] <- df
}

# Combine all species into one data frame for faceted plotting.
combined_df <- bind_rows(all_dfs)

# Determine overall plotting extent.
combined_bbox <- combined_df %>%
  summarise(xmin = min(x, na.rm = TRUE),
            xmax = max(x, na.rm = TRUE),
            ymin = min(y, na.rm = TRUE),
            ymax = max(y, na.rm = TRUE))

# Plot all species' range shifts together.
p_all <- ggplot(combined_df) +
  geom_sf(data = usa, fill = "grey70") +
  geom_tile(aes(x = x, y = y, fill = fill)) +
  geom_sf(data = usa, fill = NA, color = "white", linewidth = 0.5) +
  coord_sf(xlim = c(combined_bbox$xmin - 2, combined_bbox$xmax + 2),
           ylim = c(combined_bbox$ymin - 2, combined_bbox$ymax + 2),
           expand = FALSE) +
  scale_fill_viridis_d(name = "Pixel Change", na.value = "white") +
  facet_wrap(~ species) +
  labs(title = "Current vs Future Habitat Suitability Change",
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 10),
        strip.text = element_text(size = 9)) +
  annotation_north_arrow(
    location = "tl",
    height = unit(1, "cm"),
    width = unit(1, "cm"))

ggsave("data/06_ENM_processing/all_species_range_change_map.png",
       plot = p_all,
       width = 10,
       height = 8,
       dpi = 300)
