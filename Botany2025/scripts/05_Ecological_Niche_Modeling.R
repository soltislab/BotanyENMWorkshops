# Ecological_Niche_Modeling.R
# ---------------------------------------------------------------
# Purpose: Build ecological niche models for Diapensiaceae species
# using ENMeval and Maxent with VIF-selected environmental variables.
# ---------------------------------------------------------------

## ---- Load Required Packages ----
library(terra)         # For raster data
library(ENMeval)       # For ENM evaluation and tuning
library(predicts)      # For niche modeling tools
library(dplyr)         # For data manipulation
library(ggplot2)       # Optional for additional plots
library(sf)            # For vector spatial operations
library(rnaturalearth) # For downloading shapefiles of countries and states
library(ggspatial)     # Scale bars and north arrows

## ---- A) Load Cleaned Occurrence Data ----

# Load occurrence data
alldf <- read.csv("data/02_cleaning/maxent_ready/diapensiaceae_maxentready_2025_06_27.csv")

## ---- B) Subset for Single Species ----

# Example species: Galax urceolata
Galax_urceolata <- alldf %>%
  filter(accepted_name == "Galax urceolata")

## ---- C) Load VIF-Selected Climate Layers ----

# List environmental layers
vif_list <- list.files(
  "data/04_climate_processing/Cropped/Galax_urceolata/VIF/",
  full.names = TRUE,
  recursive = FALSE
)

# Load as raster stack
vifStack <- terra::rast(vif_list)

## ---- D) Run ENMeval to Generate Models ----

# ENMeval model generation:
# - Tests combinations of feature classes (FC)
# - Tests different regularization multipliers (RM)

eval <- ENMevaluate(
  occs = Galax_urceolata[, c("longitude", "latitude")],
  envs = vifStack,
  tune.args = list(fc = c("L", "Q"), rm = 1:2),
  partitions = "block",
  n.bg = 10000,
  parallel = FALSE,
  algorithm = 'maxent.jar',
)

# Save ENMeval object
save(eval, file = "data/05_ENMs/Galax_urceolata_ENM_eval.RData")

## ---- E) Visualize Model Predictions ----

# Extract predictions from ENMeval object
maps <- eval@predictions

# Plot first 6 prediction rasters
terra::plot(maps,
            nc = 2,
            main = names(maps))

## ---- F) Calculate Niche Overlap Between Models ----

# Calculate Schoener's D among all models
mod_overlap <- calc.niche.overlap(maps, overlapStat = "D")
print(mod_overlap)

## ---- G) Examine Overall Tuning Results ----

# Retrieve summary table
results <- eval.results(eval)
head(results)

## Plot model tuning results
evalplot.stats(
  e = eval,
  stats = c("or.10p", "auc.val"),
  color = "fc",
  x.var = "rm")

## ---- H) Select Optimal Model Based on Criteria ----
opt.seq <- results %>%
  filter(!is.na(AICc)) %>%                      # Exclude models with NA AICc
  filter(AICc == min(AICc)) %>%                 # Minimum AICc
  filter(or.10p.avg != 0) %>%                   # Exclude zero omission
  filter(or.10p.avg == min(or.10p.avg)) %>%     # Minimum omission
  filter(auc.val.avg == max(auc.val.avg))       # Maximum AUC

opt.seq

# Save optimal model summary
write.table(
  opt.seq,
  file = "data/05_ENMs/Galax_urceolata_OptModel.txt",
  sep = "\t",
  row.names = FALSE
)

## ---- I) Visualize the variable contributions  ----

# Retrieve the optimal model from ENMeval object
opt.mod <- eval.models(eval)[[opt.seq$tune.args]]

# Save optimal model as RData
save(opt.mod, file = "data/05_ENMs/Galax_urceolata_opt_mod.RData")

# Plot variables and their contribution to the model
png("data/05_ENMs/Galax_urceolata_Variable_Contribution.png",
    width = 1200,
    height = 800,
    res = 150)

# Plot variable contributions
dismo::plot(opt.mod, main = "Variable Contribution - Optimal Model")

dev.off()

## ---- J) Visualize Response Curves ----

# Plot response curves for optimal model
predicts::partialResponse(eval@models[[opt.seq$tune.args]], var = "elev")

## ---- K) Plot Optimal Model ----
opt.pred <- eval.predictions(eval)[[as.character(opt.seq$tune.args)]]

r_df <- as.data.frame(opt.pred, xy = TRUE)
colnames(r_df) <- c("x", "y", "suitability")

#  Occurrences  sf
occ_sf <- st_as_sf(eval@occs, coords = c("longitude", "latitude"), crs = 4326)

# get bounding box for plotting from occurences
bbox_vals <- sf::st_bbox(occ_sf)

# USA states
usa <- ne_states(country = "United States of America", returnclass = "sf")

# Plot
p <- ggplot() +
  geom_tile(data = r_df, aes(x = x, y = y, fill = suitability)) +
  scale_fill_viridis_c(name = "Suitability") +
  geom_sf(data = occ_sf, color = "red", size = 0.1) +
  geom_sf(data = usa, fill = NA, color = "white", linewidth = 0.5) +
  coord_sf(
    xlim = c(bbox_vals["xmin"] - 2, bbox_vals["xmax"] + 2),
    ylim = c(bbox_vals["ymin"] - 2, bbox_vals["ymax"] + 2),
    expand = FALSE) +
  labs(
    title = "Predicted Suitability (Optimal ENMeval Model)",
    x = "Longitude",
    y = "Latitude")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey80", color = NA),
    plot.background = element_rect(fill = "white", color = NA)) +
  annotation_north_arrow(location = "tl", height = unit(1, "cm"), width = unit(1, "cm"))

print(p)

# Save Outputs
ggsave(filename = "data/05_ENMs/Galax_urceolata_Optimal_Model_Map.png",
       plot = p,
       width = 8,
       height = 5,
       dpi = 300)

writeRaster(x = opt.pred,
            filename = "data/05_ENMs/Galax_urceolata_ENM_optModel.asc",
            overwrite = TRUE,
            NAflag = -9999)


## ---- L) Loop Through Other Species ----

# List all unique species
species_list <- unique(alldf$accepted_name)

# Exclude Galax urceolata (already processed)
species_list <- setdiff(species_list, "Galax urceolata")

# Loop over remaining species
for (sp in species_list) {
  cat("Processing species:", sp, "\n")

  ## Subset occurrence data
  sp_df <- alldf %>%
    filter(accepted_name == sp)

  ## List VIF-selected rasters for this species
  sp_vif_path <- paste0(
    "data/04_climate_processing/Cropped/",
    gsub(" ", "_", sp),
    "/VIF/"
  )

  vif_list <- list.files(
    sp_vif_path,
    full.names = TRUE,
    recursive = FALSE
  )

  ## Load rasters
  vifStack <- terra::rast(vif_list)

  ## Run ENMevaluate
  eval <- ENMevaluate(
    occs = sp_df[, c("longitude", "latitude")],
    envs = vifStack,
    tune.args = list(fc = c("L", "Q"), rm = 1:2),
    partitions = "block",
    n.bg = 10000,
    parallel = FALSE,
    algorithm = 'maxent.jar'
  )

  ## Save ENMeval object
  save(eval, file = paste0("data/05_ENMs/", gsub(" ", "_", sp), "_ENM_eval.RData"))

  ## Identify optimal model
  opt.seq <- results %>%
    filter(!is.na(AICc)) %>%
    filter(AICc == min(AICc)) %>%
    filter(or.10p.avg != 0) %>%
    filter(or.10p.avg == min(or.10p.avg)) %>%
    filter(auc.val.avg == max(auc.val.avg))

  ## Save optimal model summary
  write.table(
    opt.seq,
    file = paste0("data/05_ENMs/", gsub(" ", "_", sp), "_OptModel.txt"),
    sep = "\t",
    row.names = FALSE
  )

  ## Plot variable contributions
  opt.mod <- eval.models(eval)[[opt.seq$tune.args]]

  # Save optimal model as RData
  save(opt.mod, file = paste0("data/05_ENMs/", gsub(" ", "_", sp), "_opt_mod.RData"))

  png(filename = paste0("data/05_ENMs/", gsub(" ", "_", sp), "_Variable_Contribution.png"),
      width = 1200,
      height = 800,
      res = 150)

  dismo::plot(opt.mod, main = paste("Variable Contribution -", sp))

  dev.off()

  ## Plot optimal prediction raster
  opt.pred <- eval.predictions(eval)[[as.character(opt.seq$tune.args)]]

  r_df <- as.data.frame(opt.pred, xy = TRUE)
  colnames(r_df) <- c("x", "y", "suitability")

  occ_sf <- st_as_sf(eval@occs, coords = c("longitude", "latitude"), crs = 4326)
  usa <- ne_states(country = "United States of America", returnclass = "sf")

  bbox_vals <- sf::st_bbox(occ_sf)

  p <- ggplot() +
    geom_tile(data = r_df, aes(x = x, y = y, fill = suitability)) +
    scale_fill_viridis_c(name = "Suitability") +
    geom_sf(data = occ_sf, color = "red", size = 0.1) +
    geom_sf(data = usa, fill = NA, color = "white", linewidth = 0.5) +
    coord_sf(
      xlim = c(bbox_vals["xmin"] - 2, bbox_vals["xmax"] + 2),
      ylim = c(bbox_vals["ymin"] - 2, bbox_vals["ymax"] + 2),
      expand = FALSE) +
    labs(
      title = paste("Predicted Suitability -", sp),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey80", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    annotation_north_arrow(
      location = "tl",
      height = unit(1, "cm"),
      width = unit(1, "cm")
    )

  print(p)

  ## Save raster
  writeRaster(x = opt.pred,
              filename = paste0("data/05_ENMs/", gsub(" ", "_", sp), "_ENM_optModel.asc"),
              overwrite = TRUE,
              NAflag = -9999)

  ggsave(filename = paste0("data/05_ENMs/", gsub(" ", "_", sp), "_Optimal_Model_Map.png"),
         plot = p,
         width = 8,
         height = 5,
         dpi = 300)

}
