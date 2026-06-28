# Data_Exploration_Updated.R
# ---------------------------------------------------------------
# Purpose: Analyze climatic variables associated with species occurrences using PCA and ANOVA.
# Includes PCA loadings extraction and plots for top contributors to PC1 and PC2.
# ---------------------------------------------------------------

## ---- Load Required Packages ----
library(gtools)         # For mixedsort of file names
library(terra)          # For raster and spatial extraction
library(dplyr)          # For data wrangling
library(tidyr)          # For handling NAs and reshaping
library(ggplot2)        # For plotting
library(ggbiplot)       # For PCA biplot visualization
library(multcompView)   # For compact letter display (Tukey HSD)
library(gridExtra)      # For combining plots

## ---- A) Load Cleaned Data ----

# Read in cleaned occurrence data with geographic coordinates and species ID
alldf <- read.csv("data/02_cleaning/maxent_ready/diapensiaceae_maxentready_2025_06_27.csv")

## ---- B) Load and Stack Climatic Variables from WorldClim ----

# List and sort all .tif files in the WorldClim directory
list <- list.files("data/04_climate_processing/BioClim/", full.names = TRUE)
list <- gtools::mixedsort(sort(list))

# Stack the raster files
envtStack <- terra::rast(list)

## ---- C) Extract Climatic Variables at Occurrence Points ----

# Extract values of environmental rasters at occurrence points
ptExtracted <- terra::extract(envtStack, alldf[, c("longitude", "latitude")])

# Combine with species names and coordinates
ptExtracteddf <- ptExtracted %>%
  dplyr::mutate(name = as.character(alldf$accepted_name), 
                x = alldf$longitude, 
                y = alldf$latitude)

# Remove any rows with missing climate values
ptExtracteddf <- ptExtracteddf %>% tidyr::drop_na()

## ---- D) Perform PCA on Climatic Variables ----

# Subset climate variables (assumed to be columns 2 to 21)
data.bioclim <- ptExtracteddf[, 2:21]
data.species <- ptExtracteddf[, "name"]

# Scale and run PCA
pca_result <- prcomp(data.bioclim, scale. = TRUE)

# Extract variable loadings
loadings <- pca_result$rotation

# Calculate relative contributions of each variable to each PC
loadings_relative <- sweep(abs(loadings), 2, colSums(abs(loadings)), "/") * 100

## ---- E) Visualize PCA ----

# Customize theme
custom_theme <- theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line"),
  axis.text = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  text = element_text(size = 12)
)

# Define color palette for species
pal <- c("#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF")

# PCA Biplot
pca_plot <- ggbiplot(pca_result, obs.scale = 1, var.scale = 1,
                     groups = data.species, ellipse = TRUE, circle = TRUE) +
  scale_color_manual(name = '', values = pal) +
  theme(legend.direction = 'vertical', legend.position = 'bottom',
        legend.text = element_text(size = 12, face = "italic")) +
  custom_theme

print(pca_plot)

## ---- F) Identify Top Contributors to PC1 and PC2 ----

# Top 2 contributing variables to PC1 and PC2
top_PC1_vars <- rownames(loadings)[order(abs(loadings[, "PC1"]), decreasing = TRUE)][1:2]
top_PC2_vars <- rownames(loadings)[order(abs(loadings[, "PC2"]), decreasing = TRUE)][1:2]
selected_vars <- unique(c(top_PC1_vars, top_PC2_vars))  # remove duplicates

## ---- G) ANOVA + Tukey HSD + Plotting for Top Variables ----

plotlist <- list()

for (i in seq_along(selected_vars)) {
  varname <- selected_vars[i]
  
  # Run ANOVA
  model <- aov(as.formula(paste0(varname, " ~ name")), data = ptExtracteddf)
  tukey <- TukeyHSD(model)
  pvals <- tukey$`name`[, 4]
  
  # Convert p-values to significance letters
  letters_out <- multcompLetters(pvals)
  groups_df <- data.frame(
    name = names(letters_out$Letters),
    groups = letters_out$Letters
  )
  
  # Set y position for labels
  max_y <- max(ptExtracteddf[[varname]], na.rm = TRUE)
  groups_df$y_position <- max_y * 1.1
  
  # Create plot
  plotlist[[i]] <- ggplot(ptExtracteddf, aes(x = name, y = .data[[varname]])) +
    geom_jitter(aes(color = name), width = 0.3, alpha = 0.5, size = 1) +
    geom_boxplot(fill = NA, color = "black", size = 0.6) +
    geom_text(data = groups_df,
              aes(x = name, y = y_position, label = groups),
              inherit.aes = FALSE, size = 4) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    ggtitle(paste0(varname)) +
    ylab(varname) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
    ylim(0, max_y * 1.3)
}

# Combine plots
gridExtra::grid.arrange(grobs = plotlist, ncol = 2)
