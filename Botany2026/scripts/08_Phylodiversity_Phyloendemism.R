# Phylodiversity and Phyloendemism

### Purpose:
###
### - Use ecological niche models (ENMs) & a pruned phylogenetic tree on 57 Cactaceae species to analyze phylodiversity and phyloendemism in the Sonoran Basin and Range within AZ
### - Construct a P/A matrix for the Cactaceae Community of the Sonoran Basin and Range
### - Use the CANAPER (Categorization of Neo- and PaleoEndemism) R package to evaluate PD, RPD, and Phyloendemism Categories
### - Visualize Basic Species Richness, Observed PD, Significant PD, and Significant RPD measures
### - Identify what makes a spatial cell significant for phylogenetic endemism
### - Visualize the categories of phyloendemism across the Sonoran Basin and Range for Cactaceae


## Load Required Packages ----

library(terra)            # For reading and manipulating raster data (modern replacement for raster package)
library(sf)               # For working with simple features spatial data
library(ggplot2)          # For visualizations
library(ggtree)           # Visualizing the phylogenetic tree
library(ggforce)          # Provides a few plot annotations methods (arrows in this script)
library(patchwork)        # Allows easier viewing of multiple plots
library(rnaturalearth)    # For downloading shapefiles of countries and states
#library(future)           # For running CANAPE fxns in parallel processing (not needed for this demo, but load this in on a personal laptop if you want to use the plan() fxn)
library(canaper)          # Categorical Analysis of Neo- and Paleo- Endemism, also used for significant PD, RPD, PE etc
library(ape)              # Read in Phylogeny, provides pruning fxns as well
library(dplyr)            # Provides data wrangling functions that resemble SQL


## A) Area of interest (AOI) ----
# Load the Sonoran Basin and Range (EPA level 3 ecoregion)
ecoregions <- sf::st_read("data/08_phylodiversity_phyloendemism/Ecoregions/az_eco_l3/az_eco_l3.shp")
son_ecoregion <- ecoregions %>% dplyr::filter(NA_L3NAME %in% c("Sonoran Basin and Range"))
# Load in rnaturalearths sociopolitical boundaries to make it easier to understand where we are
pol_basemap <- rnaturalearth::ne_states(c("United States of America", "Mexico"), returnclass = "sf")
# Simplify to states that surround the Sonoran Basin and Range
son_states <- pol_basemap |> dplyr::filter(name %in% c("Arizona", "California", "Sonora", "Baja California"))
# Reproject to match ecoregion
son_states_reproj <- son_states |> st_transform(crs = st_crs(son_ecoregion))
# Add a few major cities to provide some context
locations <- data.frame(
  name = c("Tuscon", "Phoenix"),
  lon = c(-110.97452583202717, -112.07265839152035),
  lat = c(32.256114232723995, 33.445900912002735)
)
# Make spatial using sf (with wgs84), transform to ecoregion projection
locations_sf <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)
locations_proj <- st_transform(locations_sf, st_crs(son_ecoregion))
# Get bounding box of ecoregion and expand slightly for zoom
bbox <- st_bbox(son_ecoregion)
x_range <- bbox["xmax"] - bbox["xmin"]
y_range <- bbox["ymax"] - bbox["ymin"]
expand <- 0.15  # can be adjusted to zoom out more if we want
# Create nudging so that the names of the cities dont overlap the pts
nudge_x_val <- x_range * 0.05
nudge_y_val <- y_range * 0.05
# Plot!
ggplot() +
  geom_sf(son_states_reproj, mapping = aes(), fill = "white") +
  geom_sf(son_ecoregion, mapping = aes(color = NA_L3NAME), fill = "darkorange", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(labels = c("The Sonoran Basin and Range EPA Ecoregion"), values = c("darkred", "darkorange")) +
  ggtitle("The Sonoran Basin and Range EPA Ecoregion") +
  # add pts of cities
  geom_sf(data = locations_proj, color = "black", size = 3) +
  # add city labels nudged
  geom_sf_text(data = locations_proj, aes(label = name), nudge_x = nudge_x_val, nudge_y = nudge_y_val, color = "black", fontface = "bold", size = 4) +
  coord_sf(crs = st_crs(son_ecoregion),
           xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]))  +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_minimal()

## B) Operational Taxonomic Units (OTUs) P/A data ----
# Load in pre-cropped SDMs for Sonoran Basin and Range Cactaceae
sdm_stack_sonoran <- rast("data/08_phylodiversity_phyloendemism/cactaceae_sdm_stack_sonoran.tif")
# Convert stack to community matrix in dataframe style (cells as sites, species as columns)
comm_matrix_xy <- as.data.frame(sdm_stack_sonoran, xy = TRUE, na.rm = FALSE)
comm_matrix_xy[is.na(comm_matrix_xy)] <- 0 # make NAs 0
# Remove Empty Cells
comm_matrix_xy <- comm_matrix_xy[rowSums(comm_matrix_xy[, -(1:2)]) > 0, ]
# Save coords and assign site ids so we can join these back up later with the canaper output
cell_coords <- data.frame(
  site = paste0("cell_", seq_len(nrow(comm_matrix_xy))),
  x = comm_matrix_xy$x,
  y = comm_matrix_xy$y
)
# Remake community matrix so that coordinate values are removed and rownames are the sites (necessary for canape fxns)
comm_matrix <- comm_matrix_xy[, -(1:2)] # removes coordinate values
rownames(comm_matrix) <- paste0("cell_", seq_len(nrow(comm_matrix)))

# Take a look at the structure of our P/A Matrix?
dim(comm_matrix) # 3967 spatial cells (5x5km), 67 Cactaceae species

## C) Prune a Phylogeny ----

# Read in the Euphyllophyte phylogenetic tree: Carruthers et al. in Review. https://www.biorxiv.org/content/10.64898/2026.01.06.695000v1.full.pdf
molc_plant_tree <- ape::read.tree("data/08_phylodiversity_phyloendemism/pruned-molc-tree-12-10-2025.tre")
# Redo tip labels (currently they have higher order/family info on tip labels)
molc_tip_df <- data.frame(
  original_label = molc_plant_tree$tip.label,
  name = sub(".*_([^_]+_[^_]+)$", "\\1", molc_plant_tree$tip.label) # removes Order Family info
) |>
  mutate(name = gsub("_", " ", name))
molc_plant_tree$tip.label <- molc_tip_df$name # reassign

# Prune the Phylogeny to the taxa in our community matrix
shared_sp_molc <- intersect(molc_plant_tree$tip.label, colnames(comm_matrix)) # find shared taxa among tree and community matrix
pruned_molc_tree <- ape::keep.tip(molc_plant_tree, shared_sp_molc) # prunes the phylogeny to these shared taxa

# Now subset the community matrix so that we're only including taxa that are in the pruned tree
comm_matrix_molc <- comm_matrix[, shared_sp_molc]

# Drop sites that are completely empty post-pruning
comm_matrix_molc <- comm_matrix_molc[rowSums(comm_matrix_molc, na.rm = TRUE) > 0,]
comm_matrix_molc <- comm_matrix_molc[, colSums(comm_matrix_molc, na.rm = TRUE) > 0]

# Check whether the dimensions of sp in community matrix and pruned tree matches
dim(comm_matrix_molc) # sites, sp
length(pruned_molc_tree$tip.label) # number of tips match number of sites, so we're good to go.

## D) Visualize Cactaceae Phylogeny ----

# Prepare tree for plotting
base_tree <- ggtree(pruned_molc_tree, layout = "rectangular", linewidth = 0.25) # set up basic tree plot using ggtree
# Create the final plot with proper axis extension for tip labels
tree_plot <- base_tree +
  geom_tiplab(size = 2.5, align = TRUE, linesize = 0.2) +
  theme_tree2() +
  xlim(NA, max(base_tree$data$x, na.rm = TRUE) * 1.1) + # buffer the axis for viewing tip labels a bit better
  ggtitle("Cactaceae of the Sonoran Basin and Range Phylogeny") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.margin = margin(10, 20, 10, 10), # Increased right margin
    axis.text.x = element_blank(), # removes x-axis
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Print the plot
print(tree_plot)

## Run Phylodiversity Analysis with CANAPER ----
# In canaper's framework, we use randomization tests where we generate a set of random communities, calculate the metric of interest (e.g. PD, RPD) for each random community, then compare the observed values to the random values
# If the observed values are that of the extremes (2.5% for a two-tailed test), they can be considered more or less diverse than random

# canaper generates random communities using the vegan package. To check what randomizations you can use ?vegan::commsim()
# for each randomization algorithm used in canaper, there is two main components:
# 1) n_reps: the number of random communities to simulate (if we have n_reps = 100, we are compared the observed metrics to 100 random replicates)
# 2) n_iterations: the number of swaps used to produce each randomized community (calculation described below)


# compare the percentage of original community matrix similarity to that of the randomized community using 20,000 iterations/swaps
iter_sim_res_molc <- cpr_iter_sim(
  comm = comm_matrix_molc, # our observed community matrix
  null_model = "curveball", # curveball maintains sp richness and abundance patterns, while randomizing species identity. Its a good choice typically but feel free to try others, it can change results!
  n_iterations = 20000, # make an educated guess at number of iterations needed. The more pixels or taxa you have, the higher this usually is
  thin = 100, # frequency to record percentage similarity between original matrix and randomized matrix. Pick something that makes sense with n_iterations
  seed = 123 # seed to make reproducible.
)

# Visualize how many iterations are needed to arrive at a plateau of community dissimilarity for comparing a random community to the observed community
ggplot(iter_sim_res_molc, aes(x = iteration, y = similarity)) +
  geom_line() +
  theme_bw() +
  ggtitle("Number of Iterations to Arrive at Maximum Dissimilarity") +
  labs(x = "Num. iterations", y = "% Similarity") # Appears that maximum dissimilarity between observed community and randomized happens around ~20,000 iterations.

# in regards to n_reps (number of randomized communities to simulate), this isn't as simple. Ideally you could do a sensitivity style analysis by increasing this value and checking where results converge.
# examples of how to do this in a iterative way can be found here:  https://docs.ropensci.org/canaper/articles/canape.html
# keep in mind though that picking a very high number is safer at the expense of computation time, comparing more communities ensures you really are comparing observed vs random

# Run randomization test
#plan(multisession, workers = 7) # set up parallel computing, commented out for workshop. If you are using a personal computer, consider allocating workers to make this faster
set.seed(42)
son_rand_res_molc <- cpr_rand_test(
  comm = comm_matrix_molc, # our community matrix
  phy = pruned_molc_tree, # pruned phylogeny
  null_model = "curveball", # curveball maintains sp richness and abundance patterns, while randomizing species identity. Its a good choice typically but feel free to try others, it can change results!
  n_reps = 50, # number of communities, 50 is a bit low. typically choose 500 or more (but its slow!S)
  n_iterations = 20000, # number of swaps described above, based of cpr_iter_sim() output
  tbl_out = TRUE # creates a tabular output
)
#plan(sequential) # go back to single core mode, disabled for workshop

# Classify Significant PD, RPD, PE, RPE, and Endemism Types (two tailed tests)
son_canape <-
  cpr_classify_endem(son_rand_res_molc) |>
  cpr_classify_signif("pd") |>
  cpr_classify_signif("rpd") |>
  cpr_classify_signif("pe") |>
  cpr_classify_signif("rpe")

## F) Plot Diversity Measures (SR, Observed PD, Significant PD, Significant RPD) ----

# Join classified output w/Sonoran Desert sites for plotting purposes
son_canape_map <- son_canape |>
  full_join(cell_coords, by = "site")

# Reproj our state map & our ecoregion map to match the projection of the outputs
basemap_xy <- sf::st_transform(son_states, crs = st_crs(sdm_stack_sonoran))
son_ecoregion_xy <- sf::st_transform(son_ecoregion, crs = st_crs(sdm_stack_sonoran))

# Create a bbox to narrow in on the AOI
bbox <- st_bbox(son_ecoregion_xy)
x_range <- bbox["xmax"] - bbox["xmin"]
y_range <- bbox["ymax"] - bbox["ymin"]
expand <- 0.15 # adjust if you want this zoomed more out

# Sum across all layers of the SDMs to add up number of species per cell
species_richness <- terra::app(sdm_stack_sonoran, fun = sum, na.rm = TRUE)
# Convert to a df
richness_df <- as.data.frame(species_richness, xy = TRUE)
colnames(richness_df) <- c("x", "y", "richness") # rename columns

# Plot a species richness map
ggplot() +
  geom_tile(data = richness_df, aes(x = x, y = y, fill = richness)) +
  geom_sf(data = son_ecoregion_xy, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = son_states_reproj, fill = NA, color = "black", linewidth = 1) +
  scale_fill_viridis_c(name = "Species\nRichness", option = "plasma") +
  coord_sf(crs = st_crs(sdm_stack_sonoran),
           xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +

  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Sonoran Basin and Range Cactaceae Species Richness Map")

# Plot the observed phylodiversity map (compare against SR, you'll notice some similar trends)
ggplot() +
  geom_tile(data = son_canape_map, aes(x = x, y = y, fill = pd_obs)) +
  scale_fill_viridis_c(name = "Observed\nPhylogenetic Diversity", option = "plasma") +
  geom_sf(data = basemap_xy, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = son_ecoregion_xy, fill = NA, color = "black", linewidth = 1) +
  coord_sf(crs = st_crs(sdm_stack_sonoran),
           xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "Longitude", y = "Latitude", title = "Sonoran Basin and Range Cactaceae Observed Phylogenetic Diversity Map")

# Plot the significant cells for comparing the observed phylodiversity against the randomized communities
pd_sig_map <- ggplot() +
  geom_tile(data = son_canape_map, aes(x = x, y = y, fill = pd_signif)) + # pd_signif contains the two-tailed test output denoting when cells are significant for PD as compared to the random communities
  scale_fill_manual(values = cpr_signif_cols_2, name = "Phylogenetic diversity") +
  geom_sf(data = basemap_xy, fill = NA, color = "black") +
  geom_sf(data = son_ecoregion_xy, fill = NA, color = "black", linewidth = 1) +
  coord_sf(crs = st_crs(sdm_stack_sonoran),
           xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "Longitude", y = "Latitude", title = "Sonoran Basin and Range Cactaceae Significant Phylogenetic Diversity Map")

print(pd_sig_map)

# Plot the significant cells for comparing the alternative tree against the randomized communties
rpd_sig_map <- ggplot() +
  geom_tile(data = son_canape_map, aes(x = x, y = y, fill = rpd_signif)) + # rpd_signif contains the two-tailed test output denoting when cells are significant for RPD as compared to the random communities
  scale_fill_manual(values = cpr_signif_cols_2, name = "Relative Phylogenetic diversity") +
  geom_sf(data = basemap_xy, fill = NA, color = "black") +
  geom_sf(data = son_ecoregion_xy, fill = NA, color = "black", linewidth = 1) +
  coord_sf(crs = st_crs(sdm_stack_sonoran),
           xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "Longitude", y = "Latitude", title = "Sonoran Basin and Range Cactaceae Significant Relative Phylogenetic Diversity")

print(rpd_sig_map)

## G) Phyloendemism Analysis ----
# CANAPE (Categorization of Neo- and Paleoendemism) classifies sites/cells by endemic types (paleo, neo, or mixed) in two steps:
# Step 1) If a cell is significantly high for either the PE on the observed tree (PE_obs) or PE on the alternative tree (PE_alt_obs), we evaluate the cell for neo- or paleoendemism
# Step 2) RPE = PE_observed_tree / PE_alternative_tree. If RPE is significantly high = paleoendemism, if RPE is significantly low = neoendemism. (two tailed test)
# if a cell is significant for PE_obs or PE_alt, but RPE is non-significant, then the cell is assigned to mixed endemism
# a slightly special case of mixed is where PE_obs & PE_obs_alt are both highly significant in a cell (p < 0.01), then its called 'super'

# Plot and visualize all of the cells significant for either PE_obs or PE_alt_obs (one-tailed test)
endemism_significant <- son_canape_map %>%
  mutate(
    PE_significant =
      case_when(
        # If any of the p‑upper/lower values are missing, return NA (logical)
        is.na(pe_obs_p_upper) | is.na(pe_alt_obs_p_upper) |
          is.na(rpe_obs_p_upper) | is.na(rpe_obs_p_lower) ~ NA,
        # If either observed or alternative tree p‑upper exceeds 0.95, call it significant (1 tailed test)
        pe_obs_p_upper > 0.95 | pe_alt_obs_p_upper > 0.95 ~ TRUE,
        # All other cases: not significant
        TRUE ~ FALSE
      )
  )

# Plot sites that are significant for PE in step 1
ggplot() +
  geom_tile(data = endemism_significant, aes(x = x, y = y, fill = PE_significant)) +
  scale_fill_manual(values = c("grey", "purple")) +
  geom_sf(data = basemap_xy, fill = NA, color = "black") +
  geom_sf(data = son_ecoregion_xy, fill = NA, color = "black", linewidth = 1) +
  coord_sf(crs = st_crs(sdm_stack_sonoran),
           xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  guides(fill = guide_legend(title = "Significant Endemism Present")) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "Longitude", y = "Latitude",
       title = "Sonoran Basin and Range Cactaceae Sites that are Significant for Phyloendemism")


# Plot and visualize significant phyloendemism & categories of endemism from CANAPE
ggplot() +
  geom_tile(data = son_canape_map, aes(x = x, y = y, fill = endem_type)) + #  cpr_classify_endem() ran the two tailed test to determine endemism types, so we can just call them here
  scale_fill_manual(values = cpr_endem_cols_4) +
  geom_sf(data = basemap_xy, fill = NA, color = "black") +
  geom_sf(data = son_ecoregion_xy, fill = NA, color = "black", linewidth = 1) +
  coord_sf(crs = st_crs(sdm_stack_sonoran),
           xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  guides(fill = guide_legend(title = "Endemism Type")) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = "Longitude", y = "Latitude",
       title = "Sonoran Basin and Range Cactaceae Categorical Analysis of Neo- and Paleo-Endemism")

## H) Optional: Visualize Lineage Compositions in Significant Cells

# Write a function that visualizes lineage composition within a site, and reports both the observed and expected diversity metric...
create_significant_cell_plot <- function(metric, type, cell_id,
                                         son_canape_map, comm_matrix_molc,
                                         base_tree, bbox, sig_map) {

  # determine significance column and filter values based on metric and type
  signif_col <- ifelse(metric == "pd", "pd_signif", "rpd_signif")

  if (type == "low") {
    signif_values <- c("< 0.01", "< 0.025")
    title_prefix <- "Significantly Low"
    tree_color <- "#E31A1C"  # Red for low
    arrow_from <- "top-right"  # Arrow from top-right corner
  } else if (type == "high") {
    signif_values <- c("> 0.975", "> 0.99")
    title_prefix <- "Significantly High"
    tree_color <- "#1F78B4"  # Blue for high
    arrow_from <- "bottom-left"  # Arrow from bottom-left corner
  } else {
    stop("type must be either 'high' or 'low'")
  }

  # extract the representative cell
  representative_cell <- son_canape_map %>%
    filter(.data[[signif_col]] %in% signif_values & site == cell_id)

  if (nrow(representative_cell) == 0) {
    stop(sprintf("No cell found with %s %s significance and site ID '%s'",
                 title_prefix, metric, cell_id))
  }

  # grab cell coordinates
  cell_x <- representative_cell$x
  cell_y <- representative_cell$y

  # calculate offsets based on plot dimensions (same across both plots)
  x_range <- bbox["xmax"] - bbox["xmin"]
  y_range <- bbox["ymax"] - bbox["ymin"]
  offset_x <- 0.05 * x_range  # 5% of x range as offset
  offset_y <- 0.05 * y_range  # 5% of y range as offset

  # determine arrow start coordinates based on type
  if (arrow_from == "top-right") {
    start_x <- bbox["xmax"] + offset_x
    start_y <- bbox["ymax"] + offset_y
  } else {  # bottom-left
    start_x <- bbox["xmin"] - offset_x
    start_y <- bbox["ymin"] - offset_y
  }

  # create the map visual with arrow pointing to cell
  sig_cell_map <- sig_map +
    # first draw a thicker white arrow as an outline
    geom_segment(data = data.frame(
      x = start_x,
      y = start_y,
      xend = cell_x,
      yend = cell_y
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "white",  # outline color
    size = 1.2,       # thicker for outline effect
    arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +  # Larger arrowhead
    # angle arrow from the specified corner
    geom_segment(data = data.frame(
      x = start_x,
      y = start_y,
      xend = cell_x,
      yend = cell_y
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black",
    size = 0.8,
    arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    ggtitle(paste("Focus on", title_prefix, metric, representative_cell$site, "Pixel"))

  # extract lineages present in representative cell
  comm_matrix_molc_rep <- comm_matrix_molc[representative_cell$site, , drop = FALSE]  # extract by site
  taxa_present_rep <- colnames(comm_matrix_molc_rep)[comm_matrix_molc_rep == 1]  # grab taxa with presence at site

  # create a label to show comparison of observed vs expected metric
  obs_col <- ifelse(metric == "pd", "pd_obs", "rpd_obs")
  exp_col <- ifelse(metric == "pd", "pd_rand_mean", "rpd_rand_mean")

  metric_label <- paste0(
    toupper(substr(metric, 1, 1)), substr(metric, 2, nchar(metric)),
    " observed = ", round(representative_cell[[obs_col]], 3), "\n",
    toupper(substr(metric, 1, 1)), substr(metric, 2, nchar(metric)),
    " expected = ", round(representative_cell[[exp_col]], 3)
  )

  # create groupings using ggtree's groupOTU() to showcase branches
  if (length(taxa_present_rep) > 0) {
    # create a named list for groupOTU using setNames
    group_list <- setNames(list(taxa_present_rep), type)
    tree_grouped <- groupOTU(base_tree$data, group_list)

    # create function for abbreviating Genus + specificEpithet
    abbreviate_taxon <- function(taxon_name) {
      parts <- strsplit(taxon_name, " ")[[1]]
      if (length(parts) >= 2) {
        return(paste0(substr(parts[1], 1, 1), ". ", parts[2]))
      } else {
        return(taxon_name)
      }
    }

    # abbreviate tip labels in the grouped trees
    tree_grouped$label <- sapply(tree_grouped$label, abbreviate_taxon)

    # plot showing contributing branches + metric observed and expected
    values <- c("0" = "grey30")
    values[[type]] <- tree_color

    cell_lineages <- ggtree(tree_grouped, aes(color = group)) +
      scale_color_manual(values = values, guide = "none") +
      geom_tiplab(size = 2.5, align = TRUE, linesize = 0.2) +
      theme_tree2() +
      xlim(NA, max(base_tree$data$x, na.rm = TRUE) * 1.1) +
      ggtitle(paste(title_prefix, metric, representative_cell$site, "Present Lineages")) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.margin = margin(10, 20, 10, 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      annotate("text", x = -Inf, y = Inf, label = metric_label,
               hjust = -0.1, vjust = 1.5,
               size = 4, colour = "black", fontface = "bold")
  } else {
 # catch issues with missing taxa (not sure this really could happen at this stage but whatever)
    warning("No taxa present in the selected cell")
    cell_lineages <- ggtree(base_tree$data) +
      theme_tree2() +
      ggtitle(paste("No taxa present in", representative_cell$site)) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.margin = margin(10, 20, 10, 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  # plot both side by side
  combined_plot <- sig_cell_map + cell_lineages +
    plot_layout(widths = c(1, 2))

  return(combined_plot)
}

# Plot:
# Sig Low PD cell
low_pd_plot <- create_significant_cell_plot(
  metric = "pd",
  type = "low",
  cell_id = "cell_616",
  son_canape_map = son_canape_map,
  comm_matrix_molc = comm_matrix_molc,
  base_tree = base_tree,
  bbox = bbox,
  sig_map = pd_sig_map
)
print(low_pd_plot)

# Sig High PD cell
high_pd_plot <- create_significant_cell_plot(
  metric = "pd",
  type = "high",
  cell_id = "cell_2722",
  son_canape_map = son_canape_map,
  comm_matrix_molc = comm_matrix_molc,
  base_tree = base_tree,
  bbox = bbox,
  sig_map = pd_sig_map
)
print(high_pd_plot)

# Sig Low RPD cell
low_rpd_plot <- create_significant_cell_plot(
  metric = "rpd",
  type = "low",
  cell_id = "cell_507",
  son_canape_map = son_canape_map,
  comm_matrix_molc = comm_matrix_molc,
  base_tree = base_tree,
  bbox = bbox,
  sig_map = rpd_sig_map
)
print(low_rpd_plot)

# Sig High RPD cell
high_rpd_plot <- create_significant_cell_plot(
  metric = "rpd",
  type = "high",
  cell_id = "cell_2050",
  son_canape_map = son_canape_map,
  comm_matrix_molc = comm_matrix_molc,
  base_tree = base_tree,
  bbox = bbox,
  sig_map = rpd_sig_map
)
print(high_rpd_plot)

