# Define which scripts have interactive chunks
webr_interactive <- list(
  "00_Setup" = c(),
  "01_Data_Downloads" = c(5),  # Chunk 5 for leaflet preview
  "02_Data_Cleaning" = c(),
  "03_Data_Exploration" = c(),  
  "04_AccessibleArea_VariableSelection" = c(),
  "05_Ecological_Niche_Modeling" = c(),
  "06_Null_ENMs" = c(),
  "07_ENM_Processing" = c(),
  "08_Phylodiversity_Phyloendemism" = c()
)

# Define data files to pre-load for WebR (hidden from users)
# In webr_config.R for "01_Data_Downloads"
webr_data_files <- list(
  "01_Data_Downloads" = list(
    "rawdf" = list(
      file = "../data/Shortia_galacifolia_raw_2026_07_14.rds",  # ← Note the "../"
      type = "rds"
    )
  )
)
