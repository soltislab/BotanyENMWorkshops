# Define which scripts and which chunks should be interactive
# Format: list(script_name = c(chunk_numbers))

webr_interactive <- list(
  "00_Setup" = c(),  # No interactive chunks yet
  "01_Data_Downloads" = c(5), # Example: chunk 3 is interactive
  "02_Data_Cleaning" = c(),
  "03_Data_Exploration" = c(),  
  "04_AccessibleArea_VariableSelection" = c(),
  "05_Ecological_Niche_Modeling" = c(),
  "06_Null_ENMs" = c(),
  "07_ENM_Processing" = c(),
  "08_Phylodiversity_Phyloendemism" = c()
)

# Define data files to pre-load for WebR (hidden from users)
webr_data_files <- list(
  "01_Data_Downloads" = list(
    "rawdf" = "data/Shortia_galacifolia_raw_2026_07_14.rds"
    # Add more if needed: "other_data" = "data/other_file.rds"
  )
)
