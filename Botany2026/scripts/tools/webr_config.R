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
