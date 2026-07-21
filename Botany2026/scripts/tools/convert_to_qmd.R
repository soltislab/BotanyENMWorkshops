# -------------------------------------------------
# convert_to_qmd.R  –  converts R scripts → QMD for Quarto
# -------------------------------------------------
convert_script <- function(infile, outfile, interactive_chunks = NULL, data_files = NULL){
  lines <- readLines(infile, warn = FALSE)
  out <- character()
  in_chunk <- FALSE
  chunk <- 1
  chunk_prefix <- paste0(
    "script_",
    gsub(
      "[^A-Za-z0-9]+",
      "_",
      tools::file_path_sans_ext(basename(infile))
    )
  )
 # ----- Hidden data‑loading chunk (WebR) -----
if (!is.null(data_files) && length(data_files) > 0) {
  out <- c(
    out,
    "",
    "```{webr-r}",
    "#| include: false",
    "#| context: setup",
    "",
    "# Load pre‑saved data files"
  )
  for (var_name in names(data_files)) {
    file_info <- data_files[[var_name]]

    # ---- NEW: guard against missing `type` ----
    if (is.list(file_info) && !is.null(file_info$type)) {
      # Structured entry with explicit type
      if (identical(file_info$type, "csv")) {
        out <- c(out, paste0(var_name, " <- read.csv('", file_info$file, "')"))
      } else {   # default to RDS (or any other non‑CSV type)
        out <- c(out, paste0(var_name, " <- readRDS('", file_info$file, "')"))
      }
    } else {
      # Simple string path or list without a `type` field → assume RDS
      path <- if (is.list(file_info) && !is.null(file_info$file)) {
        file_info$file
      } else {
        as.character(file_info)   # plain character vector
      }
      out <- c(out, paste0(var_name, " <- readRDS('", path, "')"))
    }
  }

  # Close the webr chunk with a status message
 out <- c(out,
         "# Debug: check if file exists",
         paste0("cat('Looking for: ", path, "\\n')"),
         paste0("cat('File exists: ', file.exists('", path, "'), '\\n')"),
         paste0(var_name, " <- readRDS('", path, "')"))
}
  
  # ----- Chunk open/close helpers -----
  open_chunk <- function(label, chunk_num){
    is_interactive <- chunk_num %in% interactive_chunks
    if(is_interactive){
      c(paste0("```{webr-r}"), "")
    } else {
      c(paste0("```{r ", label, ", eval=FALSE, warning=FALSE, message=FALSE}"), "")
    }
  }
  close_chunk <- function(){ c("", "```", "") }
  # ----- Turn comment lines into markdown -----
  comment_to_md <- function(x){
    txt <- sub("^###\\s?", "", x)
    if(grepl("^\\s*$", txt)) return("")
    txt <- sub("^\\s*-\\s+", "- ", txt)
    txt <- sub("^\\s*([0-9]+\\.)\\s+", "\\1 ", txt)
    txt
  }
  i <- 1
  while(i <= length(lines)){
    line <- lines[i]
    # ----- Title (first line only) -----
    if(i == 1 && grepl("^# ", line)){
      out <- c(out, paste0("# ", sub("^# ", "", line)), "")
      i <- i + 1
      while(i <= length(lines) && grepl("^\\s*$", lines[i])) i <- i + 1
      while(i <= length(lines) && grepl("^###", lines[i])){
        out <- c(out, comment_to_md(lines[i]))
        i <- i + 1
      }
      out <- c(out, "")
      next
    }
    # ----- Section heading -----
    if(grepl("^## .*----$", line)){
      if(in_chunk){ out <- c(out, close_chunk()); in_chunk <- FALSE }
      heading <- sub("^##\\s*", "", line)
      heading <- sub("\\s*----$", "", heading)
      out <- c(out, paste0("## ", heading), "")
      i <- i + 1
      while(i <= length(lines) && grepl("^\\s*$", lines[i])) i <- i + 1
      while(i <= length(lines) && grepl("^###", lines[i])){
        out <- c(out, comment_to_md(lines[i]))
        i <- i + 1
      }
      out <- c(out, "")
      next
    }
    # ----- Code chunk -----
    if(!in_chunk){
      label <- paste0(chunk_prefix, "_chunk", chunk)
      out <- c(out, open_chunk(label, chunk))
      chunk <- chunk + 1
      in_chunk <- TRUE
    }
    out <- c(out, line)
    i <- i + 1
  }
  if(in_chunk) out <- c(out, close_chunk())
  writeLines(out, outfile)
}
# -------------------------------------------------
# Actually run the conversion
# -------------------------------------------------
# Load your WebR configuration (defines which chunks are interactive and which data to pre‑load)
source("Botany2026/scripts/tools/webr_config.R")
# Make sure the output folder exists
dir.create(file.path("book", "chapters"), recursive = TRUE, showWarnings = FALSE)
# Find all R scripts to convert
files <- list.files(
  "Botany2026/scripts",
  pattern = "\\.R$",
  full.names = TRUE
)
# Convert each script → .qmd (note the .qmd extension)
for(f in files){
  script_name <- tools::file_path_sans_ext(basename(f))
  # Get interactive chunks for this script (if any)
  interactive <- if(script_name %in% names(webr_interactive)){
    webr_interactive[[script_name]]
  } else { NULL }
  # Get data files to pre‑load for this script (if any)
  data_files <- if(script_name %in% names(webr_data_files)){
    webr_data_files[[script_name]]
  } else { NULL }
  convert_script(
    infile = f,
    outfile = file.path(
      "book",
      "chapters",
      paste0(script_name, ".qmd")
    ),
    interactive_chunks = interactive,
    data_files = data_files
  )
}
