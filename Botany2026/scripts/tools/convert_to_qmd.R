convert_script <- function(infile, outfile, interactive_chunks = NULL){
  lines <- readLines(infile, warn = FALSE)
  out <- character()
  in_chunk <- FALSE
  chunk <- 1
  chunk_prefix <- paste0(
    "script_",
    gsub("[^A-Za-z0-9]+", "_",
         tools::file_path_sans_ext(basename(infile)))
  )

  # Determine chunk type: heavy computation vs. exploration
  open_chunk <- function(label, chunk_num){
    is_interactive <- chunk_num %in% interactive_chunks

    if(is_interactive){
      # Interactive WebR chunk
      c(paste0("```{webr-r}"), "")
    } else {
      # Static chunk (pre-computed)
      c(paste0("```{r ", label, ", warning=FALSE, message=FALSE, results='hide'}"), "")
    }
  }

  close_chunk <- function(){
    c("", "```", "")
  }

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

    if(i == 1 && grepl("^# ", line)){
      out <- c(out, paste0("# ", sub("^# ", "", line)), "")
      i <- i + 1
      while(i <= length(lines) && grepl("^\\s*$", lines[i])) {
        i <- i + 1
      }
      while(i <= length(lines) && grepl("^###", lines[i])) {
        out <- c(out, comment_to_md(lines[i]))
        i <- i + 1
      }
      out <- c(out, "")
      next
    }

    if(grepl("^## .*----$", line)){
      if(in_chunk){
        out <- c(out, close_chunk())
        in_chunk <- FALSE
      }
      heading <- sub("^##\\s*", "", line)
      heading <- sub("\\s*----$", "", heading)
      out <- c(out, paste0("## ", heading), "")
      i <- i + 1
      while(i <= length(lines) && grepl("^\\s*$", lines[i])) {
        i <- i + 1
      }
      while(i <= length(lines) && grepl("^###", lines[i])) {
        out <- c(out, comment_to_md(lines[i]))
        i <- i + 1
      }
      out <- c(out, "")
      next
    }

    if(!in_chunk){
      label <- paste0(chunk_prefix, "_chunk", chunk)
      out <- c(out, open_chunk(label, chunk))
      chunk <- chunk + 1
      in_chunk <- TRUE
    }
    out <- c(out, line)
    i <- i + 1
  }

  if(in_chunk)
    out <- c(out, close_chunk())

  writeLines(out, outfile)
}

# Load WebR configuration
source("Botany2026/scripts/tools/webr_config.R")

dir.create(
  file.path("book", "chapters"),
  recursive = TRUE,
  showWarnings = FALSE
)

files <- list.files(
  "Botany2026/scripts",
  pattern = "\\.R$",
  full.names = TRUE
)

for(f in files){
  script_name <- tools::file_path_sans_ext(basename(f))
  
  # Get interactive chunks for this script
  interactive <- if(script_name %in% names(webr_interactive)) {
    webr_interactive[[script_name]]
  } else {
    NULL
  }
  
  convert_script(
    infile = f,
    outfile = file.path(
      "book",
      "chapters",
      paste0(script_name, ".qmd")
    ),
    interactive_chunks = interactive
  )
}
