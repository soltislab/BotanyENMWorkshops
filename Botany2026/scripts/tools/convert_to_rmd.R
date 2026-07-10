convert_script <- function(infile, outfile){

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

  open_chunk <- function(label){
    c(
      paste0("```{r ", label, ", eval=FALSE}"),
      ""
    )
  }

  close_chunk <- function(){
    c("", "```", "")
  }

  ## ---------------------------------------------
  ## Convert ### comments to Markdown
  ## ---------------------------------------------

  comment_to_md <- function(x){

    txt <- sub("^### ?", "", x)

    # blank documentation line
    if(grepl("^\\s*$", txt))
      return("")

    # bullet list
    txt <- sub("^\\s*-\\s+", "- ", txt)

    # numbered list
    txt <- sub("^\\s*([0-9]+\\.)\\s+", "\\1 ", txt)

    txt
  }

  i <- 1

  while(i <= length(lines)){

    line <- lines[i]

    ## ---------------------------------------------
    ## First title
    ## ---------------------------------------------

    if(i == 1 && grepl("^# ", line)){

      out <- c(out,
               paste0("# ", sub("^# ", "", line)),
               "")

      i <- i + 1

      while(i <= length(lines) &&
            grepl("^###", lines[i])){

        out <- c(out, comment_to_md(lines[i]))

        i <- i + 1
      }

      out <- c(out, "")

      next
    }

    ## ---------------------------------------------
    ## Section headings
    ## ---------------------------------------------

    if(grepl("^## .*----$", line)){

      if(in_chunk){

        out <- c(out, close_chunk())

        in_chunk <- FALSE
      }

      heading <- sub("^##\\s*", "", line)
      heading <- sub("\\s*----$", "", heading)

      out <- c(out,
               paste0("## ", heading),
               "")

      i <- i + 1

      while(i <= length(lines) &&
            grepl("^###", lines[i])){

        out <- c(out, comment_to_md(lines[i]))

        i <- i + 1
      }

      out <- c(out, "")

      next
    }

    ## ---------------------------------------------
    ## R code
    ## ---------------------------------------------

    if(!in_chunk){

      label <- paste0(chunk_prefix, "_chunk", chunk)

      out <- c(out, open_chunk(label))

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

  convert_script(
    infile = f,
    outfile = file.path(
      "book",
      "chapters",
      paste0(
        tools::file_path_sans_ext(basename(f)),
        ".Rmd"
      )
    )
  )
}
