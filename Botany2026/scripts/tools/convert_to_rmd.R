convert_script <- function(infile, outfile){

  lines <- readLines(infile, warn = FALSE)

  out <- character()

  in_chunk <- FALSE
  chunk <- 1

  open_chunk <- function(label){
    c(
      paste0("```{r ", label, ", eval=FALSE}"),
      ""
    )
  }

  close_chunk <- function(){
    c("", "```", "")
  }

  i <- 1

  while(i <= length(lines)){

    line <- lines[i]

    ## -------------------------------------------------
    ## First title
    ## -------------------------------------------------

    if(i == 1 && grepl("^# ", line)){

      out <- c(out, sub("^# ", "# ", line))
      out <- c(out, "")

      i <- i + 1

      while(i <= length(lines) &&
            grepl("^# ?", lines[i])){

        txt <- sub("^# ?", "", lines[i])

        out <- c(out, txt)

        i <- i + 1
      }

      out <- c(out, "")

      next
    }

    ## -------------------------------------------------
    ## Section heading
    ## -------------------------------------------------

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
            grepl("^# ", lines[i])){

        out <- c(out,
                 sub("^# ?", "", lines[i]))

        i <- i + 1
      }

      out <- c(out, "")

      next
    }

    ## -------------------------------------------------
    ## Code
    ## -------------------------------------------------

    if(!in_chunk){

      label <- paste0("chunk", chunk)

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

files <- list.files(
  "scripts",
  "\\.R$",
  full.names = TRUE
)

for(f in files){

  convert_script(
    infile = f,
    outfile = file.path(
      "book",
      paste0(
        tools::file_path_sans_ext(basename(f)),
        ".Rmd"
      )
    )
  )
}
