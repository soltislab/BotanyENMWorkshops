library(tools)

slugify <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("(^-)|(-$)", "", x)
  x
}

is_webr_safe <- function(code) {

  unsupported <- c(
    "terra",
    "sf",
    "rJava",
    "biomod2",
    "ENMeval",
    "ENMTools",
    "dismo",
    "leaflet",
    "rnaturalearth",
    "ggspatial"
  )

  !any(sapply(unsupported, grepl, x = code, fixed = TRUE))
}

convert_script <- function(infile, outdir = "../website") {

  lines <- readLines(infile, warn = FALSE)

  base <- file_path_sans_ext(basename(infile))

  outfile <- file.path(outdir, paste0(base, ".qmd"))

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  out <- character()

  i <- 1

  ## -------------------------
  ## Title
  ## -------------------------

  title <- sub("^#\\s*", "", lines[1])

  out <- c(
    out,
    "---",
    paste0("title: \"", title, "\""),
    "---",
    ""
  )

  ## -------------------------
  ## Intro text
  ## -------------------------

  i <- 2

  while(i <= length(lines) &&
        grepl("^# ?", lines[i])){

    txt <- sub("^# ?", "", lines[i])

    out <- c(out, txt)

    i <- i + 1
  }

  out <- c(out, "")

  section <- "intro"

  while(i <= length(lines)){

    line <- lines[i]

    ## -------------------------
    ## Section heading
    ## -------------------------

    if(grepl("^## .*----$", line)){

      section <- slugify(
        sub("\\s*----$", "",
            sub("^##\\s*", "", line))
      )

      out <- c(
        out,
        "",
        paste0("## ", sub("-", " ", tools::toTitleCase(section))),
        ""
      )

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

    ## -------------------------
    ## Collect code block
    ## -------------------------

    code <- character()

    while(i <= length(lines) &&
          !grepl("^## .*----$", lines[i])){

      code <- c(code, lines[i])

      i <- i + 1
    }

    if(length(code)){

      runnable <- is_webr_safe(code)

      label <- paste(
        slugify(base),
        section,
        sep = "-"
      )

      if(runnable){

        out <- c(
          out,
          paste0("```{webr-r ", label, "}"),
          code,
          "```",
          ""
        )

      } else {

        out <- c(
          out,
          paste0("```{r ", label, "}"),
          "#| eval: false",
          code,
          "```",
          ""
        )

      }
    }
  }

  writeLines(out, outfile)

}

files <- list.files(
  "scripts/",
  pattern = "\\.R$",
  full.names = TRUE
)

lapply(files, convert_script)

