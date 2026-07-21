source("Botany2026/scripts/tools/webr_config.R")

if (!exists("webr_interactive")) {
  webr_interactive <- list()
}

if (!exists("webr_csv_files")) {
  webr_csv_files <- list()
}

doi_to_record_id <- function(doi) {
  sub("^.*zenodo\\.", "", doi)
}

make_webr_setup_chunk <- function(csv_files) {
  if (is.null(csv_files) || length(csv_files) == 0) {
    return(character())
  }

  url_lines <- character()

  for (key in names(csv_files)) {
    info <- csv_files[[key]]
    record_id <- doi_to_record_id(info$doi)
    filename <- info$filename

    url <- sprintf(
      "https://zenodo.org/records/%s/files/%s?download=1",
      record_id,
      filename
    )

    url_lines <- c(
      url_lines,
      paste0("  '", key, "' = '", url, "'")
    )
  }

  url_block <- paste(url_lines, collapse = ",\n")

  c(
    "",
    "```{webr-r}",
    "#| include: false",
    "#| context: setup",
    "",
    "# Hidden WebR helper for CSV files hosted on Zenodo",
    "webr_csv_urls <- list(",
    url_block,
    ")",
    "",
    "webr_read_csv <- function(file, ...) {",
    "  key <- basename(file)",
    "",
    "  if (!key %in% names(webr_csv_urls)) {",
    "    stop(",
    "      paste0(",
    "        'No Zenodo URL registered for ', key, ",
    "        '. Add it to webr_csv_files in webr_config.R.'",
    "      ),",
    "      call. = FALSE",
    "    )",
    "  }",
    "",
    "  local_file <- file",
    "  local_dir <- dirname(local_file)",
    "",
    "  if (!dir.exists(local_dir) && local_dir != '.') {",
    "    dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)",
    "  }",
    "",
    "  if (!file.exists(local_file)) {",
    "    message('Downloading ', key, ' from Zenodo...')",
    "    download.file(",
    "      url = webr_csv_urls[[key]],",
    "      destfile = local_file,",
    "      mode = 'wb',",
    "      quiet = FALSE",
    "    )",
    "  }",
    "",
    "  read.csv(local_file, ...)",
    "}",
    "```",
    ""
  )
}

rewrite_webr_read_csv <- function(line) {
  # Rewrites simple literal-path calls:
  # read.csv("path/to/file.csv")
  # read.csv('path/to/file.csv', ...)
  #
  # into:
  # webr_read_csv("path/to/file.csv")
  gsub(
    "read\\.csv\\((['\"])([^'\"]+)\\1",
    "webr_read_csv(\\1\\2\\1",
    line,
    perl = TRUE
  )
}

convert_script <- function(
  infile,
  outfile,
  interactive_chunks = NULL,
  csv_files = NULL
) {
  lines <- readLines(infile, warn = FALSE)
  out <- character()
  in_chunk <- FALSE
  chunk <- 1
  setup_inserted <- FALSE

  has_interactive <- !is.null(interactive_chunks) &&
    length(interactive_chunks) > 0

  setup_chunk <- if (has_interactive) {
    make_webr_setup_chunk(csv_files)
  } else {
    character()
  }

  chunk_prefix <- paste0(
    "script_",
    gsub(
      "[^A-Za-z0-9]+",
      "_",
      tools::file_path_sans_ext(basename(infile))
    )
  )

  open_chunk <- function(label, chunk_num) {
    is_interactive <- chunk_num %in% interactive_chunks

    if (is_interactive) {
      c(
        "```{webr-r}",
        ""
      )
    } else {
      c(
        paste0(
          "```{r ",
          label,
          ", eval=FALSE, warning=FALSE, message=FALSE}"
        ),
        ""
      )
    }
  }

  close_chunk <- function() {
    c("", "```", "")
  }

  comment_to_md <- function(x) {
    txt <- sub("^###\\s?", "", x)

    if (grepl("^\\s*$", txt)) {
      return("")
    }

    txt <- sub("^\\s*-\\s+", "- ", txt)
    txt <- sub("^\\s*([0-9]+\\.)\\s+", "\\1 ", txt)
    txt
  }

  insert_setup_if_needed <- function(out) {
    if (!setup_inserted && length(setup_chunk) > 0) {
      setup_inserted <<- TRUE
      c(out, setup_chunk)
    } else {
      out
    }
  }

  i <- 1

  while (i <= length(lines)) {
    line <- lines[i]

    # -------------------------------------------------
    # Title
    # -------------------------------------------------
    if (i == 1 && grepl("^# ", line)) {
      out <- c(
        out,
        paste0("# ", sub("^# ", "", line)),
        ""
      )

      i <- i + 1

      while (
        i <= length(lines) &&
          grepl("^\\s*$", lines[i])
      ) {
        i <- i + 1
      }

      while (
        i <= length(lines) &&
          grepl("^###", lines[i])
      ) {
        out <- c(out, comment_to_md(lines[i]))
        i <- i + 1
      }

      out <- c(out, "")
      out <- insert_setup_if_needed(out)

      next
    }

    # If no title was found, insert setup before first content.
    out <- insert_setup_if_needed(out)

    # -------------------------------------------------
    # Section heading
    # -------------------------------------------------
    if (grepl("^## .*----$", line)) {
      if (in_chunk) {
        out <- c(out, close_chunk())
        in_chunk <- FALSE
      }

      heading <- sub("^##\\s*", "", line)
      heading <- sub("\\s*----$", "", heading)

      out <- c(
        out,
        paste0("## ", heading),
        ""
      )

      i <- i + 1

      while (
        i <= length(lines) &&
          grepl("^\\s*$", lines[i])
      ) {
        i <- i + 1
      }

      while (
        i <= length(lines) &&
          grepl("^###", lines[i])
      ) {
        out <- c(out, comment_to_md(lines[i]))
        i <- i + 1
      }

      out <- c(out, "")
      next
    }

    # -------------------------------------------------
    # Code
    # -------------------------------------------------
    if (!in_chunk) {
      label <- paste0(
        chunk_prefix,
        "_chunk",
        chunk
      )

      out <- c(out, open_chunk(label, chunk))
      chunk <- chunk + 1
      in_chunk <- TRUE
    }

    current_chunk <- chunk - 1
    is_interactive <- current_chunk %in% interactive_chunks

    if (is_interactive) {
      line <- rewrite_webr_read_csv(line)
    }

    out <- c(out, line)
    i <- i + 1
  }

  if (in_chunk) {
    out <- c(out, close_chunk())
  }

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

for (f in files) {
  script_name <- tools::file_path_sans_ext(basename(f))

  interactive <- if (script_name %in% names(webr_interactive)) {
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
    interactive_chunks = interactive,
    csv_files = webr_csv_files
  )
}
