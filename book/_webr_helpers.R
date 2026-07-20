# Helper function to download data from Zenodo for WebR
download_from_zenodo <- function(zenodo_record_id, filename, local_path = NULL) {
  if (is.null(local_path)) {
    local_path <- filename
  }
  
  url <- sprintf(
    "https://zenodo.org/records/%s/files/%s?download=1",
    zenodo_record_id,
    filename
  )
  
  download.file(url, local_path, quiet = TRUE)
  invisible(local_path)
}
