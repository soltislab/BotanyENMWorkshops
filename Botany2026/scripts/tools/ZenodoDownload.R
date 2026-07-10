#' @title Zenodo Download
#' @description This function was made to download all files needed to run workshop.
#'  It downloads a zipped file from Zenodo and unzips it to
#'  the working directory  within the workshop folder.
#'
#' @param overwrite Logical. If TRUE, the function will overwrite any existing files
#' in the data directory. Default is FALSE.
#' @examples
#' if(exists("crazy")){
#'    SetupBasicExample(overwrite = FALSE)
#' }
#' @returns This function downloads files to the users working directory and will
#' print messages as it downloads to keep you updated on the progress.
#'
#' @importFrom tools R_user_dir
#' @importFrom httr2 request req_perform resp_url req_timeout
#' @importFrom utils unzip
#' @import magrittr
#' @export

ZenodoDownload <- function(overwrite = FALSE){

  ## Setup directories
  dbdir <- getwd()
  if(grepl("BotanyENMWorkshops/Botany2026", dbdir) == FALSE){
    stop("User must set the working directory! Use setwd() to set your directory to BotanyENMWorkshops/Botany2026, or open the Rproject!", call. = TRUE)
  }else{
    message(paste0("Starting data download in ", dbdir))
  }


  ## Check if downloaded prior
  if(dir.exists(paste0(dbdir, "/data/")) == TRUE){
    if(overwrite == FALSE){
      message("The data has already been downloaded. You are ready to rumble!")
      token <- 0
    }else{
      message("The data has already been downloaded, but you have chosen to overwrite it.")
      token <- 1
    }
  }else{
    message("Downloading data for the first time. This may take a few minutes.")
    token <- 1
  }
  ## Download and unzip
  if(token == 1){
    message("Checking Zenodo DOI.")
    ### Zenodo DOI
    doi <- "10.5281/zenodo.16755492" ###UPDATE ME###

    landing_page <- tryCatch(
      httr2::request(paste0("https://doi.org/", doi)) |>
        httr2::req_perform() |>
        httr2::resp_url(),
      error = function(e) {
        stop("Could not resolve DOI ", doi, ".", call. = FALSE)
      }
    )

    # Extract record ID
    record_id <- sub(".*/records/([0-9]+).*", "\\1", landing_page)

    # Query Zenodo API
    url <- sprintf(
      "https://zenodo.org/records/%s/files/data.zip?download=1",
      record_id)

    message("Starting download of data.zip from Zenodo. This may take a few minutes.")
    zipfile <- tempfile(fileext = ".zip")

    tryCatch(
      httr2::request(url) |>
        httr2::req_timeout(seconds = 3600) |>
        httr2::req_perform(path = zipfile),
      error = function(e) {
        stop("Download from Zenodo failed:\n", conditionMessage(e),
             call. = FALSE)
      }
    )

    if (!file.exists(zipfile)) {
      stop("Download failed: zip file was not created.",
           call. = FALSE)
    }
    message("Zipfile Downloaded! Time to unzip!")
    utils::unzip(zipfile,
                 exdir = dbdir,
                 overwrite = TRUE)

    if(!dir.exists(paste0(dbdir, "/data/"))){
      stop("File did not unzip properly.",  call. = FALSE)
    }

    # Clean up
    unlink(zipfile)
    message("Done!")

  }else{
    message("Done!")
  }

}
