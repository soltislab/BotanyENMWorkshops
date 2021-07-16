# gators functions
## Functions will be part of the R package gatoRs (Geographic And Taxonomic Occurrence R-based Scrubbing).
## N Patten and ML Gaynor

# Load packages
library(dplyr)
library(dplyr) # v.0.8.5
library(tidyr) # v.1.0.2
library(plyr) # v.1.8.6
library(spocc) # v.1.0.8
library(ridigbio) # v.0.3.5
library(tibble) # v.3.0.0
library(rbison)

## fix names
filter_fix_names <- function(df, listofsynonyms, acceptedname){
  df$name <- as.character(df$name)
  optionlist <- listofsynonyms
  print("List of synonyms:")
  print(optionlist)
  # Define search
  search <- ""
  search <- paste(search, optionlist[1], sep = "")
  for(l in 2:length(optionlist)){
    search <- paste(search, optionlist[l], sep = "|")
  }
  search <- paste0(search,"")
  print(search)
  print("Looking at unique names in raw df")
  print(unique(df$name))
  df2 <- dplyr::filter(df, grepl(search, name, ignore.case =  TRUE))
  print("Looking at unique names in cleaned df")
  print(unique(df2$name))
  replacename <- acceptedname
  replacename <- gsub(replacename, pattern = "_",replacement =  " ")
  cleaneddata <- dplyr::mutate(df2, new_name = replacename)
  return(cleaneddata)
}


######################################
# Functions for Downloading Data

#' matchColClasses
#' Source: https://stackoverflow.com/questions/49215193/r-error-cant-join-on-because-of-incompatible-types
matchColClasses <- function(df1, df2) {
  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)
  for (n in sharedColNames) {
    class(df2[, n]) <- sharedColTypes[n]
  }
  return(df2)
}

######################################

# Make name list
list_of_wants <- read.csv("functions/list_of_wants.csv", stringsAsFactors=FALSE)
idigbio_fields <- as.character(list_of_wants$idigbio[1:16])
idigbio_fields_all <- as.character(list_of_wants$idigbio)
gbif_fields <- as.character(list_of_wants$gbif_occurence_raw)
bison_fields <- c(as.character(list_of_wants$bison[1:14]),as.character(list_of_wants$bison[17]))
bison_fields_all <- as.character(list_of_wants$bisonB)
new_names <- as.character(list_of_wants$rename)

# Correct class
correct_class <- function(reduceddataframe){
  reduceddataframe$name <- as.character(reduceddataframe$name)
  reduceddataframe$basis <- as.character(reduceddataframe$basis)
  reduceddataframe$date <- as.character(reduceddataframe$date)
  reduceddataframe$institutionID <- as.character(reduceddataframe$institutionID)
  reduceddataframe$collectionCode <- as.character(reduceddataframe$collectionCode)
  reduceddataframe$collectionID <- as.character(reduceddataframe$collectionID)
  reduceddataframe$country <- as.character(reduceddataframe$country)
  reduceddataframe$county <- as.character(reduceddataframe$county)
  reduceddataframe$state <- as.character(reduceddataframe$state)
  reduceddataframe$locality <- as.character(reduceddataframe$locality)
  reduceddataframe$Latitude <- as.numeric(reduceddataframe$Latitude)
  reduceddataframe$Longitude<- as.numeric(reduceddataframe$Longitude)
  reduceddataframe$ID <- as.character(reduceddataframe$ID)
  reduceddataframe$coordinateUncertaintyInMeters <- as.character(reduceddataframe$coordinateUncertaintyInMeters)
  reduceddataframe$informationWithheld <- as.character(reduceddataframe$informationWithheld)
  reduceddataframe$habitat <- as.character(reduceddataframe$habitat)
  reduceddataframe$prov <- as.character(reduceddataframe$prov)
  return(reduceddataframe)
}

# getidigbio
getidigbio <- function(synonyms_list){
  spocc_idigbio <- ridigbio::idig_search_records(rq = list(scientificname=synonyms_list),
                                                 fields = idigbio_fields)
  
  if(nrow(spocc_idigbio) > 0){
    spocc_idigbio_reduced <- spocc_idigbio  %>%
      mutate(prov = "idigbio") %>%
      dplyr::select(tidyselect::all_of(idigbio_fields_all))
    colnames(spocc_idigbio_reduced) <- new_names
    spocc_idigbio_reduced$Latitude <- as.numeric(spocc_idigbio_reduced$Latitude)
    spocc_idigbio_reduced$Longitude <- as.numeric(spocc_idigbio_reduced$Longitude)
    spocc_idigbio_reduced$coordinateUncertaintyInMeters <- as.character(spocc_idigbio_reduced$coordinateUncertaintyInMeters)
    return(spocc_idigbio_reduced)
  } else if (nrow(spocc_idigbio) == 0){
    return(spocc_idigbio)
  }
}


# getgbif
getgbif <- function(spocc_query){
  spocc_asgbif <- spocc_query$gbif
  spocc_gbif <- occ2df(spocc_asgbif)
  if(nrow(spocc_gbif) == 0){
    print(paste0("spocc finds no occurrence records for ", spocc_query$gbif$meta$opts$scientificName))
  } else if (nrow(spocc_gbif) > 0){
    spocc_gbif <- check_columns(spocc_gbif, gbif_fields)
    spocc_gbif_reduced <- dplyr::select(spocc_gbif, tidyselect::all_of(gbif_fields))
    spocc_gbif_reduced <- as.data.frame(spocc_gbif_reduced)
    colnames(spocc_gbif_reduced) <- new_names
    spocc_gbif_reduced <- correct_class(spocc_gbif_reduced)
    return(spocc_gbif_reduced)
  }}

# getbison
getbison <- function(spocc_query){
  spocc_asbison <- (spocc_query$bison)
  spocc_bison <- occ2df(spocc_asbison)
  if(nrow(spocc_bison) == 0){
    print(paste0("spocc finds no occurrence records for ", spocc_query$bison$meta$opts$scientificName))
  } else if (nrow(spocc_bison) > 0){
    spocc_bison <- check_columns(spocc_bison, bison_fields_all)
    spocc_bison_reduced <- dplyr::select(spocc_bison, tidyselect::all_of(bison_fields_all)) 
    spocc_bison_reduced <- as.data.frame(spocc_bison_reduced)
    colnames(spocc_bison_reduced) <- new_names
    spocc_bison_reduced <- correct_class(spocc_bison_reduced)
    return(spocc_bison_reduced)
  }}

# spocc_combine
spocc_combine <- function(synonyms_list, newfilename){
  bison_level <- list()
  for(i in 1:length(synonyms_list)){
    bison <-  rbison::bison_tax(query = synonyms_list[i])
    bison_level[i] <- bison$numFound
  }
  bison_level <- Reduce(`+`, bison_level)
  if(bison_level == 0){
    spocc_query <- occ(query = synonyms_list,
                       from = c('gbif', 'idigbio'),
                       has_coords = NULL)
    spocc_query_df <- occ2df(spocc_query)
    query_idigbio <- getidigbio(synonyms_list)
    if(nrow(spocc_query_df) == 0 & nrow(query_idigbio) == 0){
      print(paste0("Species ", synonyms_list, " is not found"))
    } else if(nrow(spocc_query_df) > 0){
      spocc_query_df <- spocc_query_df %>%
        dplyr::select(spocc.latitude  = latitude,
                      spocc.longitude = longitude,
                      ID = key,
                      spocc.prov = prov,
                      spocc.date = date,
                      spocc.name = name)
      # Query parts
      query_gbif <- getgbif(spocc_query)
      
      # Join
      if(nrow(query_idigbio) == 0 & class(query_gbif) != "character"){
        query_combinedA <- query_gbif
      } else if(nrow(query_idigbio) > 0 & class(query_gbif) == "character"){
        query_combinedA <- query_idigbio
      } else if(nrow(query_idigbio) > 0){
        query_combinedA <- rbind(query_idigbio, query_gbif)
      } 
      
      query_combinedB <- dplyr::left_join(query_combinedA, spocc_query_df, by = "ID" )
      # Write as csv
      write.csv(query_combinedB, newfilename, row.names = FALSE)
    }else if (nrow(query_idigbio) > 0 & nrow(spocc_query_df) == 0) {
      query_combinedA <- query_idigbio
      # Write as csv
      write.csv(query_combinedA, newfilename, row.names = FALSE)
      
    }
  } else if(bison_level > 0){
    spocc_query <- occ(query = synonyms_list,
                       from = c('gbif', 'bison', 'idigbio'),
                       has_coords = NULL, throw_warnings = FALSE)
    spocc_query_df <- occ2df(spocc_query)
    
    # Query parts
    query_idigbio <- getidigbio(synonyms_list)
    if(nrow(spocc_query_df) == 0 & nrow(query_idigbio) == 0){
      print(paste0("Species ", synonyms_list, " is not found"))
      
    } else if(nrow(spocc_query_df) > 0){
      spocc_query_df <- spocc_query_df %>%
        dplyr::select(spocc.latitude  = latitude,
                      spocc.longitude = longitude,
                      ID = key,
                      spocc.prov = prov,
                      spocc.date = date,
                      spocc.name = name)
      query_gbif <- getgbif(spocc_query)
      query_bison <- getbison(spocc_query)
      if(nrow(query_idigbio) == 0 & class(query_gbif) != "character" & class(query_bison) =="character"){
        query_combinedA <- query_gbif
        query_combinedB <- dplyr::left_join(query_combinedA, spocc_query_df, by = "ID" )
        
      } else if(nrow(query_idigbio) == 0 & class(query_gbif) != "character" & class(query_bison) =="character"){
        query_combinedA <- rbind(query_gbif, query_bison)
        query_combinedB <- dplyr::left_join(query_combinedA, spocc_query_df, by = "ID" )
        
      }else if(nrow(query_idigbio) == 0 & class(query_gbif) == "character" & class(query_bison) !="character"){
        query_combinedA <- query_bison
        query_combinedB <- dplyr::left_join(query_combinedA, spocc_query_df, by = "ID" )
        
      }else if(nrow(query_idigbio) > 0 & class(query_gbif) == "character" & class(query_gbif) == "character"){
        query_combinedB <- query_idigbio
      } else if(nrow(query_idigbio) > 0 & class(query_gbif) == "character" & class(query_gbif) != "character"){
        query_combinedA <- rbind(query_idigbio, query_bison)
        query_combinedB <- dplyr::left_join(query_combinedA, spocc_query_df, by = "ID" )
        
      } else if(nrow(query_idigbio) > 0 & class(query_gbif) != "character" & class(query_gbif) != "character"){
        query_combinedA <- rbind(query_idigbio, query_gbif)
        query_combinedB <- rbind(query_bison, query_combinedA)
        query_combinedB <- dplyr::left_join(query_combinedB, spocc_query_df, by = "ID" )
        
      } 
      # Write as csv
      write.csv(query_combinedB, newfilename, row.names = FALSE)
    } else if(nrow(query_idigbio) > 0 & nrow(spocc_query_df) == 0){
      query_combinedA <- query_idigbio
      # Write as csv
      write.csv(query_combinedA, newfilename, row.names = FALSE, )
    }
  }
}


# Check columns
check_columns <- function(spocc_name, fields){
  diff1 <- setdiff(fields, colnames(spocc_name))
  newframe <- data.frame(matrix(, nrow = 1 , ncol= as.numeric(length(diff1))))
  colnames(newframe) <- diff1
  spocc_name_new <- cbind(spocc_name, newframe)
  return(spocc_name_new)
}


