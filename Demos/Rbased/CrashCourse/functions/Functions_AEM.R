################################################## FUNCTIONS ##################################################
### AE Melton
# Run this code to generate a function for making a binary map of predicted occurrences
make.binary.map <- function(model, occ.dat){
  
  ###Extract suitability scores
  SuitabilityScores <- raster::extract(model, occ.dat[,3:2])
  ###Get rid of NAs
  SuitabilityScores <- SuitabilityScores[complete.cases(SuitabilityScores)]
  ###Reclassify the raster; set threshold to minimum suitability score at a known occurrence
  threshold <- min(SuitabilityScores)
  
  M <- c(0, threshold, 0,  threshold, 1, 1); 
  
  rclmat <- matrix(M, ncol=3, byrow=TRUE); 
  
  Dist <- raster::reclassify(model, rcl = rclmat);
}

# The following code creates a function to calculate hypervolumes
get_hypervolume <- function(binary_projection, envt) {
  dist.points <-  rasterToPoints(binary_projection)#Need predicted occurrence points (calculated from thresholded model)
  hv.dat <- raster::extract(envt, dist.points[,1:2]);
  hv.dat <- hv.dat[complete.cases(hv.dat),];
  hv.dat <- scale(hv.dat, center=TRUE, scale=TRUE)
  hv <- hypervolume(data = hv.dat, method = "box")
}

##################################ENMTOOLS functions########################################################
#' Calculates the correlation coefficient between two rasters.
#'
#' @param x Either a raster or an ENMTools model object with a suitability raster.
#' @param y Either a raster or an ENMTools model object with a suitability raster.
#' @param method The method to be used for calculating correlations.  Defaults to spearman, but can take "kendall" or "pearson" as well.
#'
#' @examples
#' data(euro.worldclim)
#' raster.cor(euro.worldclim[[1]], euro.worldclim[[2]])


raster.cor <- function(x, y, method="spearman"){

  if(inherits(x, "enmtools.model")){
    x <- x$suitability
  }

  if(inherits(y, "enmtools.model")){
    y <- y$suitability
  }


  df <- cbind(getValues(x), getValues(y))

  df <- df[complete.cases(df),]

  return(cor(df[,1], df[,2], method=method))
}

###
#' raster.standardize, standardizes all values in a raster file
#'
#' This function is used by a lot of the metrics calculated by ENMTools, in order to
#' standardize suitability scores so they sum to 1 over a geographic space.
#'
#'
#' @param x A raster or RasterLayer object, or an ENMTools model object containing a suitability raster.
#' @param verbose Controls printing of diagnostic messages
#'
#'
#' @keywords keywords
#'
#' @examples
#' data(euro.worldclim)
#' raster.standardize(euro.worldclim[[1]])


raster.standardize <- function(x, verbose=FALSE){

  if(inherits(x, "enmtools.model")){
    x <- x$suitability
  }

  if(verbose){
    print(paste("Starting standardize on", x, "at", Sys.time()))
  }

  return(x/cellStats(x, stat=sum))
}

###
#' raster.overlap, measures overlap between two ENMs
#'
#' This function measures similarity in the geographic distribution of suitability scores
#' from two ENMs.  It returns two metrics, I and D.  These metrics are described in
#' Warren et al. 2008.
#'
#' @param x A raster or RasterLayer object, or ENMTools model object containing a suitability raster.
#' @param y Another raster or RasterLayer object, or ENMTools model object containing a suitability raster.
#' @param verbose Controls printing of diagnostic messages
#'
#' @return results A vector containing the three metrics (I, D, and Spearman rank correlation)
#'
#' @keywords keywords
#'
#' @examples
#' data(iberolacerta.clade)
#' data(euro.worldclim)
#' aurelioi.glm <- enmtools.glm(iberolacerta.clade$species$aurelioi,
#' euro.worldclim, f = pres ~ bio1 + bio12)
#' aranica.glm <- enmtools.glm(iberolacerta.clade$species$aranica,
#' euro.worldclim, f = pres ~ bio1 + bio12)
#' raster.overlap(aurelioi.glm, aranica.glm)

raster.overlap <- function(x, y, verbose=FALSE){

  if(any(grepl("enmtools", class(x)))){
    x <- x$suitability
  }

  if(any(grepl("enmtools", class(y)))){
    y <- y$suitability
  }

  if(verbose){
    print(paste("Starting overlap at", Sys.time()))
  }

  x <- raster.standardize(x)
  y <- raster.standardize(y)

  D <- 1 - cellStats(abs(x - y), stat=sum)/2
  I <- 1 - cellStats((sqrt(x) - sqrt(y))^2, stat=sum)/2
  rank.cor <- raster.cor(x, y)

  results <- list(D = D, I = I, rank.cor = rank.cor)
  return(results)

}


