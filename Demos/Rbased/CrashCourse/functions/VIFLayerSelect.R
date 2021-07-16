# VIF Layer Select 
## Variable inflation factor based layer selection
## This is based on code by Mike Belitz (mbelitz/Odo_SDM_Rproj)

# Load packages
library(dplyr)
library(usdm)
library(stringr)

## Calculate variable importance from MaxEnt model 
vimportance <- function(max.model){
      m_results <- as.data.frame(as.table(max.model@results)) %>% 
                   dplyr::rename(variables = 1, rem = 2, permutation.importance = 3) %>% 
                   dplyr::select(variables, permutation.importance)
      vIMP <- m_results %>% 
              dplyr::filter(stringr::str_detect(variables, '.permutation.importance')) %>% 
              dplyr::mutate(Variables = stringr::word(variables,  sep = fixed(".")))
      return(vIMP)
}

## Remove one variable 
remove_onevariable <- function(stack, vIMP){
    # Calculate variable inflation factor for current layer stack
    pVIF <- usdm::vif(stack)
    # Left join pVIF with the vIMP from the MaxEnt model
    jdf <- suppressMessages(left_join(pVIF, vIMP))
    # select the variable with the highest VIF and lowest model importance
    lowVar <- jdf %>% 
              dplyr::filter(VIF > 10) %>% 
              dplyr::filter(VIF == sort(VIF, decreasing = TRUE)[1] |
                              VIF == sort(VIF, decreasing = TRUE)[2]) %>% 
              dplyr::filter(permutation.importance == min(permutation.importance))
    # Drop layer from stack
    vv_new <- raster::dropLayer(stack, as.character(lowVar$Variables))
    print(paste0(as.character(lowVar$Variables), " has been removed"))
    return(vv_new)
  }

# VIF layer select
VIF_layerselect <- function(clippedstack, vIMP){
      vv <- clippedstack
      while(max(usdm::vif(vv)$VIF) >= 10){
        vv <- remove_onevariable(stack = vv, vIMP)
      }
      return(vv)
}
