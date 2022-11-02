#' Filter the co-accessibility from an ArchRProject
#' 
#' This function filters the co-accessibility by the specified parameters, like correlation cutoff or percent of accessible cells. 
#' Additionally it can change the resolution of peaks to make the plots less convoluted.
#' 
#' @param coAccessibility A coAccessibility DataFrame created by one of the coAccessibility functions.
#' @param corCutOff A numeric describing the minimum numeric peak-to-peak correlation to return.
#' @param perAccess A numeric describing the minimum percent of accessible cells in both peaks in each pair.
#' @param resolution A numeric describing the bp resolution to use when returning loops. This helps with overplotting of correlated regions.
#' @export



filterCoAccessibility <- function(
    coAccessibility = NULL, 
    corCutOff = NULL, 
    perAccess = NULL,
    resolution = NULL
){
  
  if(is.null(coAccessibility)){
    
    return(NULL)
  }
  
  
  if(is.null(metadata(coAccessibility)$featureSet)){
    
    return(NULL)
    
  }
  
  if (!is.null(corCutOff)){
    coAccessibility = coAccessibility[coAccessibility$correlation >= corCutOff,,drop=FALSE]
  }
  if (!is.null(perAccess)){
    coAccessibility = coAccessibility[coAccessibility$PercAccess1 >= perAccess & coAccessibility$PercAccess2 >= perAccess,,drop=FALSE]
  }
  if (!is.null(resolution)){
    peakSummits <- resize(metadata(coAccessibility)$featureSet, 1, "center")
    if(!is.null(resolution)){
      summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
    }else{
      summitTiles <- start(peakSummits)
    }
    loops <- ArchR:::.constructGR(
      seqnames = seqnames(peakSummits[coAccessibility[,1]]),
      start = summitTiles[coAccessibility[,1]],
      end = summitTiles[coAccessibility[,2]]
    )
    coAccessibility$Loops = loops
    coAccessibility = coAccessibility[!duplicated(coAccessibility$Loops),]
  }
  return(coAccessibility)
  
}