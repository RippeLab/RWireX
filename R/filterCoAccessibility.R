#' Filter the co-accessibility from an ArchRProject
#' 
#' This function filters the co-accessibility by the specified parameters, like correlation cutoff or percent of accessible cells. 
#' Additionally it can change the resolution of peaks to make the plots less convoluted.
#' 
#' When the resolution parameter is set, smaller loops are collapsed into one and the best correlation value is kept for this loop (same way as in the original ArchR). 
#' 
#' @param coAccessibility A coAccessibility DataFrame or Loops created by one of the coAccessibility functions.
#' @param corCutOff A numeric describing the minimum numeric peak-to-peak correlation to return.
#' @param onlyPos A boolean that tells to keep only positive correlation above the specified cutoff. 
#' If set to false, then all correlation bigger or equal to abs(corCutOff) are returned.
#' @param perAccess A numeric describing the minimum percent of accessible cells in both peaks in each pair.
#' @param resolution A numeric describing the bp resolution to use when returning loops. This helps with overplotting of correlated regions.
#' @export



filterCoAccessibility <- function(
    coAccessibilityLoops, 
    corCutOff = NULL, 
    onlyPos = TRUE,
    perAccess = NULL,
    resolution = NULL
){
  ArchR:::.validInput(input = coAccessibilityLoops, name = "CoAccessibility", valid = c("list"))
  ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  ArchR:::.validInput(input = perAccess, name = "perAccess", valid = c("numeric", "null"))
  ArchR:::.validInput(input = onlyPos, name = "onlyPos", valid = c("logical"))
  ArchR:::.validInput(input = resolution, name = "resolution", valid = c("numeric", "null"))
  
  if (!class(coAccessibilityLoops)=="DFrame"){
    #if it's loops
    loops_format = TRUE
    metadata_to_save = coAccessibilityLoops@metadata
    coAccessibilityLoops = coAccessibilityLoops$CoAccessibility
    
    if(is.null(coAccessibilityLoops@elementMetadata@metadata$featureSet)){
      stop("The metadata is not stored in the CoAccessibility object.")
    }
    else{
      featureSet = coAccessibilityLoops@elementMetadata@metadata$featureSet
    }
  }
  else{
    loops_format = FALSE
    featureSet = coAccessibilityLoops@metadata$featureSet
  }
  
  if (!is.null(corCutOff)){
    if (onlyPos){
      coAccessibilityLoops = coAccessibilityLoops[coAccessibilityLoops$correlation >= corCutOff,,drop=FALSE]
    }
    else {
      coAccessibilityLoops = coAccessibilityLoops[abs(coAccessibilityLoops$correlation) >= corCutOff,,drop=FALSE]
    }
    
  }
  if (!is.null(perAccess)){
    coAccessibilityLoops = coAccessibilityLoops[coAccessibilityLoops$PercAccessMean >= perAccess,,drop=FALSE]
  }
  if (loops_format){
    if (!is.null(resolution)){
      peakSummits <- resize(featureSet, 1, "center")
      
      summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
      loops <- ArchR:::.constructGR(
        seqnames = seqnames(peakSummits[coAccessibilityLoops@elementMetadata[,1]]),
        start = summitTiles[coAccessibilityLoops@elementMetadata[,1]],
        end = summitTiles[coAccessibilityLoops@elementMetadata[,2]]
      )
      
      mcols(loops) = mcols(coAccessibilityLoops)
      loops <- loops[order(mcols(loops)$correlation, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      
      loops <- SimpleList(CoAccessibility = loops)
      metadata(loops) = metadata_to_save
      return(loops)
    }
    loops <- SimpleList(CoAccessibility = coAccessibilityLoops)
    metadata(loops) = metadata_to_save
    return(loops)
  }
  else{
    return(coAccessibilityLoops)
  }
}