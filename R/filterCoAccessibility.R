#' Filter the co-accessibility from an ArchRProject
#' 
#' This function filters the co-accessibility by the specified parameters, like correlation cutoff or percent of accessible cells. 
#' Additionally it can change the resolution of peaks to make the plots less convoluted.
#' 
#' When the resolution parameter is set, smaller loops are collapsed into one and the best correlation value is kept for this loop (same way as in the original ArchR). 
#' 
#' @param coAccessibility A coAccessibility DataFrame created by one of the coAccessibility functions.
#' @param corCutOff A numeric describing the minimum numeric peak-to-peak correlation to return.
#' @param onlyPos A boolean that tells to keep only positive correlation above the specified cutoff. 
#' If set to false, then all correlation bigger or equal to abs(corCutOff) are returned.
#' @param perAccess A numeric describing the minimum percent of accessible cells in both peaks in each pair.
#' @param resolution A numeric describing the bp resolution to use when returning loops. This helps with overplotting of correlated regions.
#' @export



filterCoAccessibility <- function(
    coAccessibilityLoops = NULL, 
    corCutOff = NULL, 
    onlyPos = TRUE,
    perAccess = NULL,
    resolution = NULL
){
  loops_format = FALSE
  if(is.null(coAccessibilityLoops)){
    
    return(NULL)
  }
  if (!class(coAccessibilityLoops)=="DFrame"){
    #if it's loops
    loops_format = TRUE
    coAccessibilityLoops = coAccessibilityLoops$CoAccessibility
    if(is.null(coAccessibilityLoops@elementMetadata@metadata$featureSet)){
      
      return(NULL)
      
    }
    else{
      featureSet = coAccessibilityLoops@elementMetadata@metadata$featureSet
    }
  }
  else{
    if(is.null(coA_df@metadata$featureSet)){
      
      return(NULL)
      
    }
    featureSet = coA_df@metadata$featureSet
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
  if (!is.null(resolution)){
    peakSummits <- resize(featureSet, 1, "center")
    if(!is.null(resolution)){
      summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
    }else{
      summitTiles <- start(peakSummits)
    }
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
    metadata(loops) = metadata(coAccessibilityLoops)
    return(loops)
  }
  if (loops_format){
    loops <- SimpleList(CoAccessibility = coAccessibilityLoops)
    return(loops)
  }
  else{
    return(coAccessibilityLoops)
  }
}