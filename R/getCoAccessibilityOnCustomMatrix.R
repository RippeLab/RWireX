#' Get Peak Co-Accessibility On Custom Matrix to an ArchRProject
#' 
#' This function combines and extends ArchR::addCoAccessibility and ArchR::getCoAccessibility and will calculate co-accessibility scores between peaks or tiles of a given ArchRProject. 
#' There are two new modes of choosing cell aggregates: unique and single_cell_resolution. If you are looking to investigate gene regulation in cis, we recommend using the "single_cell_resolution" mode. There is no aggregation in this mode, hence single cell resolution. 
#' If you are looking to investigate gene regulation in trans, you can use either "ArchR_default" or the "Unique" aggregation mode. "Unique" works similar to the default of ArchR but only allows each cell to be in one cell aggregate. 
#' Additionally, you can choose between calculating co-accessiblity for peaks or genomic tiles by setting "PeakMatrix" or "TileMatrix". Especially for gene regulation in trans, the tile matrix might be of interest.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param cellsToUse A character vector of cellNames to compute coAccessibility on if desired to run on a subset of the total cells.
#' @param AggregationMethod Choose cell aggregation method from "ArchR_default" (each cell in multiple aggregates), "unique" (each cell in only one aggregate), or "single_cell_resolution" (no cell aggregation, single-cell level).
#' @param numCellsPerAggregate The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param numAggregates The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param useMatrix The name of the data matrix to use for calculation of CoAccessibility. Options include "PeakMatrix" or non-binarized "TileMatrix" (binarize = FALSE in addTileMatrix).
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions from
#' the `featureSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended to keep track
#' of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param numPermutations The number of permutations for computing a distant group of random cells.
#' Increasing this number will not necessarily lead to a more distant group of cells.
#' @param returnLoops A boolean indicating to return the co-accessibility signal as a `GRanges` "loops" object designed for use with
#' the `ArchRBrowser()` or as an `ArchRBrowserTrack()`
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @import parallel
#' @export


getCoAccessibilityOnCustomMatrix <- function (
    ArchRProj = NULL, 
    custom_matrix,
    reducedDims = "IterativeLSI", 
    dimsToUse = 1:30, 
    scaleDims = NULL, 
    corCutOff = 0.75, 
    cellsToUse = NULL, 
    AggregationMethod = "ArchR_default", 
    numCellsPerAggregate = 100, 
    numAggregates = 500,
    useMatrix = "PeakMatrix",
    binaryMatrix = FALSE,
    overlapCutoff = 0.8, 
    maxDist = 100000, 
    scaleTo = 10^4, 
    log2Norm = TRUE, 
    seed = 1, 
    returnLoops = FALSE,
    threads = getArchRThreads(), 
    numPermutations = 1000,
    verbose = TRUE, 
    logFile = createLogFile("getCoAccessibility")
){
  
  myParam <- c(as.list(environment()))
  #.checkInputParameters(myParam)
  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), "addCoAccessibility Input-Parameters", logFile = logFile)
  
  set.seed(seed)
  
  if (is.null(cellsToUse)) {
    cells_number <- nrow(ArchRProj@cellColData)
  } else {
    cells_number <- length(cellsToUse)
  }
  
  checkedParams = .checkNumAggregatesAndCellsPerAggregate(AggregationMethod, numAggregates, numCellsPerAggregate, cells_number)
  numAggregates = checkedParams[[1]]
  numCellsPerAggregate = checkedParams[[2]]
  
  #This set can also can be constructed from tile matrix.
  featureSet <- .getSet(ArchRProj, useMatrix, binaryMatrix)
  rD <- .getFilteredReducedDimensions(ArchRProj, reducedDims, corCutOff, dimsToUse, cellsToUse)
  idx <- .selectCellSeedsForAggregation(ArchRProj, rD, AggregationMethod, numPermutations, numCellsPerAggregate, numAggregates, cellsToUse)
  
  #KNN Matrix
  ArchR:::.logDiffTime(main = "Computing KNN", t1 = tstart, verbose = verbose, logFile = logFile)
  knnObj <- .selectClosestCellsOfCellSeeds(ArchRProj, rD, idx, AggregationMethod, numAggregates, numCellsPerAggregate)
  
  #Determine Overlap and Filter
  ArchR:::.logDiffTime(main = "Identifying Non-Overlapping KNN pairs", t1 = tstart, verbose = verbose, logFile = logFile)
  knnObj <- .filterKNN(knnObj, AggregationMethod, overlapCutoff, numCellsPerAggregate, numAggregates)
  ArchR:::.logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1 = tstart, verbose = verbose, logFile = logFile)
  
  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x) {
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  
  #Check Chromosomes
  chri <- gtools::mixedsort(ArchR:::.availableChr(getArrowFiles(ArchRProj), subGroup = useMatrix))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(featureSet))))
  stopifnot(identical(chri, chrj))
  
  o <- .createPairwiseThingsToTest(featureSet, maxDist)
  
  #Peak Matrix ColSums
  cS <- .getColSumsCustomMatrix(custom_matrix)
  gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm = TRUE)))
  
  ### Add pseudo-count to gS for non-accessible cell aggregates in all regions of interest
  gS <- gS + 1
  
  for (x in seq_along(chri)) {
    
    ArchR:::.logDiffTime(sprintf("Computing Co-Accessibility %s (%s of %s)", chri[x], x, length(chri)), t1 = tstart, verbose = verbose, logFile = logFile)
    
    groupMat = .createCustomGroupMatrix(ArchRProj, custom_matrix, featureSet, knnObj, useMatrix, gS, log2Norm, chri[x], scaleTo)
    
    #Correlations
    idx <- BiocGenerics::which(o$seqnames == chri[x])
    corVals <- ArchR:::rowCorCpp(idxX = o[idx, ]$idx1, idxY = o[idx, ]$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))
    ArchR:::.logThis(head(corVals), paste0("SubsetCorVals-", x), logFile = logFile)
    
    rowVars <- as.numeric(matrixStats::rowVars(groupMat))
    
    ### Calculate percent of accessible cells / cell aggregates for each peak
    numAccessUnits <- matrixStats::rowSums2(groupMat > 0)
    numUnits <- ncol(groupMat)
    percAccessUnits <- numAccessUnits/numUnits*100
    
    o[idx, ]$correlation <- as.numeric(corVals)
    o[idx, ]$Variability1 <- rowVars[o[idx, ]$idx1]
    o[idx, ]$Variability2 <- rowVars[o[idx, ]$idx2]
    
    o[idx, ]$PercAccess1 <- percAccessUnits[o[idx, ]$idx1]
    o[idx, ]$PercAccess2 <- percAccessUnits[o[idx, ]$idx2]
    o[idx, ]$PercAccessMean <- rowMeans(cbind(o[idx, ]$PercAccess1, o[idx, ]$PercAccess2))
    
    ArchR:::.logThis(groupMat, paste0("SubsetGroupMat-", x), logFile = logFile)
    ArchR:::.logThis(o[idx, ], paste0("SubsetCoA-", x), logFile = logFile)
  }
  
  o <- o[!is.na(o$correlation), ]
  
  o$TStat <- (o$correlation / sqrt((pmax(1-o$correlation^2, 0.00000000000000001, na.rm = TRUE))/(length(knnObj)-2))) #T-statistic P-value
  o$Pval <- 2 * pt(-abs(o$TStat), length(knnObj) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o$VarQuantile1 <- ArchR:::.getQuantiles(o$Variability1)
  o$VarQuantile2 <- ArchR:::.getQuantiles(o$Variability2)
  
  
  
  mcols(featureSet) <- NULL
  o@metadata$featureSet <- featureSet
  
  o@metadata$SettingsCoAccessibility <- list(reducedDims = reducedDims, dimsToUse = dimsToUse, scaleDims = scaleDims, corCutOff = corCutOff, cellsToUse = cellsToUse,
                                             AggregationMethod = AggregationMethod, numCellsPerAggregate = numCellsPerAggregate, numAggregates = numAggregates, useMatrix = useMatrix,
                                             overlapCutoff = overlapCutoff, maxDist = maxDist, scaleTo = scaleTo, log2Norm = log2Norm)
  if (returnLoops){
    return(.createLoopsSameChromosome(o))
  }
  ArchR:::.endLogging(logFile = logFile)
  
  return(o)
}
