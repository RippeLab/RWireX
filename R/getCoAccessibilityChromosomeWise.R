#' Get Co-Accessibility to an ArchRProject
#' 
#' This function is an extended version of ArchR package that will add co-accessibility scores between ALL peaks for EACH chromosome in a given ArchRProject. 
#' There are two new modes of choosing cell aggregates: unique and single_cell_resolution. If you are looking for co-accessibility in homogeneous population,
#' where the difference in accessibility is explained by stochastic nature rather than external influence, consider "single_cell_resolution" mode. There is no aggregation in this mode, hence single cell resolution. 
#' Either all cells or only the cells specified in cellsToUse are used.
#' "Unique" mode works similar to the default of ArchR but only allows each cell to be in one aggregate. 
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

getCoAccessibilityChromosomeWise <- function (
    ArchRProj = NULL, 
    reducedDims = "IterativeLSI", 
    dimsToUse = 1:30, 
    scaleDims = NULL, 
    corCutOff = 0.75, 
    cellsToUse = NULL, 
    AggregationMethod = "ArchR_default", 
    numCellsPerAggregate = 100, 
    numAggregates = 500,
    useMatrix = "PeakMatrix",
    overlapCutoff = 0.8, 
    scaleTo = 10^4, 
    log2Norm = TRUE, 
    seed = 1, 
    threads = getArchRThreads(), 
    numPermutations = 1000,
    returnLoops = FALSE,
    verbose = TRUE, 
    logFile = createLogFile("getCoAccessibilityChromosomeWise")
){
  
  myParam <- c(as.list(environment()))
  .checkInputParameters(myParam)
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
  featureSet <- .getSet(ArchRProj, useMatrix)
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
  
  #Create DataFrame For Pairwise Things to Test
  columns = c("queryHits","subjectHits","seqnames", "idx1", "idx2", "correlation", "Variability1", "Variability2", "PercAccess1", "PercAccess2") 
  o <- DataFrame(matrix(nrow = 0, ncol = length(columns)))
  colnames(o) = columns
  
  #Peak Matrix ColSums
  cS <- ArchR:::.getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = useMatrix)
  gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm = TRUE)))
  
  ### Add pseudo-count to gS for non-accessible cell aggregates in all regions of interest
  gS <- gS + 1
  
  
  o <- .addMetadataForAggregates(ArchRProj, o, featureSet, knnObj, useMatrix, gS, log2Norm, scaleTo)
  
  for (x in seq_along(chri)) {
    ArchR:::.logDiffTime(sprintf("Computing Co-Accessibility for complete chromosome %s (%s of %s)", chri[x], x, length(chri)), t1 = tstart, verbose = verbose, logFile = logFile)
    
    peaks_in_one_chrom_idx = which(seqnames(featureSet) == chri[x])
    pairwiseComb = combn(peaks_in_one_chrom_idx, 2)
    pairwiseDF = DataFrame(queryHits=pairwiseComb[1,], subjectHits=pairwiseComb[2,])
    pairwiseDF$seqnames <- chri[x]
    pairwiseDF$idx1 <- featureSet$idx[pairwiseComb[1,]]
    pairwiseDF$idx2 <- featureSet$idx[pairwiseComb[2,]]
    pairwiseDF$correlation <- -999.999
    pairwiseDF$Variability1 <- 0.000
    pairwiseDF$Variability2 <- 0.000
    
    groupMat = .createGroupMatrix(ArchRProj, featureSet, knnObj, useMatrix, gS, log2Norm, chri[x], scaleTo)
    
    #Correlations
    corVals <- ArchR:::rowCorCpp(idxX = pairwiseDF$idx1, idxY = pairwiseDF$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))
    ArchR:::.logThis(head(corVals), paste0("SubsetCorVals-", x), logFile = logFile)
    
    rowVars <- as.numeric(matrixStats::rowVars(groupMat))
    
    ### Calculate percent of accessible cells / cell aggregates for each peak
    numAccessUnits <- matrixStats::rowSums2(groupMat > 0)
    numUnits <- ncol(groupMat)
    percAccessUnits <- numAccessUnits/numUnits*100
    
    pairwiseDF$correlation <- as.numeric(corVals)
    pairwiseDF$Variability1 <- rowVars[pairwiseDF$idx1]
    pairwiseDF$Variability2 <- rowVars[pairwiseDF$idx2]
    
    pairwiseDF$PercAccess1 <- percAccessUnits[pairwiseDF$idx1]
    pairwiseDF$PercAccess2 <- percAccessUnits[pairwiseDF$idx2]
    pairwiseDF$PercAccessMean <- rowMeans(cbind(pairwiseDF$PercAccess1, pairwiseDF$PercAccess2))
    
    ArchR:::.logThis(groupMat, paste0("SubsetGroupMat-", x), logFile = logFile)
    ArchR:::.logThis(pairwiseDF, paste0("SubsetCoA-", x), logFile = logFile)
    o = rbind(o, pairwiseDF)
  }
  
  o <- o[!is.na(o$correlation), ]
  
  o$TStat <- (o$correlation / sqrt((pmax(1-o$correlation^2, 0.00000000000000001, na.rm = TRUE))/(length(knnObj)-2))) #T-statistic P-value
  o$Pval <- 2 * pt(-abs(o$TStat), length(knnObj) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o$VarQuantile1 <- ArchR:::.getQuantiles(o$Variability1)
  o$VarQuantile2 <- ArchR:::.getQuantiles(o$Variability2)
  
  mcols(featureSet) <- NULL
  o@metadata$FeatureSet <- featureSet
  
  o@metadata$SettingsCoAccessibility <- list(reducedDims = reducedDims, dimsToUse = dimsToUse, scaleDims = scaleDims, corCutOff = corCutOff, cellsToUse = cellsToUse,
                                                              AggregationMethod = AggregationMethod, numCellsPerAggregate = numCellsPerAggregate, numAggregates = numAggregates, useMatrix = useMatrix,
                                                              overlapCutoff = overlapCutoff, scaleTo = scaleTo, log2Norm = log2Norm)
  if (returnLoops){
    return(.createLoopsSameChromosome(o))
  }
  ArchR:::.endLogging(logFile = logFile)
  return(o)
}
