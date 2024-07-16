#' getBackgroundCoAccessibility
#' 
#' This function returns background co-accessibility scores for peaks of a given ArchRProject. The background co-accessibility is determined by shuffling of features and cells.
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
#' @param useMatrix The name of the data matrix to use for calculation of BackgroundCoAccessibility. Options include "PeakMatrix" or non-binarized "TileMatrix" (binarize = FALSE in addTileMatrix).
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
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @import parallel
#' @export


getBackgroundCoAccessibility <- function(
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
    binaryMatrix = FALSE,
    overlapCutoff = 0.8, 
    maxDist = 100000,
    chromosomewise=FALSE,
    scaleTo = 10^4,
    log2Norm = TRUE,
    seed = 1, 
    threads = getArchRThreads(),
    numPermutations = 1000,
    verbose = TRUE,
    logFile = createLogFile("getBackgroundCoAccessibility")
    ){
  
    myParam <- c(as.list(environment()))
    .checkInputParameters(myParam)

    tstart <- Sys.time()
    ArchR:::.startLogging(logFile = logFile)
    ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "getBackgroundCoAccessibility Input-Parameters", logFile = logFile)

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

    #Create Ranges
    peakSummits <- resize(featureSet, 1, "center")
    peakWindows <- resize(peakSummits, maxDist, "center")
    
    if (chromosomewise){
      columns = c("queryHits","subjectHits","seqnames", "idx1", "idx2", "correlation", "Variability1", "Variability2", "PercAccess1", "PercAccess2") 
      o <- DataFrame(matrix(nrow = 0, ncol = length(columns)))
      colnames(o) = columns
    }
    else{
      o <- .createPairwiseThingsToTest(featureSet, maxDist)
    }
    
    
    o_featShuffle <- o
    o_cellShuffle <- o


    #Peak Matrix ColSums
    cS <- ArchR:::.getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = useMatrix)
    gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))

    #Add pseudo-count to gS for non-accessible cell aggregates in all regions of interest
    gS <- gS + 1


    for(x in seq_along(chri)){

        ArchR:::.logDiffTime(sprintf("Computing Background Co-Accessibility %s (%s of %s)", chri[x], x, length(chri)), t1=tstart, verbose=verbose, logFile=logFile)
        
        if (chromosomewise){
          peaks_in_one_chrom_idx = which(seqnames(featureSet) == chri[x])
          pairwiseComb = combn(peaks_in_one_chrom_idx, 2)
          pairwiseDF = DataFrame(queryHits=pairwiseComb[1,], subjectHits=pairwiseComb[2,])
          pairwiseDF$seqnames <- chri[x]
          pairwiseDF$idx1 <- featureSet$idx[pairwiseComb[1,]]
          pairwiseDF$idx2 <- featureSet$idx[pairwiseComb[2,]]
          pairwiseDF$correlation <- -999.999
          pairwiseDF$Variability1 <- 0.000
          pairwiseDF$Variability2 <- 0.000
          o = rbind(o, pairwiseDF)
        }
       
        groupMat = .createGroupMatrix(ArchRProj, featureSet, knnObj, useMatrix, gS, log2Norm, chri[x], scaleTo)
        
        ###Shuffle features for background correlations
        groupMat_featShuffle <- apply(groupMat, 2, function(x){sample(x)}) ### shuffle features
        groupMat_cellShuffle <- t(apply(groupMat, 1, function(x){sample(x)})) ### shuffle cells

        #Correlations
        idx <- BiocGenerics::which(o$seqnames==chri[x])
        
        corVals_featShuffle <- ArchR:::rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat_featShuffle), Y = as.matrix(groupMat_featShuffle))
        corVals_cellShuffle <- ArchR:::rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat_cellShuffle), Y = as.matrix(groupMat_cellShuffle))
        ArchR:::.logThis(head(corVals_featShuffle), paste0("SubsetCorVals-", x), logFile = logFile)
        ArchR:::.logThis(head(corVals_cellShuffle), paste0("SubsetCorVals-", x), logFile = logFile)

        rowVars_featShuffle <- as.numeric(matrixStats::rowVars(groupMat_featShuffle))
        rowVars_cellShuffle <- as.numeric(matrixStats::rowVars(groupMat_cellShuffle))
        
        numAccessUnits_featShuffle <- matrixStats::rowSums2(groupMat_featShuffle > 0)
        numUnits_featShuffle <- ncol(groupMat_featShuffle)
        percAccessUnits_featShuffle <- numAccessUnits_featShuffle/numUnits_featShuffle*100
        numAccessUnits_cellShuffle <- matrixStats::rowSums2(groupMat_cellShuffle > 0)
        numUnits_cellShuffle <- ncol(groupMat_cellShuffle)
        percAccessUnits_cellShuffle <- numAccessUnits_cellShuffle/numUnits_cellShuffle*100
        
        o_featShuffle[idx,]$correlation <- as.numeric(corVals_featShuffle)
        o_cellShuffle[idx,]$correlation <- as.numeric(corVals_cellShuffle)
        
        o_featShuffle[idx, ]$Variability1 <- rowVars_featShuffle[o_featShuffle[idx, ]$idx1]
        o_featShuffle[idx, ]$Variability2 <- rowVars_featShuffle[o_featShuffle[idx, ]$idx2]
        o_cellShuffle[idx, ]$Variability1 <- rowVars_cellShuffle[o_cellShuffle[idx, ]$idx1]
        o_cellShuffle[idx, ]$Variability2 <- rowVars_cellShuffle[o_cellShuffle[idx, ]$idx2]
        
        
        o_featShuffle[idx, ]$PercAccess1 <- percAccessUnits_featShuffle[o_featShuffle[idx, ]$idx1]
        o_featShuffle[idx, ]$PercAccess2 <- percAccessUnits_featShuffle[o_featShuffle[idx, ]$idx2]
        o_featShuffle[idx, ]$PercAccessMean <- rowMeans(cbind(o_featShuffle[idx, ]$PercAccess1, o_featShuffle[idx, ]$PercAccess2))
        o_cellShuffle[idx, ]$PercAccess1 <- percAccessUnits_cellShuffle[o_cellShuffle[idx, ]$idx1]
        o_cellShuffle[idx, ]$PercAccess2 <- percAccessUnits_cellShuffle[o_cellShuffle[idx, ]$idx2]
        o_cellShuffle[idx, ]$PercAccessMean <- rowMeans(cbind(o_cellShuffle[idx, ]$PercAccess1, o_cellShuffle[idx, ]$PercAccess2))

        ArchR:::.logThis(groupMat, paste0("SubsetGroupMat-", x), logFile = logFile)
        ArchR:::.logThis(o[idx,], paste0("SubsetCoA-", x), logFile = logFile)
    }

    o_featShuffle$idx1 <- NULL
    o_featShuffle$idx2 <- NULL
    o_featShuffle <- o_featShuffle[!is.na(o_featShuffle$correlation),]

    o_cellShuffle$idx1 <- NULL
    o_cellShuffle$idx2 <- NULL
    o_cellShuffle <- o_cellShuffle[!is.na(o_cellShuffle$correlation),]

    mcols(featureSet) <- NULL

    o_featShuffle@metadata$featureSet <- featureSet             
    o_cellShuffle@metadata$featureSet <- featureSet

    ArchR:::.endLogging(logFile = logFile)

    list(BackgroundCutoff = max(quantile(abs(o_featShuffle$correlation), seq(0.01,1,by=0.01), na.rm=T)[99], 
                                quantile(abs(o_cellShuffle$correlation), seq(0.01,1,by=0.01), na.rm=T)[99]),
         BackgroundCoAccessibilityQuantiles = pmax(quantile(abs(o_featShuffle$correlation), seq(0.01,1,by=0.01), na.rm=T), 
                                                   quantile(abs(o_cellShuffle$correlation), seq(0.01,1,by=0.01), na.rm=T)),
         BackgroundCoAccessibility = list(featShuffle = o_featShuffle, cellShuffle = o_cellShuffle))
}
