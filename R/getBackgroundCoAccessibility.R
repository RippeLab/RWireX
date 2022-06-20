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
#' @param AggregationMethod Choose cell aggregation method from "duplicated" (each cell in multiple aggregates), "unique" (each cell in only one aggregate), or "none" (no cell aggregation, single-cell level).
#' @param numCellsPerAggregate The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param numAggregates The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param useMatrix The name of the data matrix to use for calculation of BackgroundCoAccessibility. Options include "PeakMatrix" or non-binarized "TileMatrix" (binarize = FALSE in addTileMatrix).
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions from
#' the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended to keep track
#' of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
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
    ### ADAPTED IL ###
    #k = 100, 
    #knnIteration = 500,
    AggregationMethod = "duplicated", 
    numCellsPerAggregate = 100, 
    numAggregates = 500,
    useMatrix = "PeakMatrix",
    ###
    overlapCutoff = 0.8, 
    maxDist = 100000,
    scaleTo = 10^4,
    log2Norm = TRUE,
    seed = 1, 
    threads = getArchRThreads(),
    verbose = TRUE,
    logFile = createLogFile("getBackgroundCoAccessibility")
    ){

    ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    ArchR:::.validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
    ArchR:::.validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
    ArchR:::.validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
    ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
    ArchR:::.validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
    ### ADAPTED IL ###
    #.validInput(input = k, name = "k", valid = c("integer"))
    #.validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
    ArchR:::.validInput(input = AggregationMethod, name = "AggregationMethod", valid = c("character"))
    ArchR:::.validInput(input = numCellsPerAggregate, name = "numCellsPerAggregate", valid = c("integer"))
    ArchR:::.validInput(input = numAggregates, name = "numAggregates", valid = c("integer"))
    ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
    ###
    ArchR:::.validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
    ArchR:::.validInput(input = maxDist, name = "maxDist", valid = c("integer"))
    ArchR:::.validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
    ArchR:::.validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
    ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
    ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
    ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))

    tstart <- Sys.time()
    ArchR:::.startLogging(logFile = logFile)
    ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "getBackgroundCoAccessibility Input-Parameters", logFile = logFile)

    set.seed(seed)

    #Get Peak Set
    ### ADAPTED IL ###
    if (useMatrix == "PeakMatrix"){
        peakSet <- getPeakSet(ArchRProj)
    } else if (useMatrix == "TileMatrix"){
        tileSet <- getMatrixFromProject(ArchRProj, useMatrix = "TileMatrix")@elementMetadata
        peakSet <- GRanges(seqnames = tileSet$seqnames, 
                           ranges = IRanges(start = tileSet$start, width = tileSet$start[2]-tileSet$start[1]))
        peakSet$idx <- tileSet$idx
        peakSet$id <- 1:length(peakSet)
    }
    ###

    #Get Reduced Dims
    rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
    if(!is.null(cellsToUse)){
        rD <- rD[cellsToUse, ,drop=FALSE]
    }

    #Subsample
    ### ADAPTED IL ###
    #idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
    ### Select "cell seeds" for aggregation
    if(AggregationMethod == "unique"){
      ### Select most distant (randomly distributed) cells in low dimensional embedding
      # From: https://stackoverflow.com/questions/22152482/choose-n-most-distant-points-in-r
      bestavg <- 0
      bestSet <- NA
      for (i in 1:nrow(rD)){
          subset <- rD[sample(1:nrow(rD),numAggregates),]
          avg <- mean(dist(subset))
          if (avg > bestavg) {
              bestavg <- avg
              bestSet <- subset
          }
      }
      idx <- match(rownames(bestSet), rownames(rD))
    } else if (AggregationMethod == "duplicated"){
      idx <- sample(seq_len(nrow(rD)), numAggregates, replace = !nrow(rD) >= numAggregates)
    } else if (AggregationMethod == "none"){
      numCellsPerAggregate <- 1
      if (is.null(cellsToUse)) {
          numAggregates <- nrow(ArchRProj@cellColData)
      } else {
          numAggregates <- length(cellsToUse)
      }
      idx <- numCellsPerAggregate:numAggregates
    }
    ###

    #KNN Matrix
    ArchR:::.logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)

    ### ADAPTED IL ###
    # knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)
    ### Select closest cells of cell seeds
    if(AggregationMethod == "unique"){
      ### Select every cell only once
      query1 <- rD[idx[1], ] %>% as.matrix(.) %>% t(.)
      rownames(query1) <- rownames(rD)[1]
      y <- ArchR:::.computeKNN(data = rD[-idx,], query = query1, k = numCellsPerAggregate-1)
      names <- rownames(rD[-idx,])[as.integer(y)]
      knnObj <- c(idx[1], match(names, rownames(rD))) %>% as.matrix(.) %>% t(.)

      for (i in 2:numAggregates){
          query1 <- rD[idx[i], ] %>% as.matrix(.) %>% t(.)
          rownames(query1) <- rownames(rD)[i]
          y <- ArchR:::.computeKNN(data = rD[-unique(c(idx, knnObj %>% as.integer())),], query = query1, k = numCellsPerAggregate-1)
          names <- rownames(rD[-unique(c(idx, knnObj %>% as.integer())),])[as.integer(y)]
          knnObj <- rbind(knnObj, c(idx[i], match(names, rownames(rD))))
      }
    } else if (AggregationMethod %in% c("duplicated", "none")){
      knnObj <- ArchR:::.computeKNN(data = rD, query = rD[idx, ], k = numCellsPerAggregate)
    }
    ###

    #Determine Overlap
    ArchR:::.logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
    keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * numCellsPerAggregate))

    #Keep Above Cutoff
    knnObj <- knnObj[keepKnn==0,]
    ArchR:::.logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)

    ### ADAPTED IL ###
    ### Convert integer to matrix:array
    if (AggregationMethod == "none"){
      knnObj <- matrix(knnObj, nrow = numAggregates)
    }
    ###

    #Convert To Names List
    knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
    }) %>% SimpleList

    #Check Chromosomes
    ### ADAPTED IL ###
    #chri <- gtools::mixedsort(ArchR:::.availableChr(getArrowFiles(ArchRProj), subGroup = "PeakMatrix"))
    #chrj <- gtools::mixedsort(unique(paste0(seqnames(getPeakSet(ArchRProj)))))
    chri <- gtools::mixedsort(ArchR:::.availableChr(getArrowFiles(ArchRProj), subGroup = useMatrix))
    chrj <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
    ###
    stopifnot(identical(chri,chrj))

    #Create Ranges
    peakSummits <- resize(peakSet, 1, "center")
    peakWindows <- resize(peakSummits, maxDist, "center")

    #Create Pairwise Things to Test
    o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
    o <- o[o[,1] != o[,2],]
    o$seqnames <- seqnames(peakSet)[o[,1]]
    o$idx1 <- peakSet$idx[o[,1]]
    o$idx2 <- peakSet$idx[o[,2]]
    o$correlation <- -999.999
    o$Variability1 <- 0.000
    o$Variability2 <- 0.000
    ### ADAPTED IL ###
    o$PercAccess1 <- 0.000
    o$PercAccess2 <- 0.000
    
    o_featShuffle <- o
    o_cellShuffle <- o
    ###

    #Peak Matrix ColSums
    ### ADAPTED IL ###
    #cS <- ArchR:::.getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = "PeakMatrix")
    cS <- ArchR:::.getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = useMatrix)
    ###
    gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))

    ### ADAPTED IL ###
    ### Add pseudo-count to gS for non-accessible cell aggregates in all regions of interest
    gS <- gS + 1
    ###

    for(x in seq_along(chri)){

        ArchR:::.logDiffTime(sprintf("Computing Background Co-Accessibility %s (%s of %s)", chri[x], x, length(chri)), t1=tstart, verbose=verbose, logFile=logFile)

        #Features
        featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == chri[x]),]
        featureDF$seqnames <- chri[x]

        #Group Matrix
        groupMat <- ArchR:::.getGroupMatrix(
            ArrowFiles = getArrowFiles(ArchRProj), 
            featureDF = featureDF, 
            groupList = knnObj, 
            ### ADAPTED IL ###
            #useMatrix = "PeakMatrix", 
            useMatrix = useMatrix,
            ###
            threads = threads,
            verbose = FALSE
        )

        #Scale
        groupMat <- t(t(groupMat) / gS) * scaleTo
        
        if (log2Norm) {
            groupMat <- log2(groupMat + 1)
        }

        ### ADAPTED IL ###
        ###Shuffle features for background correlations
        groupMat_featShuffle <- apply(groupMat, 2, function(x){sample(x)}) ### shuffle features
        groupMat_cellShuffle <- t(apply(groupMat, 1, function(x){sample(x)})) ### shuffle cells
        ###

        #Correlations
        idx <- BiocGenerics::which(o$seqnames==chri[x])
        ### ADAPTED IL ###
        #corVals <- ArchR:::rowCorCpp(idxX = o[idx, ]$idx1, idxY = o[idx, ]$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))
        #ArchR:::.logThis(head(corVals), paste0("SubsetCorVals-", x), logFile = logFile)
        corVals_featShuffle <- ArchR:::rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat_featShuffle), Y = as.matrix(groupMat_featShuffle))
        corVals_cellShuffle <- ArchR:::rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat_cellShuffle), Y = as.matrix(groupMat_cellShuffle))
        ArchR:::.logThis(head(corVals_featShuffle), paste0("SubsetCorVals-", x), logFile = logFile)
        ArchR:::.logThis(head(corVals_cellShuffle), paste0("SubsetCorVals-", x), logFile = logFile)
        ###

        ### ADAPTED IL ###
        #rowVars <- as.numeric(matrixStats::rowVars(groupMat))
        rowVars_featShuffle <- as.numeric(matrixStats::rowVars(groupMat_featShuffle))
        rowVars_cellShuffle <- as.numeric(matrixStats::rowVars(groupMat_cellShuffle))
        
        numAccessUnits_featShuffle <- matrixStats::rowSums2(groupMat_featShuffle > 0)
        numUnits_featShuffle <- ncol(groupMat_featShuffle)
        percAccessUnits_featShuffle <- numAccessUnits_featShuffle/numUnits_featShuffle*100
        numAccessUnits_cellShuffle <- matrixStats::rowSums2(groupMat_cellShuffle > 0)
        numUnits_cellShuffle <- ncol(groupMat_cellShuffle)
        percAccessUnits_cellShuffle <- numAccessUnits_cellShuffle/numUnits_cellShuffle*100
        
        #o[idx, ]$correlation <- as.numeric(corVals)
        o_featShuffle[idx,]$correlation <- as.numeric(corVals_featShuffle)
        o_cellShuffle[idx,]$correlation <- as.numeric(corVals_cellShuffle)
        
        #o[idx, ]$Variability1 <- rowVars[o[idx, ]$idx1]
        #o[idx, ]$Variability2 <- rowVars[o[idx, ]$idx2]
        o_featShuffle[idx, ]$Variability1 <- rowVars_featShuffle[o_featShuffle[idx, ]$idx1]
        o_featShuffle[idx, ]$Variability2 <- rowVars_featShuffle[o_featShuffle[idx, ]$idx2]
        o_cellShuffle[idx, ]$Variability1 <- rowVars_cellShuffle[o_cellShuffle[idx, ]$idx1]
        o_cellShuffle[idx, ]$Variability2 <- rowVars_cellShuffle[o_cellShuffle[idx, ]$idx2]
        
        o_featShuffle[idx, ]$PercAccess1 <- percAccessUnits_featShuffle[o_featShuffle[idx, ]$idx1]
        o_featShuffle[idx, ]$PercAccess2 <- percAccessUnits_featShuffle[o_featShuffle[idx, ]$idx2]
        o_cellShuffle[idx, ]$PercAccess1 <- percAccessUnits_cellShuffle[o_cellShuffle[idx, ]$idx1]
        o_cellShuffle[idx, ]$PercAccess2 <- percAccessUnits_cellShuffle[o_cellShuffle[idx, ]$idx2]
        ###

        ArchR:::.logThis(groupMat, paste0("SubsetGroupMat-", x), logFile = logFile)
        ArchR:::.logThis(o[idx,], paste0("SubsetCoA-", x), logFile = logFile)
    }

    ### ADAPTED IL ###
    #o$idx1 <- NULL
    #o$idx2 <- NULL
    #o <- o[!is.na(o$correlation), ]

    o_featShuffle$idx1 <- NULL
    o_featShuffle$idx2 <- NULL
    o_featShuffle <- o_featShuffle[!is.na(o_featShuffle$correlation),]

    o_cellShuffle$idx1 <- NULL
    o_cellShuffle$idx2 <- NULL
    o_cellShuffle <- o_cellShuffle[!is.na(o_cellShuffle$correlation),]

    #o$TStat <- (o$correlation / sqrt((pmax(1-o$correlation^2, 0.00000000000000001, na.rm = TRUE))/(length(knnObj)-2))) #T-statistic P-value
    #o$Pval <- 2 * pt(-abs(o$TStat), length(knnObj) - 2)
    #o$FDR <- p.adjust(o$Pval, method = "fdr")
    #o$VarQuantile1 <- ArchR:::.getQuantiles(o$Variability1)
    #o$VarQuantile2 <- ArchR:::.getQuantiles(o$Variability2)
    ###

    mcols(peakSet) <- NULL

    ### ADAPTED IL ###
    #o@metadata$peakSet <- peakSet
    o_featShuffle@metadata$peakSet <- peakSet             
    o_cellShuffle@metadata$peakSet <- peakSet

    #metadata(ArchRProj@peakSet)$CoAccessibility <- o
    ###

    ArchR:::.endLogging(logFile = logFile)

    ### ADAPTED IL ###
    #ArchRProj
    list(BackgroundCutoff = max(quantile(abs(o_featShuffle$correlation), seq(0.01,1,by=0.01), na.rm=T)[99], 
                                quantile(abs(o_cellShuffle$correlation), seq(0.01,1,by=0.01), na.rm=T)[99]),
         BackgroundCoAccessibilityQuantiles = pmax(quantile(abs(o_featShuffle$correlation), seq(0.01,1,by=0.01), na.rm=T), 
                                                   quantile(abs(o_cellShuffle$correlation), seq(0.01,1,by=0.01), na.rm=T)),
         BackgroundCoAccessibility = list(featShuffle = o_featShuffle, cellShuffle = o_cellShuffle))
    ##                   
}

                      
### Adapted from addCoAccessibility function from https://github.com/GreenleafLab/ArchR/blob/968e4421ce7187a8ac7ea1cf6077412126876d5f/R/IntegrativeAnalysis.R on 22/03/2022