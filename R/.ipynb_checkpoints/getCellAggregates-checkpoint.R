#' getCellAggregates
#' 
#' This function returns cell aggregates and Peak matrix of cell aggregates for further inspection.
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
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
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


getCellAggregates <- function (
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
    ###
    overlapCutoff = 0.8, 
    scaleTo = 10^4,
    log2Norm = TRUE, 
    seed = 1, 
    threads = getArchRThreads(), 
    verbose = TRUE, 
    logFile = createLogFile("getCellAggregates")
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
    ###
    ArchR:::.validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
    ArchR:::.validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
    ArchR:::.validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
    ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
    ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
    ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
    
    tstart <- Sys.time()
    ArchR:::.startLogging(logFile = logFile)
    ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), "getCellAggregates Input-Parameters", logFile = logFile)
    
    set.seed(seed)
    
    #Get Peak Set
    peakSet <- getPeakSet(ArchRProj)
    
    #Get Reduced Dims
    rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
    if (!is.null(cellsToUse)) {
        rD <- rD[cellsToUse, , drop = FALSE]
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
    }
    ###

    #KNN Matrix
    ArchR:::.logDiffTime(main = "Computing KNN", t1 = tstart, verbose = verbose, logFile = logFile)
    
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
    } else if (AggregationMethod == "duplicated"){
        knnObj <- ArchR:::.computeKNN(data = rD, query = rD[idx, ], k = numCellsPerAggregate)
    }
    ###

    #Determine Overlap
    ArchR:::.logDiffTime(main = "Identifying Non-Overlapping KNN pairs", t1 = tstart, verbose = verbose, logFile = logFile)
    keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * numCellsPerAggregate))

    #Keep Above Cutoff
    knnObj <- knnObj[keepKnn==0,]
    ArchR:::.logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1 = tstart, verbose = verbose, logFile = logFile)

    #Convert To Names List
    knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
        rownames(rD)[knnObj[x, ]]
    }) %>% SimpleList

    #Check Chromosomes
    chri <- gtools::mixedsort(ArchR:::.availableChr(getArrowFiles(ArchRProj), subGroup = "PeakMatrix"))
    chrj <- gtools::mixedsort(unique(paste0(seqnames(getPeakSet(ArchRProj)))))
    stopifnot(identical(chri, chrj))
    
    #Peak Matrix ColSums
    cS <- ArchR:::.getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = "PeakMatrix")
    gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm = TRUE)))
                        
    ### ADAPTED IL ###
    ### Add pseudo-count to gS for non-accessible cell aggregates in all regions of interest
    gS <- gS + 1
    ###

    #Features
    ### ADAPTED IL ###
    #featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == chri[x]), ]
    #featureDF$seqnames <- chri[x]
    featureDF <- mcols(peakSet)
    featureDF$seqnames <- peakSet@seqnames
    ###

    #Group Matrix
    groupMat <- ArchR:::.getGroupMatrix(
        ArrowFiles = getArrowFiles(ArchRProj), 
        featureDF = featureDF, 
        groupList = knnObj, 
        useMatrix = "PeakMatrix",
        verbose = FALSE
    )

    #Scale
    groupMat <- t(t(groupMat) / gS) * scaleTo
                        
    if (log2Norm) {
        groupMat <- log2(groupMat + 1)
    }

    return(list(Aggregates = knnObj,
                AggregatePeakMatrix = groupMat))
}
                        
### Adapted from addCoAccessibility function from https://github.com/GreenleafLab/ArchR/blob/968e4421ce7187a8ac7ea1cf6077412126876d5f/R/IntegrativeAnalysis.R on 22/03/2022