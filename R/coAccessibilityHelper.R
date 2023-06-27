# This file contains helper functions for the getCoAccessibility and getBackgroundCoAccessibility methods. 
# The separation into smaller functions allows for less repetition, cleaner code and easier further development.
# Majority of the functions mirror the existing functionality of addCoAccessibility in ArchR. 
# Some functions are adapted to include the additional functionality in RWireX.

`%notin%` <- Negate(`%in%`)

#' Input parameter check.
#' 
#' @description The function takes the list of possible parameters and checks the input type.
#' @keywords internal
#' @export
.checkInputParameters <- function(parameterList){
  validInput = list(ArchRProj = c("ArchRProj"),
                    chromsToUse = c("character", "null"),
                    reducedDims= c("character"),
                    dimsToUse = c("numeric", "null"),
                    scaleDims = c("boolean", "null"),
                    cellsToUse = c("character", "null"),
                    AggregationMethod = c("character"),
                    numCellsPerAggregate = c("integer"),
                    numAggregates = c("integer"),
                    useMatrix = c("character"),
                    binaryMatrix = c("boolean"),
                    overlapCutoff = c("numeric"),
                    maxDist = c("integer"),
                    chromosomewise = c("boolean"),
                    scaleTo = c("numeric"),
                    log2Norm = c("boolean"),
                    seed = c("integer"),
                    threads = c("integer"),
                    numPermutations = c("integer"),
                    verbose = c("boolean"),
                    returnLoops = c("boolean"),
                    corCutOff = c("numeric"),
                    logFile = c("character"))
  for(parameter in names(parameterList)){
    if (parameter %notin% names(validInput)){
      stop("No valid input defined for ", parameter,"!")
    }
    ArchR:::.validInput(input = parameterList[[parameter]], name = parameter, valid = validInput[[parameter]])
  }
}

#' Cell number parameter check.
#' 
#' @description The function checks whether there are enough cells and adapts parameters if needed.
#' @keywords internal
#' @export
.checkNumAggregatesAndCellsPerAggregate <- function(AggregationMethod, numAggregates, numCellsPerAggregate, cells_number){
  if (AggregationMethod == "single_cell_resolution"){
    numCellsPerAggregate <- 1
    numAggregates <- cells_number
    return(list(numAggregates, numCellsPerAggregate))
  }
  else if (AggregationMethod == "unique"){
    if (cells_number < numAggregates*numCellsPerAggregate){
      stop("Not enough cells! Product of numAggregates and  numCellsPerAggregate must be less or equal to the number of cells!")
    }
  }

  return(list(numAggregates, numCellsPerAggregate))
}



#' Feature Set getter
#' 
#' @description The function returns feature set based on the specified matrix. 
#' @keywords internal
#' @export
.getSet <- function(ArchRProj, useMatrix, binaryMatrix){
  if (useMatrix == "PeakMatrix"){
    set <- getPeakSet(ArchRProj)
  } else if (useMatrix == "TileMatrix"){
    tileSet <- getMatrixFromProject(ArchRProj, useMatrix = "TileMatrix", binarize = binaryMatrix)@elementMetadata
    set <- GRanges(seqnames = tileSet$seqnames, 
                       ranges = IRanges(start = tileSet$start, width = tileSet$start[2]-tileSet$start[1]))
    set$idx <- tileSet$idx
    set$id <- 1:length(set)
  }
  else{
    stop("This analysis only works with Peak or Tile Matrix sofar.")
  }
  return(set)
}

#' Reduced Dimensions getter.
#' 
#' @description The function gets the reduced dimensions and, if specified, filters the cells.
#' @keywords internal
#' @export
.getFilteredReducedDimensions <- function(ArchRProj, reducedDims, corCutOff, dimsToUse, cellsToUse){
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if (!is.null(cellsToUse)) {
    rD <- rD[cellsToUse, , drop = FALSE]
  }
  return(rD)
}


#' Cell seeds selection.
#' 
#' @description The function creates the list of cells on which the co-accessibility analysis is performed.
#' Depends on the AggregationMethod parameter.
#' @keywords internal
#' @export
.selectCellSeedsForAggregation <-function(ArchRProj, reducedDimensions, AggregationMethod, numPermutations, numCellsPerAggregate, numAggregates, cellsToUse){
  ### Select "cell seeds" for aggregation
  if(AggregationMethod == "unique"){
    ### Different groups of cell are chosen at random "numPermutations" times and from that the cell group with most distance among cells is selected 
    ### for low dimensional embedding.
    # From: https://stackoverflow.com/questions/22152482/choose-n-most-distant-points-in-r
    bestavg <- 0
    bestSet <- NA
    for (i in 1:numPermutations){
      subset <- reducedDimensions[sample(1:nrow(reducedDimensions),numAggregates),]
      avg <- mean(dist(subset))
      if (avg > bestavg) {
        bestavg <- avg
        bestSet <- subset
      }
    }
    idx <- match(rownames(bestSet), rownames(reducedDimensions))
  } else if (AggregationMethod == "ArchR_default"){
    idx <- sample(seq_len(nrow(reducedDimensions)), numAggregates, replace = !nrow(reducedDimensions) >= numAggregates)
  } else if (AggregationMethod == "single_cell_resolution"){
    idx <- numCellsPerAggregate:numAggregates
  }
  return(idx)
}

#' Closest cells to seeds selection.
#' 
#' @description The function performs KNN to find the closest cells to seeds.
#' @keywords internal
#' @export
.selectClosestCellsOfCellSeeds <- function(ArchRProj, reducedDimensions, idx, AggregationMethod, numAggregates, numCellsPerAggregate){
  expand = 10
  if(AggregationMethod == "unique"){
    cell_number = nrow(reducedDimensions)
    if (numCellsPerAggregate*expand > cell_number-numAggregates){
      error("Not enough cells to consider for possible unique aggregates. Consider lowering numCellsPerAggregate or choosing another algorithm.")
    }
    used_cells = unique(c(idx))
    seed_cells_rD <- reducedDimensions[idx, ] %>% as.matrix(.)
    closest_cells <- ArchR:::.computeKNN(data = reducedDimensions[-used_cells,], query = seed_cells_rD, k = numCellsPerAggregate*expand)
    
    knnObj = matrix(, nrow = 0, ncol = numCellsPerAggregate)
    ### Select every cell only once (including in the aggregates around seed cells)
    for (i in 1:numAggregates){
      names <- rownames(reducedDimensions[-idx,])[as.integer(closest_cells[i,])]
      idx_closest_cells = match(names, rownames(reducedDimensions))
      neighbours = idx_closest_cells[idx_closest_cells %notin% used_cells]
      if  (length(neighbours)>=numCellsPerAggregate){
        seed_closest_cells = neighbours[1:numCellsPerAggregate-1]
        
        knnObj <- rbind(knnObj, c(idx[i], seed_closest_cells))
        used_cells = unique(c(idx, knnObj %>% as.integer()))
      }
    }
  } else if (AggregationMethod %in% c("ArchR_default")){
    knnObj <- ArchR:::.computeKNN(data = reducedDimensions, query = reducedDimensions[idx, ], k = numCellsPerAggregate)
  } else if (AggregationMethod == "single_cell_resolution"){
    knnObj <- matrix(idx, ncol = 1)
  }
  return(knnObj)
}

#' Filter KNN object.
#' 
#' @description The function calculates the overlap between KNNs of each seed cell and keeps those with more than overlapCutoff unique cells in KNN.
#' It does not apply to single cell level mode of aggregation, because we want to keep all cells. 
#' @keywords internal
#' @export
.filterKNN <- function(knnObj, AggregationMethod, overlapCutoff, numCellsPerAggregate, numAggregates){
  if (AggregationMethod == "ArchR_default"){
    #Determine Overlap
    keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * numCellsPerAggregate))
    #Filter out
    knnObj <- knnObj[keepKnn == 0, ]
  }
  else if (AggregationMethod == "single_cell_resolution"){
    knnObj <- matrix(knnObj, nrow = numAggregates)
  }
  else if(AggregationMethod == "unique"){
    return(knnObj)
  }
  else{
    stop("The aggregation method is not supported!")
  }
  return(knnObj)
}

#' Pairwise list of features.
#' 
#' @description The function creates the pairwise list of features for the co-accessibility analysis.
#' @keywords internal
#' @export
.createPairwiseThingsToTest <- function(peakSet, maxDist){
  #Create Ranges
  peakSummits <- resize(peakSet, 1, "center")
  peakWindows <- resize(peakSummits, 2*maxDist + 1, "center")
  
  #Create Pairwise Things to Test
  o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
  o <- o[o[, 1] != o[, 2], ]
  o$seqnames <- seqnames(peakSet)[o[, 1]]
  o$idx1 <- peakSet$idx[o[, 1]]
  o$idx2 <- peakSet$idx[o[, 2]]
  o$correlation <- -999.999
  o$Variability1 <- 0.000
  o$Variability2 <- 0.000
  o$PercAccess1 <- 0.000
  o$PercAccess2 <- 0.000
  o$PercAccessMean <- 0.000
  o$peaksetLookUpIndex1 <- o[, 1]
  o$peaksetLookUpIndex2 <- o[, 2]
  return(o)
}
  
#' Adding matrix metadata.
#' 
#' @description The function adds closest cells to seeds and normalized matrix as metadata to the analysis result.
#' @keywords internal
#' @export  
.addMetadataForAggregates <- function(ArchRProj, o, peakSet, knnObj, useMatrix, gS, log2Norm, scaleTo){
  ### Additional metadata for aggregates
  featureDF <- mcols(peakSet)
  featureDF$seqnames <- peakSet@seqnames
  
  #Group Matrix
  groupMat <- ArchR:::.getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF, 
    groupList = knnObj, 
    useMatrix = useMatrix,
    verbose = FALSE
  )
  
  #Scale
  groupMat <- t(t(groupMat) / gS) * scaleTo
  
  if (log2Norm) {
    groupMat <- log2(groupMat + 1)
  }
  
  o@metadata$Aggregates <- knnObj
  o@metadata$AggregatePeakMatrix <- groupMat
  return(o)
}

#' Get and normalize Group Matrix.
#' 
#' @description The function extracts the matrix for the given features and scales/normalizes it. 
#' @keywords internal
#' @export
.createGroupMatrix <- function(ArchRProj, peakSet, knnObj, useMatrix, gS, log2Norm, chr, scaleTo){
  #Features
  featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == chr), ]
  featureDF$seqnames <- chr
  
  #Group Matrix
  groupMat <- ArchR:::.getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF, 
    groupList = knnObj, 
    useMatrix = useMatrix,
    verbose = FALSE
  )
  
  #Scale
  groupMat <- t(t(groupMat)/gS) * scaleTo
  
  if (log2Norm) {
    groupMat <- log2(groupMat + 1)
  }
  return(groupMat)
}

#' Get and normalize Group Matrix.
#' 
#' @description The function extracts the matrix for the given features and scales/normalizes it. 
#' @keywords internal
#' @export
.createCustomGroupMatrix <- function(ArchRProj, customMatrix, peakSet, knnObj, useMatrix, gS, log2Norm, chr, scaleTo){
  #Features
  featureDF_index = BiocGenerics::which(seqnames(peakSet) == chr)
  #Group Matrix
  customMatrixFiltered = customMatrix[featureDF_index, ]
  customMatrixFiltered = as.matrix(customMatrixFiltered)

  mat <- matrix(0, nrow = length(featureDF_index), ncol = length(knnObj))
  colnames(mat) <- names(knnObj)
  for(z in seq_along(knnObj)){
    cellsGroupz <- knnObj[[z]]
    idx <- BiocGenerics::which(colnames(customMatrixFiltered) %in% cellsGroupz)
    mat[,z] <- Matrix::rowSums(as.matrix(customMatrixFiltered[,idx,drop=FALSE]))
      
  }
  #Scale
  groupMat <- t(t(mat)/gS) * scaleTo
  
  if (log2Norm) {
    groupMat <- log2(groupMat + 1)
  }
  return(groupMat)
  
    
}

#' Get ColSums for custom matrix.
#' 
#' @description The function substitutes .getColSums from ArchR. 
#' @keywords internal
#' @export
.getColSumsCustomMatrix <- function(customMatrix){
  cS = colSums2(customMatrix)
  names(cS) = colnames(customMatrix)
  return(cS)
}

#' Create Loops between Peaks.
#' 
#' @description The function creates loops between pairs of peaks from peak indices.
#' @keywords internal
#' @export
.createLoopsSameChromosome <- function(coA){
  peakSummits <- resize(metadata(coA)$featureSet, 1, "center")
  summitTiles <- start(peakSummits)
  
  loops <- ArchR:::.constructGR(
    seqnames = seqnames(peakSummits[coA[,1]]),
    start = summitTiles[coA[,1]],
    end = summitTiles[coA[,2]]
  )
  
  mcols(loops) <- coA

  loops <- loops[order(mcols(loops)$correlation, decreasing=TRUE)]
  loops <- unique(loops)
  loops <- loops[width(loops) > 0]
  loops <- sort(sortSeqlevels(loops))
  
  loops <- SimpleList(CoAccessibility = loops)
  metadata(loops) = metadata(coA)
  return(loops)
}
