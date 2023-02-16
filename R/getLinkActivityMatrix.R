#' Get Link Activity Matix
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param coAccessibility The co-accessibility data produced by RWireX method getCoAccessibility.
#' @param useMatrix Same matrix type that was used for the accessibility.
#' @param peakSubset If specified, then the matrix is subsetted by the links with the first peak in the peakSubset.
#' @param rangeForSubset If specified, then the matrix is subsetted by the links where the first peak is in this range.
#' @param correlationCoefficient multiply links by correlation.
#' @param binarize if matrix needs to be binarized. Have to be the same as the matrix used for coAccessibility calculation.
#' 
#' @export


getLinkActivityMatrix <- function (ArchRProj, coAccessibility, useMatrix="PeakMatrix", peakSubset=NULL, rangeForSubset=NULL, binarize=TRUE, correlationCoefficient=FALSE){
  matr = getMatrixFromProject(ArchRProj, useMatrix=useMatrix, binarize = binarize)
  if (!(identical(matr@rowRanges@ranges@start, coAccessibility@metadata$featureSet@ranges@start)&&identical(as.character(matr@rowRanges@seqnames), as.character(coAccessibility@metadata$featureSet@seqnames)))){
    stop("The order of peaks in matrix and coAccessibility Feature Set is not the same. The lookup method is not implemented yet.")
  }
  matrix_data = matr@assays@data$PeakMatrix
  if (!is.null(peakSubset)&&!is.null(rangeForSubset)){
    stop("Please choose only one subsetting option: either provide a peak list or a range.")
  }
  if (!is.null(peakSubset)){
    index_of_subset = match(peakSubset, matr@rowRanges)
    coAccessibility = coAccessibility[coAccessibility$peaksetLookUpIndex1 %in% index_of_subset, ]
  }
  if (!is.null(rangeForSubset)){
    index_of_subset = which(!is.na(findOverlaps(matr@rowRanges, rangeForSubset, select="first")))
    coAccessibility = coAccessibility[coAccessibility$peaksetLookUpIndex1 %in% index_of_subset, ]
  }
  link_matrix = matrix_data[coAccessibility$peaksetLookUpIndex1,]+matrix_data[coAccessibility$peaksetLookUpIndex2,]
  if (correlationCoefficient){
    link_matrix = link_matrix * coAccessibility$correlation
  }
  return(link_matrix)
}