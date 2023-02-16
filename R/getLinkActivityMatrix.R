#' Get Link Activity Matix
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param coAccessibility The co-accessibility data produced by RWireX method getCoAccessibility.
#' @param useMatrix Same matrix type that was used for the accessibility.
#' @param correlationCoefficient multiply links by correlation.
#' @param binarize if matrix needs to be binarized. Have to be the same as the matrix used for coAccessibility calculation.
#' 
#' @export


getLinkActivityMatrix <- function (ArchRProj, coAccessibility, useMatrix="PeakMatrix", binarize=TRUE, correlationCoefficient=FALSE){
  matr = getMatrixFromProject(ArchRProj, useMatrix=useMatrix, binarize = binarize)
  if (!(identical(matr@rowRanges@ranges@start, coAccessibility@metadata$featureSet@ranges@start)&&identical(matr@rowRanges@seqnames, coAccessibility@metadata$featureSet@seqnames))){
    stop("The order of peaks in matrix and coAccessibility Feature Set is not the same. The lookup method is not implemented yet.")
  }
  matrix_data = matr@assays@data$PeakMatrix
  link_matrix = matrix_data[coAccessibility$peaksetLookUpIndex1,]+matrix_data[coAccessibility$peaksetLookUpIndex2,]
  if (correlationCoefficient){
    link_matrix = link_matrix * coAccessibility$correlation
  }
  return(link_matrix)
}