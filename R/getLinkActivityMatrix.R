#' Get Link Activity Matix
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param CoAccessibility The co-accessibility data produced by RWireX method getCoAccessibility.
#' @param useMatrix Same matrix type that was used for the accessibility.
#' @param peakSubset If specified, then the matrix is subsetted by the links with the first peak in the peakSubset.
#' @param regionSubset If specified, then the matrix is subsetted by the links where the first peak is in this range.
#' @param correlationCoefficient multiply links by correlation.
#' @param binarize if matrix needs to be binarized. Have to be the same as the matrix used for CoAccessibility calculation.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging RWireX output.
#' 
#' @export


getLinkActivityMatrix <- function (ArchRProj, CoAccessibility, useMatrix="PeakMatrix", peakSubset=NULL, regionSubset=NULL, binarize=TRUE, correlationCoefficient=FALSE, logFile = createLogFile("getLinkActivityMatrix"), verbose=TRUE){
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = CoAccessibility, name = "CoAccessibility", valid = c("list"))
  ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  ArchR:::.validInput(input = peakSubset, name = "peakSubset", valid = c("list", "null"))
  ArchR:::.validInput(input = regionSubset, name = "regionSubset", valid = c("GRanges", "GRangesList", " IntegerRanges", "null"))
  ArchR:::.validInput(input = binarize, name = "binarize", valid = c("logical"))
  ArchR:::.validInput(input = correlationCoefficient, name = "correlationCoefficient", valid = c("logical"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("logical"))

  if (class(CoAccessibility) == "DFrame"){
    CoAccessibility <- .createLoopsSameChromosome(CoAccessibility)
  }
  
  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), "getLinkActivityMatrix Input-Parameters", logFile = logFile)
  
  matr = getMatrixFromProject(ArchRProj, useMatrix=useMatrix, binarize = binarize, logFile = logFile)
  if (!(identical(matr@rowRanges@ranges@start, CoAccessibility@metadata$featureSet@ranges@start)&&identical(as.character(matr@rowRanges@seqnames), as.character(CoAccessibility@metadata$featureSet@seqnames)))){
    stop("The order of peaks in matrix and CoAccessibility Feature Set is not the same. The lookup method is not implemented yet.")
  }
  matrix_data = matr@assays@data$PeakMatrix
  if (!is.null(peakSubset)&&!is.null(regionSubset)){
    stop("Please choose only one subsetting option: either provide a peak list or a range.")
  }
  
  CoAccessibility <- CoAccessibility$CoAccessibility
  if (!is.null(peakSubset)){
    index_of_subset <- S4Vectors::match(peakSubset, matr@rowRanges)
    CoAccessibility <- CoAccessibility[CoAccessibility$peaksetLookUpIndex1 %in% index_of_subset, ]
  }
  if (!is.null(regionSubset)){
    index_of_subset <- which(!is.na(findOverlaps(matr@rowRanges, regionSubset, select="first")))
    CoAccessibility <- CoAccessibility[CoAccessibility$peaksetLookUpIndex1 %in% index_of_subset, ]
  }
  
  ArchR:::.logDiffTime(main = "Generating Link Matrix", t1 = tstart, verbose = verbose, logFile = logFile)
  link_matrix <- matrix_data[CoAccessibility$peaksetLookUpIndex1,] * matrix_data[CoAccessibility$peaksetLookUpIndex2,]
  link_matrix@Dimnames[[1]] <- paste0(CoAccessibility$seqnames, ":", CoAccessibility$idx1, "-",
                                      CoAccessibility$seqnames, ":", CoAccessibility$idx2)
  if (correlationCoefficient){
    link_matrix <- link_matrix * CoAccessibility$correlation
  }
  
  ArchR:::.endLogging(logFile = logFile)
  return(link_matrix)
}
