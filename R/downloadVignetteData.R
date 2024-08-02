#' Download Data for RWireX Vignette
#' 
#' This function will download data for the vignette and store it in the specified location.
#' @param directory_to_store directory where ArchRProject will be stored
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging RWireX output.
#' @export
downloadVignetteData <- function(directory_to_store = getwd(),
                                 logFile = createLogFile("downloadVignetteData"), 
                                 verbose = TRUE){
  
  ArchR:::.validInput(input = directory_to_store, name = "directory_to_store", valid = c("character"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("logical"))
  
  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), "downloadVignetteData Input-Parameters", logFile = logFile)
  #Make Sure URL doesnt timeout
  oldTimeout <- getOption('timeout')
  options(timeout=150000)
  
  ArchR:::.logDiffTime(main = "Downloading PeakSet", t1 = tstart, verbose = verbose, logFile = logFile)
  ### Get peak set
  url_peaks = "https://zenodo.org/records/13142236/files/PeakSet.rds"
  
  download.file(
    url = url_peaks, 
    destfile = file.path(directory_to_store, "PeakSet.rds")
  ) 
  
  ArchR:::.logDiffTime(main = "Downloading project", t1 = tstart, verbose = verbose, logFile = logFile)
  ### Get ArchR project
  url_proj = "https://zenodo.org/records/13142236/files/ArchRProj_RWireX_Vignette.zip"
  
  download.file(
    url = url_proj, 
    destfile = file.path(directory_to_store, "ArchRProj_RWireX_Vignette.zip")
  )

  ArchR:::.logDiffTime(main = "Unzipping project", t1 = tstart, verbose = verbose, logFile = logFile)
  unzip(
    file.path(directory_to_store, "ArchRProj_RWireX_Vignette.zip"), 
    exdir = directory_to_store
  )

  ### Return local paths to ArchR project and peak set
  paths <- list(proj = file.path(directory_to_store, "ArchRProj_RWireX_Vignette"),
    peaks = file.path(directory_to_store, "PeakSet.rds"))
  
  ArchR:::.endLogging(logFile = logFile)
  return(paths)
}
