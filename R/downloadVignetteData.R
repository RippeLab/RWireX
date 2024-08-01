#' Download Data for RWireX Vignette
#' 
#' This function will download data for the vignette and store it in the specified location.
#' @param directory_to_store directory where ArchRProject will be stored
#' @export
downloadVignetteData <- function(directory_to_store){
  #Make Sure URL doesnt timeout
  oldTimeout <- getOption('timeout')
  options(timeout=150000)
  
  ### Get peak set
  url_peaks = "https://hub.dkfz.de/s/peey26JC9X7drAb/download?path=%2F&files=PeakSet.rds"
  
  download.file(
    url = url_peaks, 
    destfile = file.path(directory_to_store, "PeakSet.rds")
  ) 
  
  ### Get ArchR project
  url_proj = "https://hub.dkfz.de/s/peey26JC9X7drAb/download?path=%2F&files=ArchRProj_RWireX_Vignette.zip"
  
  download.file(
    url = url_proj, 
    destfile = file.path(directory_to_store, "ArchRProj_RWireX_Vignette.zip")
  )

  unzip(
    file.path(directory_to_store, "ArchRProj_RWireX_Vignette.zip"), 
    exdir = directory_to_store
  )

  ### Return local paths to ArchR project and peak set
  paths <- list(proj = file.path(directory_to_store, "ArchRProj_RWireX_Vignette"),
    peaks = file.path(directory_to_store, "PeakSet.rds"))

  return(paths)
}
