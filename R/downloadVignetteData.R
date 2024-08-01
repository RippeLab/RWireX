#' Download Data for RWireX Vignette
#' 
#' This function will download data for the vignette and store it in the specified location.
#' @param directory_to_store directory where ArchRProject will be stored
#' @export
downloadVignetteData <- function(directory_to_store){
  #Make Sure URL doesnt timeout
  oldTimeout <- getOption('timeout')
  options(timeout=150000)
  
  project_dir = paste0(directory_to_store, "/Project")
  dir.create(project_dir, showWarnings = FALSE)
  
  url_peaks = "https://hub.dkfz.de/s/wGRBSbgATPFmmPA/download?path=%2F&files=PeakSet.rds"
  
  download.file(
    url = url_peaks, 
    destfile = file.path(directory_to_store, "PeakSet.Rds")
  ) 
  
  arrow_url = "https://hub.dkfz.de/s/wGRBSbgATPFmmPA/download?path=%2FExemplary_ArchR_Project&files=ArrowFiles"
  download.file(
    url = arrow_url, 
    destfile = file.path(project_dir, "ArrowFiles")
  ) 
  
  embedding_url = "https://hub.dkfz.de/s/wGRBSbgATPFmmPA/download?path=%2FExemplary_ArchR_Project&files=Embeddings"
  download.file(
    url = save_project_url, 
    destfile = file.path(project_dir, "Embeddings")
  ) 
  
  iterative_url = "https://hub.dkfz.de/s/wGRBSbgATPFmmPA/download?path=%2FExemplary_ArchR_Project&files=IterativeLSI"
  download.file(
    url = save_project_url, 
    destfile = file.path(project_dir, "IterativeLSI")
  ) 
  
  save_project_url = "https://hub.dkfz.de/s/wGRBSbgATPFmmPA/download?path=%2FExemplary_ArchR_Project&files=Save-ArchR-Project.rds"
  download.file(
    url = save_project_url, 
    destfile = file.path(project_dir, "Save-ArchR-Project.rds")
  ) 

}
