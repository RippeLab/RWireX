% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloadVignetteData.R
\name{downloadVignetteData}
\alias{downloadVignetteData}
\title{Download Data for RWireX Vignette}
\usage{
downloadVignetteData(
  directory_to_store,
  logFile = createLogFile("downloadVignetteData"),
  verbose = TRUE
)
}
\arguments{
\item{directory_to_store}{directory where ArchRProject will be stored}

\item{logFile}{The path to a file to be used for logging RWireX output.}

\item{verbose}{A boolean value that determines whether standard output should be printed.}
}
\description{
This function will download data for the vignette and store it in the specified location.
}
