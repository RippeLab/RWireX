% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCoAccessibleDomains.R
\name{getCoAccessibleDomains}
\alias{getCoAccessibleDomains}
\title{getCoAccessibleDomains}
\usage{
getCoAccessibleDomains(
  CoAccessibility,
  enrich_cutoff = 0.9,
  logFile = createLogFile("getCoAccessibleDomains"),
  verbose = TRUE
)
}
\arguments{
\item{CoAccessibility}{A GRanges object with co-accessible links from getCoAccessibility function.}

\item{enrich_cutoff}{Quantile cutoff to determine co-accessibility enrichment.}

\item{logFile}{The path to a file to be used for logging RWireX output.}

\item{verbose}{A boolean value that determines whether standard output should be printed.}
}
\value{
GRanges object
}
\description{
This function returns a GRanges object of co-accessibly enriched domains.
}
\examples{
Coming soon.
}
\keyword{archr}
\keyword{co-accessibility}
\keyword{domain}
\keyword{rwirex}
