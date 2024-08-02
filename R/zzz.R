RWireXDependency <- c(
  "parallel",
  "ggExtra", 
  "plotgardener"
)

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is version ", packageVersion(pkgname), 
                        " of ", pkgname, ".")
  pkgs <- RWireXDependency
  for(i in seq_along(pkgs)){
    packageStartupMessage("\tLoading Package : ", pkgs[i], " v", packageVersion(pkgs[i]))
    tryCatch({
      suppressPackageStartupMessages(require(pkgs[i], character.only=TRUE))
    }, error = function(e){
      packageStartupMessage("\tFailed To Load Package : ", pkgs[i], " v", packageVersion(pkgs[i]))
    })
  }
}