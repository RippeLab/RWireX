RWireXDependency <- c(
  "parallel",
  "ggExtra", 
  "plotgardener"
)

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is version ", packageVersion(pkgname), 
                        " of ", pkgname, ".\n
                        Please be mindful that the package is under development.\n 
                        If you encounter errors or bugs, PLEASE REPORT them to anastasiya.vladimirova@zoho.com or on slack.")
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