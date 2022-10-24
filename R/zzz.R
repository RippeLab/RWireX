.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is version ", packageVersion(pkgname), 
                        " of ", pkgname, ".\n
                        Please be mindful that the package is under development.\n 
                        If you encounter errors or bugs, PLEASE REPORT them to anastasiya.vladimirova@zoho.com or on slack.")
}