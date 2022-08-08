.onAttach <- function(libname, pkgname) {

    ver <- "2.0"

    loadmsg <- paste0("\nLoading the 'orchaRd' package (version ", ver, "). For an\nintroduction and vignette to the package please see: https://daniel1noble.github.io/orchaRd/\n")

    packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)

  }
