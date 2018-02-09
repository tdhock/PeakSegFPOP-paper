### Write down what package versions work with your R code, and
### attempt to download and load those packages. The first argument is
### the version of R that you used, e.g. "3.0.2" and then the rest of
### the arguments are package versions. For
### CRAN/Bioconductor/R-Forge/etc packages, write
### e.g. RColorBrewer="1.0.5" and if RColorBrewer is not installed
### then we use install.packages to get the most recent version, and
### warn if the installed version is not the indicated version. For
### GitHub packages, write "user/repo@commit"
### e.g. "tdhock/animint@f877163cd181f390de3ef9a38bb8bdd0396d08a4" and
### we use install_github to get it, if necessary.
works_with_R <- function(Rvers,...){
  local.lib <- file.path(getwd(), "library")
  old.path.vec <- .libPaths()
  if(! local.lib %in% old.path.vec){
    dir.create(local.lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(local.lib)
  }
  pkg_ok_have <- function(pkg,ok,have){
    stopifnot(is.character(ok))
    if(!as.character(have) %in% ok){
      warning("works with ",pkg," version ",
              paste(ok,collapse=" or "),
              ", have ",have)
    }
  }
  pkg_ok_have("R",Rvers,getRversion())
  pkg.vers <- list(...)
  for(pkg.i in seq_along(pkg.vers)){
    vers <- pkg.vers[[pkg.i]]
    pkg <- if(is.null(names(pkg.vers))){
      ""
    }else{
      names(pkg.vers)[[pkg.i]]
    }
    if(pkg == ""){# Then it is from GitHub.
      ## suppressWarnings is quieter than quiet.
      if(!suppressWarnings(require(requireGitHub))){
        ## If requireGitHub is not available, then install it using
        ## devtools.
        if(!suppressWarnings(require(devtools))){
          install.packages("devtools")
          require(devtools)
        }
        install_github("tdhock/requireGitHub")
        require(requireGitHub)
      }
      requireGitHub(vers)
    }else{# it is from a CRAN-like repos.
      if(!suppressWarnings(require(pkg, character.only=TRUE))){
        install.packages(pkg)
      }
      pkg_ok_have(pkg, vers, packageVersion(pkg))
      library(pkg, character.only=TRUE)
    }
  }
}
suppressPackageStartupMessages({
  works_with_R(
    "3.4.3",
    tikzDevice="0.10.1",
    LambertW="0.6.4",
    microbenchmark="1.4.2.1",
    Segmentor3IsBack="2.0",
    lpSolveAPI="5.5.2.0.9",
    quadmod="2013.8.23",
    geometry="0.3.6",
    data.table="1.10.5",
    "tdhock/memtime@8b80d7d2b151cf5e877b27138c2791eca365e7b1",
    "tdhock/coseg@4a996702eb6988ad253fdfd2510df3feb392ff3e",
    "tdhock/PeakError@b0f0b4edc413176ebb183fc68f1504c9d86e3ef7",
    "tdhock/PeakSegDP@4e476f1ddf8b6252179b73dfec8c4cd616e9b5ad",
    "faizan-khan-iit/ggplot2@5fb99d0cece13239bbbc09c6b8a7da7f86ac58e2",
    "tdhock/animint@c0db9f34c525bec35c797ccdf8be9564b67c578c",
    "tdhock/directlabels@dcf34672129bf99a79ddfaceaef73236ae0f696d")
})
options(
  tikzDocumentDeclaration=paste(
    "\\documentclass[12pt]{article}",
    "\\usepackage{amsmath,amssymb,amsthm}"),
  tikzMetricsDictionary="tikzMetrics")

