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
  works_with_R("3.2.2",
               tikzDevice="0.9",
               LambertW="0.6.2",
               microbenchmark="1.4.2",
               Segmentor3IsBack="2.0",
               geometry="0.3.6",
               "tdhock/memtime@8b80d7d2b151cf5e877b27138c2791eca365e7b1",
               ##"tdhock/coseg@073dcf152355e6de32236f9556d97dda29574bf5",
               "tdhock/PeakError@b0f0b4edc413176ebb183fc68f1504c9d86e3ef7",
               "tdhock/PeakSegDP@4e476f1ddf8b6252179b73dfec8c4cd616e9b5ad",
               ##"tdhock/animint@03735869af84629d269556442345b2ea506ab42a",
               "tdhock/directlabels@8a1f6f3501d5badf061d15abd23e4e42d5d32bbe",
               "Rdatatable/data.table@7515fbe6c6f60114da72067db44fbe78ecdbd8fb")
})
library(coseg)
options(
  tikzDocumentDeclaration=paste(
    "\\documentclass{article}",
    "\\usepackage{amsmath,amssymb,amsthm}"),
  tikzMetricsDictionary="tikzMetrics")
