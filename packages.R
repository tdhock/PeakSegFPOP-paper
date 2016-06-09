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
               "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
               "tdhock/PeakSegDP@b9ae720406f07346150905541635d67c7c3c6e70",
               "hadley/scales@2c3edf45de56d617444dc38e47e0404173817886",
               "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
               "tdhock/animint@6b1c9e588b03f632cd39cdec9bbcfa730db9e889",
               "tdhock/directlabels@7b4b08a5dd0ab86e0b90902b3a233903ddd42311",
               data.table="1.9.6")
})
