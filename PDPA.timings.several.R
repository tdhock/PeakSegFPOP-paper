proj.dir <- getwd()
counts.RData.vec <- Sys.glob("../chip-seq-paper/chunks/H*/*/counts.RData")
path.to.R <- R.home(file.path("bin","R"))
commands <-
  sprintf("%s --vanilla --args '%s' < %s",
          path.to.R,
          normalizePath(counts.RData.vec, mustWork=TRUE),
          normalizePath("PDPA.timings.one.R", mustWork=TRUE))

system(commands[1]) ## TEST if one of these works without qsub!
for(cmd in commands){
  f <- tempfile()
  qsubcmd <- sprintf("qsub %s",f)
  writeLines(cmd,f)
  system(qsubcmd)
}
