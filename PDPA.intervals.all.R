source("packages.R")

library(data.table)
library(coseg)

model.file.vec <- Sys.glob("../chip-seq-paper/chunks/H*/*/PDPA.model.RData")
PDPA.intervals.list <- list()
PDPA.intervals.raw <- list()
file.i <- 42
sample.i <- 24
file.i.vec <- seq_along(model.file.vec)
for(file.i in file.i.vec){
  model.file <- model.file.vec[[file.i]]
  chunk.id.dir <- dirname(paste(model.file))
  chunk.id <- basename(chunk.id.dir)
  set.dir <- dirname(chunk.id.dir)
  set.name <- basename(set.dir)
  cat(sprintf("%4d / %4d %s\n", file.i, length(model.file.vec), model.file))
  if(!file.exists(model.file)){
    stop("run PDPA.timings.R first")
  }else{
    load(model.file) #counts
    for(sample.id in names(PDPA.model)){
      model.list <- PDPA.model[[sample.id]]
      q.vec <- quantile(model.list$intervals.mat[, -(1:19)])
      PDPA.intervals.list[[paste(model.file, sample.id)]] <- data.table(
        n.data=model.list$n.data, t(q.vec))
      ## PDPA.intervals.raw[[paste(model.file, sample.id)]] <- with(model.list, {
      ##   data.table(
      ##     n.data,
      ##     intervals=as.integer(intervals.mat),
      ##     as.integer(row(intervals.mat)
    }
  }
}
PDPA.intervals.all <- do.call(rbind, PDPA.intervals.list)

save(PDPA.intervals.all, file="PDPA.intervals.all.RData")
