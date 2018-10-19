source("packages.R")

library(penaltyLearning)

load("PDPA.peaks.error.RData")

selection.list <- list()
model.files <- Sys.glob("../chip-seq-paper/chunks/H*/*/PDPA.model.RData")
i.vec <- seq_along(model.files)
##i.vec <- 1:2
selected <- data.table(model.file.i=i.vec)[, {
  model.file <- model.files[[model.file.i]]
  chunk.name <- sub("../chip-seq-paper/chunks/", "", dirname(model.file))
  cat(sprintf("%4d / %4d %s\n", model.file.i, length(model.files), chunk.name))
  load(model.file)
  data.table(sample.id=names(PDPA.model))[, {
    seg.vec <- seq(1, 19, by=2)
    sample.model <- PDPA.model[[sample.id]]
    is.feasible <- sapply(seg.vec, function(p){
      seg.mean.vec <- sample.model$mean.mat[p, 1:p]
      all(diff(seg.mean.vec) != 0)
    })    
    model.dt <- data.table(
      chunk.name,
      peaks=(seg.vec-1)/2,
      loss=sample.model$loss.vec[seg.vec]
      )[is.feasible]
    modelSelection(model.dt, complexity="peaks")
  }, by=sample.id]
}, by=model.file.i]

error.dt <- data.table(PDPA.peaks.error)[, list(
  possible.fn=sum(possible.tp),
  possible.fp=sum(possible.fp),
  fp=sum(fp),
  fn=sum(fn),
  errors=sum(fn+fp),
  labels=.N
  ), by=list(chunk.name, sample.id, peaks=as.numeric(peaks))]

selected.errors <- error.dt[selected, on=list(chunk.name, sample.id, peaks)]
(selected.missing <- selected.errors[apply(is.na(selected.errors), 1, any)])
stopifnot(nrow(selected.missing)==0)

PDPA.targets <- targetIntervals(selected.errors, c("chunk.name", "sample.id"))

save(PDPA.targets, file="PDPA.targets.RData")

