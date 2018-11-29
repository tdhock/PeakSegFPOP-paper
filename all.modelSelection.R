source("packages.R")

load("PDPA.infeasible.RData")
load("PDPA.infeasible.error.RData")
##load("Segmentor.peaks.error.RData")
load("Segmentor.infeasible.error.RData")
load("dp.peaks.RData")
load("dp.peaks.error.RData")

CDPA.error.list <- list()
for(chunk.name in names(dp.peaks.error)){
  chunk.error <- dp.peaks.error[[chunk.name]]
  for(cell.type in names(chunk.error)){
    type.dt <- data.table(
      chunk.name,
      chunk.error[[cell.type]])
    names(type.dt)[3] <- "peaks"
    CDPA.error.list[[paste(chunk.name, cell.type)]] <- type.dt
  }
}
CDPA.error <- do.call(rbind, CDPA.error.list) 

col.name.vec <- c(
  "chunk.name", "sample.id", "peaks", "segments", "chromStart", 
  "chromEnd", "annotation", "tp", "possible.tp", "fp", "possible.fp", 
  "fp.status", "fn", "fn.status", "status")
CDPA.error[, segments := as.integer(paste(peaks))*2+1]
all.error <- rbind(
  data.table(algo="CDPA", CDPA.error[, ..col.name.vec]),
  data.table(algo="PDPA", Segmentor.infeasible.error[rule=="rm", ..col.name.vec]),
  data.table(algo="GPDPA", PDPA.infeasible.error[rule=="remove", ..col.name.vec]))

all.totals <- all.error[, list(
  total.fp=sum(fp),
  total.fn=sum(fn),
  total.errors=sum(fp+fn),
  possible.fp=sum(possible.fp),
  possible.fn=sum(possible.tp),
  labels=.N
  ), by=list(algo, chunk.name, sample.id, segments, peaks)]

all.loss <- rbind(
  PDPA.loss[, data.table(
    set.name, chunk.id, chunk.name, sample.id,
    algo="GPDPA", segments, peaks, loss=PoissonLoss)],
  Segmentor.loss[, data.table(
    set.name, chunk.id, chunk.name, sample.id,
    algo="PDPA", segments, peaks, loss)],
  dp.loss[, data.table(
    set.name, chunk.id, chunk.name, sample.id,
    algo="CDPA", segments, peaks, loss=error)])[all.totals, on=list(
      algo, chunk.name, sample.id, segments)]

all.modelSelection <- all.loss[, {
  cm <- cummin(loss)
  modelSelection(
    .SD[loss==cm],
    complexity="segments")
}, by=list(set.name, chunk.id, chunk.name, sample.id, algo)]

class.vec <- sapply(all.modelSelection, class)
if(any(class.vec=="list")){
  print(class.vec)
  stop("list column!")
}

save(all.modelSelection, file="all.modelSelection.RData")
