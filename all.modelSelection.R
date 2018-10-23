source("packages.R")

load("PDPA.infeasible.RData")
load("PDPA.infeasible.error.RData")
load("Segmentor.peaks.error.RData")
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

all.error <- rbind(
  data.table(algo="CDPA", CDPA.error),
  data.table(algo="PDPA", Segmentor.peaks.error),
  data.table(algo="GPDPA", PDPA.infeasible.error))

all.totals <- all.error[, list(
  total.fp=sum(fp),
  total.fn=sum(fn),
  total.errors=sum(fp+fn),
  possible.fp=sum(possible.fp),
  possible.fn=sum(possible.tp),
  labels=.N
  ), by=list(algo, chunk.name, sample.id, segments=as.integer(paste(peaks))*2+1)]

cbind(PDPA.loss[chunk.name=="H3K4me3_TDH_immune/1" & sample.id=="McGill0322", PoissonLoss],
Segmentor.loss[chunk.name=="H3K4me3_TDH_immune/1" & sample.id=="McGill0322", loss],
dp.loss[chunk.name=="H3K4me3_TDH_immune/1" & sample.id=="McGill0322", error])

all.loss <- rbind(
  PDPA.loss[, data.table(
    set.name, chunk.id, chunk.name, sample.id,
    algo="GPDPA", segments, loss=PoissonLoss)],
  Segmentor.loss[, data.table(
    set.name, chunk.id, chunk.name, sample.id,
    algo="PDPA", segments, loss)],
  dp.loss[, data.table(
    set.name, chunk.id, chunk.name, sample.id,
    algo="CDPA", segments, loss=error)])[all.totals, on=list(
      algo, chunk.name, sample.id, segments)]

all.modelSelection <- all.loss[, {
  cm <- cummin(loss)
  modelSelection(
    .SD[loss==cm],
    complexity="segments")
}, by=list(set.name, chunk.id, chunk.name, sample.id, algo)]

save(all.modelSelection, file="all.modelSelection.RData")
