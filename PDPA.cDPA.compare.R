source("packages.R")

PDPA.RData.vec <- Sys.glob("../chip-seq-paper/chunks/H3K*/*/PDPA.model.RData")
loss.list <- list()
pdpa.list <- list()
all.loss.list <- list()
dp.list <- list()
seg.vec <- seq(1, 19, by=2)
peaks.vec <- (seg.vec - 1)/2
for(file.i in seq_along(PDPA.RData.vec)){
  PDPA.RData <- PDPA.RData.vec[[file.i]]
  (objs <- load(PDPA.RData))
  chunk.dir <- dirname(PDPA.RData)
  chunk.name <- sub("^data/", "", chunk.dir)
  load(file.path(chunk.dir, "dp.model.reverse.RData"))
  dp.model.reverse <- dp.model
  load(file.path(chunk.dir, "dp.model.RData"))
  load(file.path(chunk.dir, "Segmentor.model.RData"))
  load(file.path(chunk.dir, "counts.RData"))
  counts.by.sample <- split(counts, counts$sample.id)
  for(sid in names(dp.model)){
    sample.counts <- counts.by.sample[[sid]]
    PDPA <- PDPA.model[[sid]]
    Seg <- Segmentor.model[[sid]]
    PDPA.feasible <- sapply(seg.vec, function(i){
      mean.vec <- PDPA$mean.mat[i, 1:i]
      all(diff(mean.vec) != 0)
    })
    Seg.feasible <- sapply(Seg$segments, function(df){
      all(cumsum(sign(diff(df$mean))) %in% c(0,1))
    })
    err.dt <- data.table(
      segments=seg.vec,
      peaks=peaks.vec,
      PDPA=with(PDPA, cost.mat[seg.vec, n.data]),
      PDPA.feasible,
      Segmentor=Seg$Likelihood[seg.vec],
      Seg.feasible=Seg.feasible[seg.vec],
      dp.fwd=NA_real_,
      dp.rev=NA_real_)
    setkey(err.dt, segments)
    dp <- dp.model[[sid]]
    err.dt[J(dp$error$segments), dp.fwd := dp$error$error]
    r <- dp.model.reverse[[sid]]
    err.dt[J(r$error$segments), dp.rev := r$error$error]
    better.dt <- err.dt[PDPA + 1e-6 < dp.fwd & PDPA.feasible,]
    cat(sprintf("%4d / %4d chunks %s %d better peaks\n",
                file.i, length(PDPA.RData.vec),
                sid, nrow(better.dt)))
    all.loss.list[[paste(chunk.name, sid)]] <- data.table(
      chunk.name, sample.id=sid, err.dt)
    for(row.i in seq_along(better.dt$segments)){
      one <- better.dt[row.i, ]
      change.vec <- PDPA$ends.mat[one$segments, 2:one$segments]
      first <- c(1, change.vec+1)
      last <- c(change.vec, PDPA$n.data)
      seg.dt <- data.table(
        first,
        last,
        chromStart=sample.counts$chromStart[first],
        chromEnd=sample.counts$chromEnd[last],
        status=rep(c("background", "peak"), l=one$segments))
      PDPA.peaks <- seg.dt[status=="peak",]
      dp.peaks <- dp$peaks[[paste(one$peaks)]]
      stopifnot(identical(nrow(PDPA.peaks), nrow(dp.peaks)))
      PDPA.peaks[, diffStart := chromStart-dp.peaks$chromStart]
      PDPA.peaks[, diffEnd := chromEnd-dp.peaks$chromEnd]
      diffBases <- PDPA.peaks[, sum(abs(diffStart)+abs(diffEnd))]
      loss.list[[paste(chunk.name, sid, row.i)]] <- data.table(
        chunk.name, sample.id=sid, one, diffBases)
      pdpa.list[[paste(chunk.name, sid, row.i)]] <- data.table(
        chunk.name, sample.id=sid, one, PDPA.peaks)
      dp.list[[paste(chunk.name, sid, row.i)]] <- data.table(
        chunk.name, sample.id=sid, one, dp.peaks)
    }
  }
}

PDPA.cDPA.compare <- list(
  all.loss=do.call(rbind, all.loss.list),
  loss=do.call(rbind, loss.list),
  pdpa=do.call(rbind, pdpa.list),
  dp=do.call(rbind, dp.list))
save(PDPA.cDPA.compare, file="PDPA.cDPA.compare.RData")


