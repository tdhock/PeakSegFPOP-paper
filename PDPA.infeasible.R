source("packages.R")

files <- Sys.glob("../chip-seq-paper/chunks/H*/*/PDPA.model.RData")

pattern <- paste0(
  "chunks/",
  "(?<set_name>.+?)",
  "/",
  "(?<chunk_id>[0-9]+)")
matched <- str_match_named(files, pattern)
PDPA.infeasible <- list()
PDPA.loss.list <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  regions.str <- paste0(set.name, "/", chunk.id)
  f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), f))
  load(f)
  count.file <- sub("PDPA.model", "counts", f)
  load(count.file)
  counts.by.sample <- split(counts, counts$sample.id)
  for(sample.id in names(PDPA.model)){
    fit <- PDPA.model[[sample.id]]
    count.df <- counts.by.sample[[sample.id]]
    seg.vec <- seq(1, 19, by=2)
    loss.df <- data.frame(
      segments=seg.vec,
      peaks=(seg.vec-1)/2,
      PoissonLoss=fit$loss.vec[seg.vec])
    PDPA.loss.list[[paste(regions.str, sample.id)]] <- data.table(
      set.name, chunk.id, chunk.name=regions.str, sample.id,
      loss.df)
    peaks.list <- list()
    for(n.segments in loss.df$segments){
      n.peaks <- (n.segments-1)/2
      peaks.list[[paste(n.peaks)]] <- if(n.segments==1){
        data.table()
      }else{
        mean.vec <- fit$mean.mat[n.segments, 1:n.segments]
        diff.vec <- diff(mean.vec)
        break.vec <- if(n.segments==1){
          c()
        }else{
          fit$ends.mat[n.segments, 2:n.segments]
        }
        first <- c(1, break.vec+1)
        last <- c(break.vec, nrow(count.df))
        status.str <- rep(c("background", "peak"), l=n.segments)
        peak.dt <- data.table(
          mean=mean.vec,
          first,
          last,
          is.peak=(seq_along(mean.vec)-1) %% 2,
          chromStart=count.df$chromStart[first],
          chromEnd=count.df$chromEnd[last],
          diff.before=c(Inf, diff.vec),
          diff.after=c(diff.vec, Inf),
          peaks=n.peaks,
          segments=n.segments)
        remove.dt <- peak.dt[is.peak==1 & diff.before != 0 & diff.after != 0]
        peak.dt[is.peak==0 & (diff.after==0|diff.before==0), is.peak := 1]
        peak.dt[, peak.i := cumsum(is.peak==0)]
        join.dt <- peak.dt[is.peak==1, list(
          chromStart=min(chromStart),
          chromEnd=max(chromEnd)
        ), by=list(peak.i, peaks, segments)]
        rbind(
          join.dt[, data.table(
            rule="join", segments, chromStart, chromEnd)],
          remove.dt[, data.table(
            rule="remove", segments, chromStart, chromEnd)])
      }
    }
    PDPA.infeasible[[regions.str]][[sample.id]] <- peaks.list
  }
}
PDPA.loss <- do.call(rbind, PDPA.loss.list)

save(PDPA.loss, PDPA.infeasible, file="PDPA.infeasible.RData")
