source("packages.R")

files <- Sys.glob("../chip-seq-paper/chunks/H*/*/Segmentor.model.RData")
pattern <-
  paste0("../chip-seq-paper/chunks/",
         "(?<set_name>.+?)",
         "/",
         "(?<chunk_id>[0-9]+)")
(matched <- str_match_named(files, pattern))

Segmentor.infeasible.error.list <- list()
Segmentor.model.list <- list()
Segmentor.loss.list <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  chunk.name <- paste0(set.name, "/", chunk.id)
  f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), f))
  load(f)
  Segmentor.model.list[[chunk.name]] <- Segmentor.model
  regions.RData <- sprintf("../chip-seq-paper/chunks/%s/regions.RData", chunk.name)
  load(regions.RData)
  regions.by.sample <- split(regions, regions$sample.id, drop=TRUE)
  for(sample.id in names(Segmentor.model)){
    sample.regions <- regions.by.sample[[sample.id]]
    fit <- Segmentor.model[[sample.id]]
    seg.vec <- seq(1, 19, by=2)
    Segmentor.loss.list[[paste(chunk.name, sample.id)]] <- data.table(
      set.name, chunk.id, chunk.name, sample.id,
      loss=fit$Likelihood[seg.vec],
      segments=seg.vec)
    for(seg.str in paste(seg.vec)){
      seg.dt <- data.table(fit$segments[[seg.str]])
      seg.dt[, change.after := c(diff(mean), Inf)]
      seg.dt[, change.before := c(-Inf, diff(mean))]
      seg.dt[, st.join := ifelse(
        change.before<0 & 0<change.after,
        "background", "peak")]
      seg.dt[, st.rm := ifelse(
        0<change.before & change.after<0,
        "background", "peak")]
      for(rule in c("join", "rm")){
        peak.dt <- if(seg.str=="1"){
          Peaks()
        }else{
          col.name <- paste0("st.", rule)
          status.vec <- seg.dt[[col.name]]
          seg.dt[, peak.i := cumsum(status.vec=="background")]
          seg.dt[status=="peak", list(
            chromStart=min(chromStart),
            chromEnd=max(chromEnd)
          ), by=list(peak.i)]
        }
        err.df <- PeakErrorChrom(peak.dt, sample.regions)
        id.str <- paste(rule, chunk.name, sample.id, seg.str)
        Segmentor.infeasible.error.list[[id.str]] <- data.table(
          rule, chunk.name, sample.id,
          segments=as.integer(seg.str),
          peaks=nrow(peak.dt), err.df)
      }
    }
  }
}
Segmentor.infeasible.error <- do.call(rbind, Segmentor.infeasible.error.list)
Segmentor.loss <- do.call(rbind, Segmentor.loss.list)

save(Segmentor.loss, Segmentor.infeasible.error, file="Segmentor.infeasible.error.RData")
