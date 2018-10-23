source("packages.R")

files <- Sys.glob("../chip-seq-paper/chunks/H*/*/Segmentor.model.RData")
pattern <-
  paste0("../chip-seq-paper/chunks/",
         "(?<set_name>.+?)",
         "/",
         "(?<chunk_id>[0-9]+)")
(matched <- str_match_named(files, pattern))

Segmentor.peaks.error.list <- list()
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
      segments=seg.vec)[as.integer(names(fit$peaks))+1L]
    for(peaks.str in names(fit$peaks)){
      peak.df <- fit$peaks[[peaks.str]]
      if(is.null(peak.df)){
        cat(paste(chunk.name, sample.id, peaks.str, "peaks\n"))
        peak.df <- Peaks()
      }
      err.df <- PeakErrorChrom(peak.df, sample.regions)
      Segmentor.peaks.error.list[[paste(chunk.name, sample.id, peaks.str)]] <- 
        data.table(chunk.name, sample.id, peaks=peaks.str, err.df)
    }
  }
}
Segmentor.peaks.error <- do.call(rbind, Segmentor.peaks.error.list)
Segmentor.loss <- do.call(rbind, Segmentor.loss.list)

save(Segmentor.loss, Segmentor.model.list, Segmentor.peaks.error, file="Segmentor.peaks.error.RData")
