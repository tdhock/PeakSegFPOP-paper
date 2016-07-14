source("packages.R")

load("PDPA.peaks.RData")

PDPA.peaks.error.list <- list()
for(file.i in seq_along(PDPA.peaks)){
  chunk.name <- names(PDPA.peaks)[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(PDPA.peaks), chunk.name))
  regions.RData <- sprintf("data/%s/regions.RData", chunk.name)
  load(regions.RData)
  regions.by.sample <- split(regions, regions$sample.id, drop=TRUE)
  peaks.by.sample <- PDPA.peaks[[chunk.name]]
  for(sample.id in names(regions.by.sample)){
    sample.regions <- regions.by.sample[[sample.id]]
    sample.peaks <- peaks.by.sample[[sample.id]]
    for(peaks.str in names(sample.peaks)){
      peak.df <- sample.peaks[[peaks.str]]
      if(is.null(peak.df)){
        cat(paste(chunk.name, sample.id, peaks.str, "peaks\n"))
        peak.df <- Peaks()
      }
      err.df <- PeakErrorChrom(peak.df, sample.regions)
      PDPA.peaks.error.list[[paste(chunk.name, sample.id, peaks.str)]] <- 
        data.table(chunk.name, sample.id, peaks=peaks.str, err.df)
    }
  }
}
PDPA.peaks.error <- do.call(rbind, PDPA.peaks.error.list)
save(PDPA.peaks.error, file="PDPA.peaks.error.RData")
