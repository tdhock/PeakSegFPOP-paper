source("packages.R")

load("PDPA.infeasible.RData")

PDPA.infeasible.error.list <- list()
for(file.i in seq_along(PDPA.infeasible)){
  chunk.name <- names(PDPA.infeasible)[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(PDPA.infeasible), chunk.name))
  regions.RData <- sprintf("../chip-seq-paper/chunks/%s/regions.RData", chunk.name)
  load(regions.RData)
  regions.by.sample <- split(regions, regions$sample.id, drop=TRUE)
  peaks.by.sample <- PDPA.infeasible[[chunk.name]]
  for(sample.id in names(regions.by.sample)){
    sample.regions <- regions.by.sample[[sample.id]]
    sample.peaks <- peaks.by.sample[[sample.id]]
    for(peaks.str in names(sample.peaks)){
      peak.df <- sample.peaks[[peaks.str]]
      if(nrow(peak.df)==0){
        ##cat(paste(chunk.name, sample.id, peaks.str, "peaks\n"))
        peak.df <- Peaks()
      }
      err.df <- PeakErrorChrom(peak.df, sample.regions)
      PDPA.infeasible.error.list[[paste(chunk.name, sample.id, peaks.str)]] <- 
        data.table(
          chunk.name, sample.id,
          peaks=nrow(peak.df),
          segments=as.integer(peaks.str)*2+1,
          err.df)
    }
  }
}
PDPA.infeasible.error <- do.call(rbind, PDPA.infeasible.error.list)
save(PDPA.infeasible.error, file="PDPA.infeasible.error.RData")
