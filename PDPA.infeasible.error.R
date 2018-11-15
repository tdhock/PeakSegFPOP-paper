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
      for(rule in c("join", "remove")){
        peak.dt <- sample.peaks[[peaks.str]]
        rule.peaks <- if(nrow(peak.dt)==0){
          Peaks()
        }else{
          peak.dt[rule, on=list(rule)]
        }
        err.dt <- PeakErrorChrom(rule.peaks, sample.regions)
        PDPA.infeasible.error.list[[paste(chunk.name, sample.id, peaks.str,rule)]] <- 
          data.table(
            rule,
            chunk.name, sample.id,
            peaks=nrow(rule.peaks),
            segments=as.integer(peaks.str)*2+1,
            err.dt)
      }
    }
  }
}
PDPA.infeasible.error <- do.call(rbind, PDPA.infeasible.error.list)
save(PDPA.infeasible.error, file="PDPA.infeasible.error.RData")
