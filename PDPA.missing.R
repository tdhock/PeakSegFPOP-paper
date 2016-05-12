source("PeakSegPDPA.R")

load("dp.peaks.NA.RData")

options(warn=2)

for(row.i in 1:nrow(dp.peaks.NA)){
  sample.info <- dp.peaks.NA[row.i,]
  model.RData <- sample.info[, paste0(
    "data/", chunk.name,
    "/PDPA.missing/", sample.id,
    ".RData")]
  if(!file.exists(model.RData)){
    PDPA.missing.dir <- dirname(model.RData)
    dir.create(PDPA.missing.dir, showWarnings=FALSE)
    counts.RData <- paste0("data/", sample.info$chunk.name, "/counts.RData")
    (objs <- load(counts.RData))
    counts$bases <- with(counts, chromEnd-chromStart)
    counts$count <- counts$coverage
    sample.list <- split(counts, counts$sample.id)
    sample.id <- paste(sample.info$sample.id)
    compressed <- data.table(sample.list[[sample.id]])
    bases <- sum(compressed$bases)
    n.data <- nrow(compressed)
    cat(sprintf("%4d / %4d sample %s %d bases %d data\n",
                row.i, nrow(dp.peaks.NA),
                sample.id,
                bases, n.data))
    result <- list()
    seconds <- system.time({
      model.list <- PeakSegPDPAchrom(compressed, maxPeaks=9L)
    })[["elapsed"]]
    result$model <- model.list
    result$timing <- 
      data.frame(seconds, data=n.data, bases)
    save(result, file=model.RData)
  }
}
