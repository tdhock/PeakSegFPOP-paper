library(parallel)

maxPeaks <- 9L
argv <- "/home/tdhock/projects/PeakSegFPOP-paper/data/H3K36me3_AM_immune/11/counts.RData"
argv <- commandArgs(trailingOnly=TRUE)

counts.RData <- argv[1]
print(counts.RData)
chunk.id.dir <- dirname(paste(counts.RData))
chunk.id <- basename(chunk.id.dir)
set.dir <- dirname(chunk.id.dir)
set.name <- basename(set.dir)
data.dir <- dirname(set.dir)
proj.dir <- dirname(data.dir)
source(file.path(proj.dir, "packages.R"))
source(file.path(proj.dir, "PeakSegPDPA.R"))

time.f <- sub("counts", "PDPA.timing", counts.RData)
model.f <- sub("counts", "PDPA.model", counts.RData)
load(counts.RData) #counts
counts$bases <- with(counts, chromEnd-chromStart)
counts$count <- counts$coverage
sample.list <- split(counts, counts$sample.id, drop=TRUE)
sample.ids <- names(sample.list)

sample.results <- mclapply(seq_along(sample.ids), function(sample.i){
  sample.id <- sample.ids[[sample.i]]
  compressed <- data.table(sample.list[[sample.id]])
  bases <- sum(compressed$bases)
  n.data <- nrow(compressed)
  cat(sprintf("%4d / %4d sample %s %d bases %d data\n",
              sample.i, length(sample.list),
              sample.id,
              bases, n.data))
  result <- list()
  seconds <- system.time({
    model.list <- PeakSegPDPAchrom(compressed, maxPeaks=maxPeaks)
  })[["elapsed"]]
  result$model <- model.list
  result$timing <- 
    data.frame(set.name, chunk.id, sample.id,
               seconds, data=n.data, bases)
  result
})

names(sample.results) <- sample.ids
PDPA.model <- lapply(sample.results, "[[", "model")
timing.list <- lapply(sample.results, "[[", "timing")
PDPA.timing <- do.call(rbind, timing.list)
## ggplot()+
##   geom_point(aes(data, seconds),
##              data=dp.timing, pch=1)+
##   scale_x_log10()+
##   scale_y_log10()
save(PDPA.model, file=model.f)
save(PDPA.timing, file=time.f)
