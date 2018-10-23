source("packages.R")

counts.RData.vec <- Sys.glob(file.path(
  "../chip-seq-paper/chunks",
  "H*", "*",
  "counts.RData"))

all.features <- list()
for(file.i in seq_along(counts.RData.vec)){
  counts.RData <- counts.RData.vec[[file.i]]
  load(counts.RData)
  chunk.name <- sub(".*chunks/", "", dirname(counts.RData))
  all.features[[chunk.name]] <- featureMatrix(counts, "sample.id", "coverage")
}

save(all.features, file="all.features.RData")
