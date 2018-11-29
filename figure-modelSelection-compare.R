source("packages.R")

chunk.name <- "H3K36me3_TDH_immune/1"
chunk.dir <- file.path(
  "../chip-seq-paper/chunks",
  chunk.name)
Segmentor.model.RData <- file.path(chunk.dir, "Segmentor.model.RData")
##load(Segmentor.model.RData)
sample.id <- "McGill0004"
counts.RData <- file.path(chunk.dir, "counts.RData")
load(counts.RData)
sample.counts <- data.table(counts)[sample.id, on=list(sample.id)]
uncompressed.vec <- sample.counts[, rep(coverage, chromEnd-chromStart)]
fit <- Segmentor(uncompressed.vec, Kmax=19, 
Cr <- SelectModel(model.list, penalty="oracle", keep=TRUE)


load("all.modelSelection.RData")
load("unsupervised.Segmentor.RData")

all.modelSelection[chunk.name== & segments==19 & algo=="PDPA", list(sample.id, min.lambda, max.lambda)]

oracle.Segmentor[["H3K36me3_TDH_immune/1"]]
