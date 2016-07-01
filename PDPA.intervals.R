source("packages.R")

load("data/H3K36me3_AM_immune/2/counts.RData")

counts$count <- counts$coverage
by.sample <- split(counts, counts$sample.id)
n.data.vec <- sapply(by.sample, nrow)
one <- by.sample[[which.max(n.data.vec)]]
count.vec <- one$count
weight.vec <- with(one, chromEnd-chromStart)
max.segments <- 19L

fit <- PeakSegPDPA(count.vec, weight.vec, max.segments)

PDPA.intervals <- data.frame(
  segments=as.numeric(row(fit$intervals.mat)),
  data=as.numeric(col(fit$intervals.mat)),
  intervals=as.numeric(fit$intervals.mat))

save(PDPA.intervals, fit, file="PDPA.intervals.RData")
