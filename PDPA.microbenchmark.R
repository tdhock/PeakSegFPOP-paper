source("packages.R")

load("../chip-seq-paper/chunks/H3K36me3_AM_immune/2/counts.RData")
counts$count <- counts$coverage
by.sample <- split(counts, counts$sample.id)
sapply(by.sample, nrow)
one <- by.sample[[1]]
count.vec <- one$count
weight.vec <- with(one, chromEnd-chromStart)
max.segments <- 19L

time.list <- list()
for(n.data in seq(500, 6000, by=500)){
  some.data <- count.vec[1:n.data]
  some.weights <- weight.vec[1:n.data]
  some.times <- microbenchmark(
    coseg=PeakSegPDPA(some.data, some.weights, max.segments),
    cghseg=cghseg:::segmeanCO(some.data, max.segments),
    ##Segmentor3IsBack=stop("TODO"),
    PeakSegDP=cDPA(some.data, some.weights, max.segments),
    times=3)
  time.list[[paste(n.data)]] <- data.frame(n.data, some.times)
}
PDPA.microbenchmark <- do.call(rbind, time.list)

save(PDPA.microbenchmark, file="PDPA.microbenchmark.RData")
