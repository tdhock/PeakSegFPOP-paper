##  47 /  129 chunks   26 /   27 sample McGill0091 34798 bases 1472 data
library(data.table)
chunk.name <- "H3K4me3_PGP_immune/24"
chunk.dir <- file.path("data", chunk.name)
load(file.path(chunk.dir, "counts.RData"))
H3K4me3_PGP_immune_chunk24 <- counts
## save(H3K4me3_PGP_immune_chunk24, file="~/R/coseg/data/H3K4me3_PGP_immune_chunk24.RData")
## prompt(H3K4me3_PGP_immune_chunk24, file="~/R/coseg/man/H3K4me3_PGP_immune_chunk24.Rd")
sample.id <- "McGill0009"
counts.by.sample <- split(counts, counts$sample.id)
sample.counts <- data.table(counts.by.sample[[sample.id]])[1:223,]
ggplot()+
  geom_step(aes(chromStart/1e3, coverage),
            data=sample.counts)
library(PeakSegDP)
library(coseg)
data.vec <- sample.counts$coverage
weight.vec <- sample.counts[, chromEnd-chromStart]
max.segments <- 19L
cdpa <- cDPA(data.vec, weight.vec, max.segments)
pdpa <- PeakSegPDPA(data.vec, weight.vec, max.segments)
compare.dt <- data.table(
  cdpa=as.numeric(cdpa$loss),
  pdpa=as.numeric(pdpa$cost.mat),
  segments=as.integer(row(cdpa$loss)),
  data.i=as.integer(col(cdpa$loss)))
compare.dt[cdpa+1e-6<pdpa,]

sample.counts$count <- sample.counts$coverage
pdpa.list <- PeakSegPDPAchrom(sample.counts, 9L)
cdpa.list <- PeakSegDP(sample.counts, 9L)
peaks.str <- "6"
cdpa.list$peaks[[peaks.str]]
subset(pdpa.list$segments, peaks==peaks.str & status=="peak")
