source("packages.R")

load("../PeakSeg-paper/dp.peaks.error.RData")
load("PDPA.peaks.error.RData")
load("Segmentor.peaks.error.RData")

names(dp.peaks.error)
pdpa <- data.table(PDPA.peaks.error)
names(pdpa)[3] <- "param.name"
pdpa$algorithm <- "coseg"
Seg <- data.table(Segmentor.peaks.error)
names(Seg)[3] <- "param.name"
Seg$algorithm <- "Segmentor"
error.regions.list <- list(
  coseg=pdpa,
  Segmentor=Seg)
for(chunk.name in names(dp.peaks.error)){
  one.chunk <- dp.peaks.error[[chunk.name]]
  error.regions.list[[chunk.name]] <- data.table(
    chunk.name, do.call(rbind, one.chunk),
    algorithm="PeakSegDP")
}
error.regions <- do.call(rbind, error.regions.list)

error.counts <- error.regions[, list(
  errors=sum(fp+fn),
  regions=.N
  ), by=.(algorithm, chunk.name, sample.id, param.name)]

min.train.error <-
  error.counts[, .SD[which.min(errors),], by=.(algorithm, chunk.name, sample.id)]

total.train.error <- min.train.error[, list(
  errors=sum(errors),
  problems=.N
  ), by=.(algorithm)]

train.error.wide <-
  dcast(min.train.error, chunk.name + sample.id ~ algorithm, value.var="errors")
train.error.wide[PeakSegDP < coseg & coseg < Segmentor,]
## nice example to show the difference between algos.
##  8:  H3K4me3_TDH_immune/4 McGill0005         1         4     2
table(train.error.wide[, coseg-PeakSegDP])
table(train.error.wide[, coseg-Segmentor])
train.error.wide[order(PeakSegDP-coseg),]
sort(train.error.wide[PeakSegDP<coseg, table(chunk.name)])
train.error.wide[PeakSegDP<coseg & grepl("TDH", chunk.name),]
error.counts[chunk.name=="H3K4me3_PGP_immune/15" & sample.id=="McGill0079",]
error.counts[chunk.name=="H3K4me3_PGP_immune/15" & algorithm=="coseg",]
error.counts[chunk.name=="H3K4me3_TDH_immune/1" & sample.id=="McGill0011",]
