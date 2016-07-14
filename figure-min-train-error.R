source("packages.R")

load("../PeakSeg-paper/dp.peaks.error.RData")
load("PDPA.peaks.error.RData")

names(dp.peaks.error)
pdpa <- data.table(PDPA.peaks.error)
names(pdpa)[3] <- "param.name"
pdpa$algorithm <- "PDPA"
error.regions.list <- list(PDPA=pdpa)
for(chunk.name in names(dp.peaks.error)){
  one.chunk <- dp.peaks.error[[chunk.name]]
  error.regions.list[[chunk.name]] <- data.table(
    chunk.name, do.call(rbind, one.chunk),
    algorithm="cDPA")
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
table(train.error.wide[, PDPA-cDPA])
train.error.wide[order(cDPA-PDPA),]
sort(train.error.wide[cDPA<PDPA, table(chunk.name)])
train.error.wide[cDPA<PDPA & grepl("TDH", chunk.name),]
error.counts[chunk.name=="H3K4me3_PGP_immune/15" & sample.id=="McGill0079",]
error.counts[chunk.name=="H3K4me3_PGP_immune/15" & algorithm=="PDPA",]
error.counts[chunk.name=="H3K4me3_TDH_immune/1" & sample.id=="McGill0011",]
