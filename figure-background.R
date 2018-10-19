source("packages.R")

regions.RData.vec <- Sys.glob("../chip-seq-paper/chunks/*/*/regions.RData")
all.dt.list <- list()
for(regions.RData in regions.RData.vec){
  set.name <- basename(dirname(dirname(regions.RData)))
  load(regions.RData)
  regions.dt <- data.table(regions)
  counts.RData <- sub("regions.RData", "counts.RData", regions.RData)
  load(counts.RData)
  counts.dt <- data.table(counts)
  neg.regions <- regions.dt[annotation=="noPeaks", {
    data.table(
      sample.id,
      regionStart=chromStart,
      regionEnd=chromEnd,
      regionBases=chromEnd-chromStart)
  }]
  setkey(counts.dt, sample.id, chromStart, chromEnd)
  setkey(neg.regions, sample.id, regionStart, regionEnd)
  over.dt <- foverlaps(counts.dt, neg.regions, nomatch=0L)
  over.dt[chromStart < regionStart, chromStart := regionStart]
  over.dt[regionEnd < chromEnd, chromEnd := regionEnd]
  over.dt[, bases := chromEnd-chromStart]
  region.counts <- over.dt[, list(
    total.coverage=sum((chromEnd-chromStart)*coverage)
    ), by=list(sample.id, regionStart, regionEnd, regionBases)]
  region.counts[, mean.coverage := total.coverage/regionBases]
  all.dt.list[[regions.RData]] <- data.table(set.name, region.counts)
}
(all.dt <- do.call(rbind, all.dt.list))

all.dt[, show.set := gsub("_", "\n", set.name)]
all.dt[, sample := gsub("McGill0", "", sample.id)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(show.set ~ ., scales="free", space="free")+
  geom_point(aes(mean.coverage, sample), data=all.dt, shape=1)+
  scale_x_log10("mean coverage in genomic regions with noPeaks label", breaks=10^seq(-5, 5))
pdf("figure-background-all.pdf")
print(gg)
dev.off()

one <- all.dt[set.name=="H3K36me3_AM_immune"]
one.mean <- one[, list(
  mean=mean(mean.coverage)
  ), by=list(sample.id)][order(mean)]
one[, sample.fac := factor(sample.id, one.mean$sample.id)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(mean.coverage, sample.fac), data=one, shape=1)+
  scale_x_log10("mean coverage in genomic regions with noPeaks label", breaks=10^seq(-5, 5), limits=c(0.1, 10))+
  ylab("sample")
pdf("figure-background.pdf", 6, 3)
print(gg)
dev.off()

