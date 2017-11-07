library(PeakSegJoint)
data(chr7.peaks, envir=environment())
library(ggplot2)
ggplot()+
  geom_segment(aes(
    chromStart/1e3, sample.id,
    xend=chromEnd/1e3, yend=sample.id),
    data=chr7.peaks)

clustered <- multiClusterPeaks(chr7.peaks)
library(data.table)
clusters <- data.table(clustered)[, list(
  clusterStart=as.integer(median(chromStart)),
  clusterEnd=as.integer(median(chromEnd))
), by=cluster]
gg <- ggplot()+
  geom_segment(aes(
    chromStart/1e3, sample.id,
    color=factor(cluster),
    xend=chromEnd/1e3, yend=sample.id),
    data=clustered)+
  geom_segment(aes(
    clusterStart/1e3, "clusters",
    color=factor(cluster),
    xend=clusterEnd/1e3, yend="clusters"),
    data=clusters)+
  xlab("position on chromosome (kb = kilo bases)")

pdf("figure-peak-clustering.pdf", 7, 5)
print(gg)
dev.off()
