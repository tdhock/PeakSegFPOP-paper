source("packages.R")

load("PDPA.microbenchmark.RData")

log.with.legend <- ggplot()+
  geom_point(aes(n.data, time/1e9, color=expr),
             shape=1,
             data=PDPA.microbenchmark)+
  scale_y_log10("seconds")+
  scale_x_log10("number of data to segment")
log.with.labels <- direct.label(log.with.legend, "last.polygons")+
  scale_y_log10("seconds", breaks=c(0.01, 0.1, 1, 10), limits=c(0.01, 10))+
  scale_x_log10(
    "number of data to segment",
    breaks=c(500, 1000, 5000, 6000),
    limits=c(500, 8000))
pdf("figure-PDPA-microbenchmark-log.pdf", h=3)
print(log.with.labels)
dev.off()

with.legend <- ggplot()+
  geom_point(aes(n.data, time/1e9, color=expr),
             shape=1,
             data=PDPA.microbenchmark)+
  ylab("seconds")+
  xlab("number of data to segment")
with.labels <- direct.label(with.legend, "last.polygons")+
  scale_x_continuous(
    "number of data to segment",
    breaks=c(500, 1000, 5000, 6000),
    limits=c(500, 6800))

pdf("figure-PDPA-microbenchmark.pdf", h=3)
print(with.labels)
dev.off()
