library(ggplot2)
library(directlabels)

load("PDPA.microbenchmark.RData")

ggplot()+
  geom_point(aes(log10(n.data), log10(time/1e9), color=expr),
             shape=1,
             data=PDPA.microbenchmark)

with.legend <- ggplot()+
  geom_point(aes(n.data, time/1e9, color=expr),
             shape=1,
             data=PDPA.microbenchmark)+
  ylab("seconds")+
  xlab("number of data to segment")

##with.labels <- direct.label(with.legend, "last.polygons")

pdf("figure-PDPA-microbenchmark.pdf", h=3)
print(with.legend)
dev.off()
