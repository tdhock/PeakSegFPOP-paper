library(ggplot2)
library(data.table)

load("PDPA.intervals.all.RData")

gg <- ggplot()+
  ggtitle(paste(
    "2752 segmentation problems",
    "max segments = 19",
    sep=", "))+
  geom_text(aes(3, 200, label="max(intervals)"),
            color="blue")+
  geom_text(aes(4, 0, label="median and inter-quartile range"),
            color="black")+
  geom_point(aes(log10(n.data), `100%`),
             color="blue",
             shape=21,
             data=PDPA.intervals.all)+
  geom_ribbon(aes(log10(n.data), ymin=`25%`, ymax=`75%`),
              alpha=0.5,
              data=PDPA.intervals.all)+
  geom_line(aes(log10(n.data), `50%`),
            data=PDPA.intervals.all)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_x_continuous(
    "log10(data points to segment)")+
  scale_y_continuous(
    "intervals stored by the algorithm")

PDPA.intervals.all[, log10.data := log10(n.data)]
PDPA.intervals.all[, log10.data1 := log10(n.data)]
log.range.vec <- range(PDPA.intervals.all$log10.data)
n.boxes <- 20
n.edges <- n.boxes+1
box.edge.vec <- seq(log.range.vec[1], log.range.vec[2], l=n.edges)
box.dt <- data.table(
  left=box.edge.vec[-n.edges],
  right=box.edge.vec[-1])
box.dt[, mid := (left+right)/2]
setkey(box.dt, left, right)
setkey(PDPA.intervals.all, log10.data, log10.data1)
over.dt <- foverlaps(PDPA.intervals.all, box.dt)
stopifnot(nrow(over.dt)==nrow(PDPA.intervals.all))
median.dt <- over.dt[, list(
  problems=.N,
  q25=quantile(`50%`, 0.25),
  min=quantile(`50%`, 0),
  max=quantile(`50%`, 1),
  median=median(`50%`),
  q75=quantile(`50%`, 0.75)
  ), by=mid]

ggplot()+
  ggtitle(paste(
    "2752 segmentation problems",
    "max segments = 19",
    sep=", "))+
  geom_ribbon(aes(mid, ymin=min, ymax=max),
              alpha=0.5,
              data=median.dt)+
  geom_line(aes(mid, median),
            data=median.dt)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_x_continuous(
    "log10(data points to segment)")+
  scale_y_continuous("median(intervals)")

pdf("figure-PDPA-intervals-all.pdf", w=7, h=5)
print(gg)
dev.off()

