source("packages.R")

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

gg <- ggplot()+
  geom_text(aes(3, 200, label="max"),
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
    "intervals stored")
tikz("figure-PDPA-intervals-small.tex", w=3.3, h=2.5)
print(gg)
dev.off()

gg <- ggplot()+
  geom_text(aes(10^3, 10^2.3, label="max"),
            color="blue")+
  geom_text(aes(10^4, 10^0.6, label="median, inter-quartile range"),
            color="black")+
  geom_point(aes(n.data, `100%`),
             color="blue",
             shape=21,
             data=PDPA.intervals.all)+
  geom_ribbon(aes(n.data, ymin=`25%`, ymax=`75%`),
              alpha=0.5,
              data=PDPA.intervals.all)+
  geom_line(aes(n.data, `50%`),
            data=PDPA.intervals.all)+
  theme_bw()+
  theme(
    plot.margin=grid::unit(c(6, 12, 6, 6), "pt"),
    panel.margin=grid::unit(0, "lines"))+
  scale_x_log10(
    "data points to segment (log scale)",
    breaks=c(
      1e3, 1e4, 
      range(PDPA.intervals.all$n.data))
    )+
  scale_y_log10(
    "intervals stored (log scale)",
    breaks=c(
      10, 100, 
      max(PDPA.intervals.all[["100%"]]),
      min(PDPA.intervals.all[["25%"]])
      )
    )
tikz("figure-PDPA-intervals-log-log.tex", w=3.3, h=2.5)
print(gg)
dev.off()
pdf("figure-PDPA-intervals-log-log.pdf", w=3.3, h=2.5)
print(gg)
dev.off()
theme_bw()$plot.margin

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

pdf("figure-PDPA-intervals-all.pdf", w=5, h=4)
print(gg)
dev.off()

