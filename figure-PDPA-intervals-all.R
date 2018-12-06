source("packages.R")

load("PDPA.intervals.all.RData")

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
  q25=quantile(mean.intervals, 0.25),
  min=quantile(mean.intervals, 0),
  max=max(max.intervals),
  median=median(mean.intervals),
  q75=quantile(mean.intervals, 0.75)
), by=list(mid, n.data=10^mid)]

text.dt <- median.dt[which.max(mid)][over.dt, on=list(mid), nomatch=0L][order(abs(mean.intervals-median))][1]

max.color <- "black"
med.color <- "black"
gg <- ggplot()+
  geom_text(aes(3, 200, label="max"),
            color=max.color)+
  geom_text(aes(4, 0, label="median and inter-quartile range"),
            color=med.color)+
  geom_point(aes(log10(n.data), max.intervals),
             color=max.color,
             shape=21,
             data=PDPA.intervals.all)+
  geom_line(aes(log10(n.data), mean.intervals),
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

text.color <- "red"
gg <- ggplot()+
  geom_text(aes(10^3, 10^2.3, label="max"),
            color=max.color)+
  geom_text(aes(10^4, 8, label="median, quartiles"),
            color=med.color)+
  geom_line(aes(n.data, max),
             color=max.color,
             data=median.dt)+
  geom_ribbon(aes(n.data, ymin=q25, ymax=q75),
              alpha=0.5,
              fill=med.color,
              data=median.dt)+
  geom_line(aes(n.data, median),
            color=med.color,
            data=median.dt)+
  theme_bw()+
  geom_point(aes(
    n.data, mean.intervals),
    shape=1,
    color=text.color,
    data=text.dt)+
  geom_text(aes(
    n.data, mean.intervals,
    label=sprintf("%.0f intervals", mean.intervals)),
    color=text.color,
    vjust=-1,
    hjust=1,
    data=text.dt)+
  theme(
    plot.margin=grid::unit(c(6, 12, 6, 6), "pt"),
    panel.margin=grid::unit(0, "lines"))+
  scale_x_log10(
    "$n$ = data points to segment",
    breaks=c(
      1e3, 1e4, 
      range(PDPA.intervals.all$n.data))
    )+
  scale_y_log10(
    "$I$ = intervals stored",
    breaks=c(
      10, 100, 
      max(PDPA.intervals.all$max.intervals),
      min(PDPA.intervals.all$mean.intervals)
      )
  )
print(gg)

tikz("figure-PDPA-intervals-log-log.tex", w=3, h=1.8)
print(gg)
dev.off()
pdf("figure-PDPA-intervals-log-log.pdf", w=3.3, h=1.8)
print(gg)
dev.off()
theme_bw()$plot.margin

gg <- ggplot()+
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

