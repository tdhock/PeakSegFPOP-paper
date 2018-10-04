source("jss-packages.R")

fit.dt <- readRDS("jss.disk.memory.rds")

fit.dt[, bench.seconds := time/1e9 ]

bench.stats <- fit.dt[, list(
  median=median(bench.seconds),
  q25=quantile(bench.seconds, 0.25),
  q75=quantile(bench.seconds, 0.75)
), by=list(bedGraph.lines, storage=expr)]

wide <- dcast(bench.stats, bedGraph.lines ~ storage, value.var="median")

ggplot()+
  geom_point(aes(
    log10(bedGraph.lines), disk/memory),
    data=wide)

storage.colors <- c(
  "#E41A1C",#red
  "#377EB8",#blue
  disk="#4DAF4A",#green
  memory="#984EA3",#purple
  "#FF7F00",#orange
  "#FFFF33", #yellow
  "#A65628",#brown
  "#F781BF",#pink
  "#999999")#grey

details.dt <- fit.dt[bedGraph.lines==106569]
leg <- ggplot()+
  theme_bw()+
  geom_point(aes(
    penalty, bench.seconds, color=expr),
    shape=1,
    data=details.dt)+
  scale_color_manual(values=storage.colors)+
  scale_x_log10(
    "lambda=penalty (log scale)",
    limits=c(NA, 10^5.75)
  )+
  scale_y_log10("seconds (log scale)")
dl <- direct.label(leg, list("last.qp", dl.trans(x=x+0.1)))
pdf("jss-figure-disk-memory-compare-speed-penalty.pdf", 3.3, 3)
print(dl)
dev.off()

range.dt <- details.dt[, list(
  min.seconds=min(bench.seconds),
  max.seconds=max(bench.seconds),
  bedGraph.lines=bedGraph.lines[1]
  )]
leg <- ggplot()+
  theme_bw()+
  scale_color_manual(values=storage.colors)+
  scale_fill_manual(values=storage.colors)+
  ## geom_segment(aes(
  ##   bedGraph.lines, min.seconds,
  ##   xend=bedGraph.lines, yend=max.seconds),
  ##   data=range.dt)+
  geom_ribbon(aes(
    bedGraph.lines, ymin=q25, ymax=q75, fill=storage),
    alpha=0.5,
    data=bench.stats)+
  geom_line(aes(
    bedGraph.lines, median, color=storage),
    data=bench.stats)+
  scale_y_log10("seconds (log scale)")+
  scale_x_log10(
    "N = number of data to segment (log scale)",
    limits=c(NA, 10^6.2))
dl <- direct.label(leg, list("last.qp", dl.trans(x=x+0.1)))
pdf("jss-figure-disk-memory-compare-speed.pdf", 3.3, 3)
print(dl)
dev.off()
