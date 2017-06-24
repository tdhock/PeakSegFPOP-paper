source("packages.R")

load("PDPAInf.intervals.RData")
load("PDPAInf.timings.RData")
load("PDPA.timings.RData")

all.timings <- rbind(
  data.table(algo.type="PDPAInf", PDPAInf.timings),
  data.table(algo.type="PDPA", PDPA.timings))

time.int <- all.timings[PDPAInf.intervals, on=list(algo.type, set.name, chunk.id, sample.id)]
some.missing <- time.int[apply(is.na(time.int), 1, any)]
stopifnot(nrow(some.missing)==0)

ggplot(time.int)+
  geom_line(aes(data, max.intervals, color=algo.type))+
  scale_y_log10()+
  scale_x_log10()

wide.dt <- dcast(
  time.int,
  set.name + chunk.id + sample.id ~ algo.type,
  value.var=c("mean.intervals", "max.intervals", "seconds"))
wide.dt[mean.intervals_PDPA==mean.intervals_PDPAInf]
wide.dt[max.intervals_PDPA==max.intervals_PDPAInf]

ggplot()+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="grey", slope=1, intercept=0)+
  geom_point(aes(
    mean.intervals_PDPA, mean.intervals_PDPAInf),
             shape=1,
             data=wide.dt)+
  coord_equal()

ggplot()+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="grey", slope=1, intercept=0)+
  geom_point(aes(
    max.intervals_PDPA, max.intervals_PDPAInf),
             shape=1,
             data=wide.dt)+
  coord_equal()

ggplot()+
  theme_bw()+
  geom_abline(color="grey", slope=1, intercept=0)+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(
    seconds_PDPA, seconds_PDPAInf),
             shape=1,
             data=wide.dt)+
  coord_equal()

time.int[, log10.data := log10(data)]
time.int[, log10.data1 := log10(data)]
log.range.vec <- range(time.int$log10.data)
n.boxes <- 20
n.edges <- n.boxes+1
box.edge.vec <- seq(log.range.vec[1], log.range.vec[2], l=n.edges)
box.dt <- data.table(
  left=box.edge.vec[-n.edges],
  right=box.edge.vec[-1])
box.dt[, mid := (left+right)/2]
setkey(box.dt, left, right)
setkey(time.int, log10.data, log10.data1)
over.dt <- foverlaps(time.int, box.dt)
stopifnot(nrow(over.dt)==nrow(time.int))
median.dt <- over.dt[, list(
  problems=.N,
  q25.intervals=quantile(mean.intervals, 0.25),
  max.intervals=max(max.intervals),
  median.intervals=median(mean.intervals),
  q75.intervals=quantile(mean.intervals, 0.75),
  q25.seconds=quantile(seconds, 0.25),
  median.seconds=median(seconds),
  q75.seconds=quantile(seconds, 0.75)
  ), by=list(algo.type, mid)]
median.dt[, data := 10^mid]

range(time.int$data)
range(time.int$max.intervals)
max.dt <- median.dt[, .SD[which.max(max.intervals)], by=list(algo.type)]

lab.df <- data.frame(
  y="seconds",
  seconds=c(1, 60, 60*2),
  label=c("1 second", "1 minute", "2 minutes"))
gg <- ggplot()+
  ggtitle("Infinite limits increase intervals and runtime by a constant")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(y ~ ., scales="free")+
  geom_hline(aes(yintercept=seconds),
             data=lab.df,
             color="grey")+
  geom_text(aes(100, seconds, label=label),
            data=lab.df,
            ##size=3,
            color="grey",
            hjust=0,
            vjust=1.5)+
  geom_line(aes(
    data, max.intervals, color=algo.type),
            data=data.table(median.dt, y="intervals"))+
  geom_text(aes(
    data, max.intervals, color=algo.type, label=paste(" ", max.intervals, "intervals")),
            vjust=0,
            hjust=0,
             data=data.table(max.dt, y="intervals"))+
  geom_point(aes(
    data, max.intervals, color=algo.type),
             data=data.table(max.dt, y="intervals"))+
  geom_ribbon(aes(
    data, ymin=q25.intervals, ymax=q75.intervals, fill=algo.type),
              alpha=0.5,
            data=data.table(median.dt, y="intervals"))+
  geom_line(aes(
    data, median.intervals, color=algo.type),
            data=data.table(median.dt, y="intervals"))+
  scale_y_log10("")+
  scale_x_log10("data", breaks=c(87, 1000, 10000, 100000, 263169), limits=c(87, 3e5))+
  geom_ribbon(aes(
    data, ymin=q25.seconds, ymax=q75.seconds, fill=algo.type),
              alpha=0.5,
            data=data.table(median.dt, y="seconds"))+
  geom_line(aes(
    data, median.seconds, color=algo.type),
            data=data.table(median.dt, y="seconds"))
gg.lab <- direct.label(gg)
pdf("figure-PDPAInf-intervals.pdf", 10, 10)
print(gg.lab)
dev.off()

