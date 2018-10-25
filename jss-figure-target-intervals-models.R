source("jss-packages.R")

bench.models <- fread("jss.bench.models.csv")
bench.models[, gigabytes := megabytes/1024]
prob.dir.vec <- c(
  ##"Most bases"=bench.models[which.max(bases), prob.dir],
  "Most weighted data"=bench.models[which.max(bedGraph.lines), prob.dir],
  "Largest mean intervals"=bench.models[which.max(mean.intervals), prob.dir],
  "Largest max intervals"=bench.models[which.max(max.intervals), prob.dir],
  ##"Most computation time"=bench.models[which.max(seconds), prob.dir],
  "Most megabytes stored"=bench.models[which.max(megabytes), prob.dir])
t(t(prob.dir.vec))
prob.label.dt <- data.table(
  prob.label=gsub(" ", "\n", names(prob.dir.vec)),
  prob.dir=prob.dir.vec)
one.prob <- bench.models[prob.dir==prob.dir.vec[["Largest max intervals"]] ]

ggplot()+
  geom_point(aes(
    log(penalty), log10(megabytes)),
    shape=1,
    data=one.prob)

one.prob.intervals <- melt(
  one.prob[0 < megabytes], measure.vars=c("mean.intervals", "max.intervals"))
one.prob.intervals[, stat := sub(".intervals", "", variable)]

leg <- ggplot()+
  geom_point(aes(
    log(penalty), log10(value), color=stat),
    shape=1,
    data=one.prob.intervals)+
  ylab("log10(Intervals stored in the functional cost)")+
  scale_x_continuous("log(penalty)", limits=c(NA, 10.5))
direct.label(leg, "last.polygons")

some.probs <- bench.models[prob.label.dt, on=list(prob.dir)]
some.probs.intervals <- melt(
  some.probs[0 < megabytes],
  measure.vars=c("mean.intervals", "max.intervals"))
some.probs.intervals[, stat := sub(".intervals", "", variable)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(
    log(penalty), value, color=stat),
    shape=1,
    data=some.probs.intervals)+
  scale_y_log10("Intervals stored in the functional cost")+
  scale_x_continuous("log(penalty)")+
  facet_grid(. ~ prob.label)

some.probs[, minutes := seconds/60]
some.probs[, hours := minutes/60]
some.probs[, gigabytes := megabytes/1024]
some.probs.other <- melt(
  some.probs,
  measure.vars=c(
    ##"minutes",
    "gigabytes"))

some.probs.both <- rbind(
  some.probs.intervals[, list(
    prob.label, penalty, variable="intervals", value, stat)],
  some.probs.other[, list(
    prob.label, penalty, variable, value, stat="total")])

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(
    log(penalty), value, color=stat),
    shape=1,
    data=some.probs.both[penalty<Inf])+
  scale_y_log10("")+
  scale_x_continuous("log(penalty)")+
  facet_grid(variable ~ prob.label, scales="free_y")
##print(gg)
pdf("jss-figure-target-intervals-models-penalty.pdf", 7, 3)
print(gg)
dev.off()

gigabyte.ranges <- bench.models[0 < gigabytes, list(
  min.gigabytes=min(gigabytes),
  max.gigabytes=max(gigabytes),
  models=.N
  ), by=list(bedGraph.lines, prob.dir)]

ggplot()+
  geom_segment(aes(
    log10(bedGraph.lines), min.gigabytes,
    xend=log10(bedGraph.lines), yend=max.gigabytes),
    data=gigabyte.ranges)+
  scale_y_log10()

bench.models[, minutes := seconds/60]
bench.models.tall <- melt(
  bench.models[0 < gigabytes],
  measure.vars=c(
    "gigabytes",
    "mean.intervals",
    "minutes"
  ))
bench.models.tall[, var := ifelse(
  variable=="mean.intervals", "intervals", paste(variable))]
bench.models.tall[, stat := ifelse(
  variable=="mean.intervals", "mean", "total")]
bench.models.tall.ranges <- bench.models.tall[, list(
  min.value=min(value),
  max.value=max(value),
  mean.value=mean(value)
), by=list(bedGraph.lines, prob.dir, variable, var, stat)]
bench.models.tall.ranges[, mid.value := (min.value+max.value)/2]

ggplot()+
  geom_point(aes(
    mid.value, mean.value),
    data=bench.models.tall.ranges)+
  coord_equal()+
  theme_bw()+
  geom_abline(intercept=0, slope=1, color="grey")

bench.models.max.intervals <- bench.models[, list(
  max.intervals=max(max.intervals),
  stat="max",
  var="intervals"
), by=list(bedGraph.lines, prob.dir)]

## this data point min value stands out -- probably an optimization
## error.
bench.models.tall.ranges[6 < log10(bedGraph.lines) & min.value < 0.1]
show.segments <- bench.models.tall.ranges[!(6 < log10(bedGraph.lines) & min.value < 0.1)]
show.max <- rbind(show.segments[variable != "mean.intervals", {
  .SD[which.max(max.value)]
}, by=list(var, stat)][, list(
  var, stat, bedGraph.lines, prob.dir, max.value
)], bench.models.max.intervals[which.max(max.intervals), list(
  var, stat, bedGraph.lines, prob.dir, max.value=max.intervals
  )])
blank.dt <- data.table(
  var=c("gigabytes", "minutes", "intervals"),
  y=c(10^3.2, 10^3.4, 10^3.2))
max.color <- "black"
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(var ~ ., scales="free")+
  geom_blank(aes(
    3, y),
    data=blank.dt)+
  geom_hline(aes(
    yintercept=yint),
    color="grey50",
    data=data.table(yint=1, var="gigabytes"))+
  geom_text(aes(
    x, y, label=label),
    color="grey50", 
    data=data.table(x=3, y=1, label="1 gigabyte", var="gigabytes"),
    vjust=-0.5)+
  geom_line(aes(
    log10(bedGraph.lines), max.intervals, color=stat),
    data=bench.models.max.intervals)+
  geom_segment(aes(
    log10(bedGraph.lines), min.value,
    color=stat,
    xend=log10(bedGraph.lines), yend=max.value),
    data=show.segments)+
  geom_point(aes(
    log10(bedGraph.lines), max.value),
    shape=1,
    color=max.color,
    data=show.max)+
  geom_text(aes(
    log10(bedGraph.lines)-0.05, max.value,
    label=paste(
      format(bedGraph.lines, big.mark=","),
      "data,",
      round(max.value), var)),
    color=max.color,
    hjust=1,
    vjust=0,
    data=show.max)+
  scale_y_log10("")+
  scale_x_continuous("log10(number of data after compression = lines in bedGraph file)")
print(gg)

log10.range <- log10(range(show.segments$bedGraph.lines))
box.dt <- data.table(
  box.mid=10^seq(log10.range[1], log10.range[2], l=7))
(diff.vec <- diff(log10(box.dt$box.mid)))
box.w <- diff.vec[1]/2
box.dt[, box.min := 10^(log10(box.mid)-box.w) ]
box.dt[, box.max := 10^(log10(box.mid)+box.w) ]
box.dt[, box.i := 1:.N]
box.models <- box.dt[bench.models, on=list(
  box.min < bedGraph.lines,
  box.max > bedGraph.lines)]
stopifnot(nrow(box.models) == nrow(bench.models))
box.segments <- box.dt[show.segments, on=list(
  box.min < bedGraph.lines,
  box.max > bedGraph.lines)]
stopifnot(nrow(box.segments) == nrow(show.segments))
box.segments.stats <- box.segments[, list(
  median=median(mid.value),
  q95=quantile(mid.value, 0.95),
  q05=quantile(mid.value, 0.05),
  min=min(mid.value),
  max=max(mid.value),
  models=.N
), by=list(var, stat, box.mid)]
box.max <- box.dt[bench.models.max.intervals, on=list(
  box.min < bedGraph.lines,
  box.max > bedGraph.lines)]
stopifnot(nrow(box.max) == nrow(bench.models.max.intervals))
box.max.stats <- box.max[, list(
  max.intervals=max(max.intervals)
), by=list(var, stat, box.mid)]

box.segments.stats[var=="intervals" & stat=="mean" & box.mid==max(box.mid)]

##dput(RColorBrewer::brewer.pal(Inf, "Set1"))
stat.colors <- c(
  mean="#E41A1C",
  max="#377EB8",
  total="#4DAF4A",
  "#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999")
line.size <- 1
line.dt <- rbind(box.max.stats[, data.table(
  var, stat, box.mid, line.value=max.intervals
)], box.segments.stats[, data.table(
  var, stat, box.mid, line.value=median
  )])
(show.point <- show.segments[bedGraph.lines ==max(bedGraph.lines) & variable=="mean.intervals"])
hline.dt <- rbind(
  data.table(x=10^3.5, y=1, vjust=-0.5, label="1 gigabyte", var="gigabytes"),
  data.table(x=10^3.5, y=60, vjust=1.25, label="1 hour", var="minutes"))

ref.dt <- data.table(N.data=10^seq(log10.range[1], log10.range[2], l=100))
fun.list <- list(
  "log(N)"=log,
  "N log(N)"=function(x)x*log(x),
  "N"=identity)
one.var <- "gigabytes"
one.line <- line.dt[var==one.var]
first.row <- one.line[order(box.mid)][1]
for(fun.name in names(fun.list)){
  fun <- fun.list[[fun.name]]
  first.y <- fun(first.row$box.mid)
  ref.dt[[fun.name]] <- fun(ref.dt$N.data)/first.y*first.row$line.value
}
ref.tall <- melt(ref.dt, id.vars="N.data")
leg <- ggplot()+
  geom_line(aes(
    log10(N.data), log10(value), color=variable),
    size=2,
    data=ref.tall)+
  geom_ribbon(aes(
    log10(box.mid), ymin=log10(q05), ymax=log10(q95)),
    alpha=0.5,
    data=box.segments.stats[var==one.var])+
  geom_line(aes(
    log10(box.mid), log10(line.value)),
    data=one.line)+
  scale_x_continuous(
    "log10(N = number of data to segment)",
    limits=c(NA, 7.5)
  )
dl <- direct.label(leg, "last.polygons")
pdf("jss-figure-target-intervals-models-NlogN.pdf")
print(dl)
dev.off()


ref.dt <- data.table(N.data=10^seq(log10.range[1], log10.range[2], l=100))
fun.list <- list(
  "log(N)"=log,
  "loglog(N)"=function(x)log(log(x)),
  "sqrt(log(N))"=function(x)sqrt(log(x)),
  "sqrt(N)"=sqrt)
one.var <- "intervals"
one.line <- line.dt[var==one.var & stat=="mean"]
first.row <- one.line[order(box.mid)][1]
for(fun.name in names(fun.list)){
  fun <- fun.list[[fun.name]]
  first.y <- fun(first.row$box.mid)
  ref.dt[[fun.name]] <- fun(ref.dt$N.data)/first.y*first.row$line.value
}
ref.tall <- melt(ref.dt, id.vars="N.data")
leg <- ggplot()+
  geom_line(aes(
    log10(N.data), log10(value), color=variable),
    size=2,
    data=ref.tall)+
  geom_ribbon(aes(
    log10(box.mid), ymin=log10(q05), ymax=log10(q95)),
    alpha=0.5,
    data=box.segments.stats[var==one.var])+
  geom_line(aes(
    log10(box.mid), log10(line.value)),
    data=one.line)+
  ylab("log10(mean number of intervals)")+
  scale_x_continuous(
    "log10(N = number of data to segment)",
    limits=c(NA, 7.75)
  )
dl <- direct.label(leg, "last.polygons")
pdf("jss-figure-target-intervals-models-logN.pdf")
print(dl)
dev.off()

leg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(var ~ ., scales="free")+
  geom_blank(aes(
    10^3.5, y),
    data=blank.dt)+
  geom_hline(aes(
    yintercept=y),
    color="grey50",
    data=hline.dt)+
  geom_text(aes(
    x, y, label=label),
    vjust=-0.5,
    color="grey50", 
    data=hline.dt)+
  geom_ribbon(aes(
    box.mid, ymin=q05, ymax=q95, fill=stat),
    alpha=0.5,
    data=box.segments.stats)+
  geom_line(aes(
    box.mid, line.value, color=stat),
    size=line.size,
    data=line.dt)+
  geom_point(aes(
    bedGraph.lines, max.value),
    shape=1,
    color=max.color,
    data=show.max)+
  geom_text(aes(
    bedGraph.lines, max.value,
    label=paste(
      format(bedGraph.lines, big.mark=","),
      "data,",
      round(max.value), var)),
    color=max.color,
    hjust=1,
    vjust=-0.5,
    data=show.max)+
  ## show median.
  geom_point(aes(
    bedGraph.lines, mean.value),
    shape=1,
    color=max.color,
    data=show.point)+
  geom_text(aes(
    bedGraph.lines, mean.value,
    label=paste(
      format(bedGraph.lines, big.mark=","),
      "data,",
      round(mean.value), var)),
    color=max.color,
    hjust=1,
    vjust=-1,
    data=show.point)+
  scale_y_log10("(log scales)", labels=function(chr){
    paste(chr)
  })+
  scale_color_manual(values=stat.colors)+
  scale_fill_manual(values=stat.colors, guide=FALSE)+
  scale_x_log10(
    "N = number of data to segment (log scale)",
    limits=c(NA, 10^7.25),
    labels=paste
  )
gg <- direct.label(leg, "last.polygons")
pdf("jss-figure-target-intervals-models-all.pdf", 3, 3)
print(gg)
dev.off()

## Left and right figures. Right plot disk usage and time.
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(var ~ ., scales="free")+
  geom_blank(aes(
    10^4.5, y),
    data=blank.dt[var!="intervals"])+
  geom_hline(aes(
    yintercept=y),
    color="grey50",
    data=hline.dt)+
  geom_text(aes(
    10^4, y, label=label, vjust=vjust),
    color="grey50", 
    data=hline.dt)+
  geom_ribbon(aes(
    box.mid, ymin=q05, ymax=q95),
    alpha=0.5,
    data=box.segments.stats[var!="intervals"])+
  geom_line(aes(
    box.mid, line.value),
    size=line.size,
    data=line.dt[var!="intervals"])+
  geom_point(aes(
    bedGraph.lines, max.value),
    shape=1,
    color=max.color,
    data=show.max[var!="intervals"])+
  geom_text(aes(
    bedGraph.lines, max.value,
    label=paste(
      format(bedGraph.lines, big.mark=","),
      "data,",
      round(max.value), var)),
    color=max.color,
    hjust=1,
    vjust=-0.5,
    data=show.max[var!="intervals"])+
  scale_y_log10("Computational requirements\n(log scales)", labels=function(chr){
    paste(chr)
  })+
  scale_x_log10(
    "N = data to segment (log scale)",
    labels=paste
  )
pdf("jss-figure-target-intervals-models-computation.pdf", 3.5, 3)
print(gg)
dev.off()

## Left plot intervals.
leg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_ribbon(aes(
    box.mid, ymin=q05, ymax=q95, fill=stat),
    alpha=0.5,
    data=box.segments.stats[var=="intervals"])+
  geom_line(aes(
    box.mid, line.value, color=stat),
    size=line.size,
    data=line.dt[var=="intervals"])+
  geom_point(aes(
    bedGraph.lines, max.value),
    shape=1,
    color=max.color,
    data=show.max[var=="intervals"])+
  geom_text(aes(
    1e7, max.value,
    label=paste(
      format(bedGraph.lines, big.mark=","),
      "data,",
      round(max.value), var)),
    color=max.color,
    hjust=1,
    vjust=-0.5,
    data=show.max[var=="intervals"])+
  ## show median.
  geom_point(aes(
    bedGraph.lines, mean.value),
    shape=1,
    color=max.color,
    data=show.point[var=="intervals"])+
  geom_text(aes(
    3e7, mean.value,
    label=paste(
      format(bedGraph.lines, big.mark=","),
      "data,",
      round(mean.value), var)),
    color=max.color,
    hjust=1,
    vjust=-1.75,
    data=show.point)+
  scale_y_log10(
    "Intervals/candidate changepoints\n(log scale)",
    limits=c(NA, 1000),
    labels=function(chr){
      paste(chr)
    })+
  scale_color_manual(values=stat.colors)+
  scale_fill_manual(values=stat.colors, guide=FALSE)+
  scale_x_log10(
    "N = data to segment (log scale)",
    limits=c(NA, 10^7.75),
    labels=paste
  )
dl <- direct.label(leg, list("last.qp", dl.trans(x=x+0.1)))

gg <- leg+
  geom_line(aes(
    N.data, value, group=variable),
    color="grey",
    data=ref.tall)+
  geom_text(aes(
    N.data, value, label=variable),
    hjust=0,
    color="grey",
    data=ref.tall[N.data==max(N.data)])+

pdf("jss-figure-target-intervals-models.pdf", 3.5, 3)
print(dl)
dev.off()
