source("jss-packages.R")

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)][log(bedGraph.lines) < penalty & penalty < bedGraph.lines & 1000 < bedGraph.lines]
bench.models[, minutes := seconds/60]
bench.models[, hours := minutes/60]
bench.models[, gigabytes := megabytes/1024]

min.err.ranges <- bench.models[, .SD[errors==min(errors), list(
  min.penalty=min(penalty),
  max.penalty=max(penalty),
  min.peaks=min(peaks),
  max.peaks=max(peaks)
)], by=list(bedGraph.lines, prob.dir)][order(bedGraph.lines)]

no.zero <- min.err.ranges[0 < min.peaks]
seg.dt <- no.zero[min.peaks < max.peaks]
point.dt <- no.zero[min.peaks == max.peaks]
ggplot()+
  geom_segment(aes(
    bedGraph.lines, min.peaks,
    xend=bedGraph.lines, yend=max.peaks),
    data=seg.dt)+
  geom_point(aes(
    bedGraph.lines, min.peaks),
    color="red",
    data=point.dt)+
  scale_x_log10()+
  scale_y_log10()

penalty.ranges <- bench.models[0 < penalty & penalty < Inf, list(
  min.penalty=min(penalty),
  max.penalty=max(penalty)
), by=list(bedGraph.lines, prob.dir)]

ref.dt <- data.table(N=10^seq(2, 7, l=100))
fun.list <- list(
  log=log,
  linear=identity,
  sqrt=sqrt)
for(fun.name in names(fun.list)){
  fun <- fun.list[[fun.name]]
  ref.dt[[fun.name]] <- fun(ref.dt$N)
}
ref.tall <- melt(ref.dt, id.vars="N")

leg <- ggplot()+
  geom_segment(aes(
    bedGraph.lines, min.penalty,
    xend=bedGraph.lines, yend=max.penalty),
    data=penalty.ranges)+
  geom_segment(aes(
    bedGraph.lines, min.penalty,
    xend=bedGraph.lines, yend=max.penalty),
    color="red",
    data=min.err.ranges)+
  geom_line(aes(
    N, value, color=variable),
    data=ref.tall)+
  scale_x_log10()+
  scale_y_log10()
direct.label(leg, "last.polygons")

## Plot only models with penalty between log N and N.
log.linear.models <- bench.models[log(bedGraph.lines) < penalty & penalty < bedGraph.lines]
log10.range <- log.linear.models[, range(bedGraph.lines)]
ref.dt <- data.table(N.data=10^seq(log10.range[1], log10.range[2], l=100))
fun.list <- list(
  "log(N)"=log,
  "N"=identity)
for(fun.name in names(fun.list)){
  fun <- fun.list[[fun.name]]
  ref.dt[[fun.name]] <- fun(ref.dt$N.data)
}
ref.tall <- melt(ref.dt, id.vars="N.data")
leg <- ggplot()+
  geom_point(aes(
    bedGraph.lines, penalty),
    shape=1,
    data=log.linear.models)+
  geom_line(aes(
    N.data, value, color=variable),
    data=ref.tall)+
  scale_x_log10("N = number of data to segment (log scale)")+
  scale_y_log10()
dl <- direct.label(leg, "last.polygons")

## Plot a dot at the middle number of peaks.
min.err.ranges[, mid.peaks := (min.peaks+max.peaks)/2]
ggplot()+
  geom_point(aes(
    bedGraph.lines, mid.peaks),
    data=min.err.ranges)+
  scale_x_log10()+
  scale_y_log10()

## do it by experiment type.
min.err.ranges[, experiment := sub("_.*", "", prob.dir)]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("experiment")+
  geom_point(aes(
    bedGraph.lines, mid.peaks),
    data=min.err.ranges)+
  scale_x_log10()+
  scale_y_log10()

## Make boxes with the median and quartiles of the number of peaks, by
## experiment
log10.range <- log10(range(min.err.ranges$bedGraph.lines))
box.dt <- data.table(
  box.mid=10^seq(log10.range[1], log10.range[2], l=8))
(diff.vec <- diff(log10(box.dt$box.mid)))
box.w <- diff.vec[1]/2
box.dt[, box.min := 10^(log10(box.mid)-box.w) ]
box.dt[, box.max := 10^(log10(box.mid)+box.w) ]
box.models <- box.dt[min.err.ranges, on=list(
  box.min < bedGraph.lines,
  box.max > bedGraph.lines)]
stopifnot(nrow(box.models) == nrow(min.err.ranges))
box.models.stats <- box.models[, list(
  median=median(mid.peaks),
  q75=quantile(mid.peaks, 0.75),
  q25=quantile(mid.peaks, 0.25),
  min=min(mid.peaks),
  max=max(mid.peaks),
  models=.N
), by=list(box.mid, experiment)][order(box.mid)]
ref.dt <- data.table(N.data=10^seq(log10.range[1], log10.range[2], l=100))
log10.range <- log10(range(box.models.stats$box.mid))
N.data <- 10^seq(log10.range[1], log10.range[2], l=100)
fun.list <- list(
  "O(N)"=identity,
  "O(log N)"=log,
  "loglog(N)"=function(x)log(log(x)),
  "O(sqrt N)"=sqrt)
ref.line.list <- list(
  ##OP=list(y=9, lines=c("log(N)", "sqrt(N)", "loglog(N)")),
  SN=list(y=13.5, lines=c("O(N)", "O(log N)", "O(sqrt N)")))
ref.tall.list <- list()
for(ref.name in names(ref.line.list)){
  ref.info <- ref.line.list[[ref.name]]
  for(fun.name in ref.info$lines){
    fun <- fun.list[[fun.name]]
    first.y <- fun(min(N.data))
    ref.tall.list[[paste(fun.name, ref.name)]] <- data.table(
      N.data,
      ref.name,
      fun.name,
      value=fun(N.data)/first.y*ref.info$y)
  }
}
ref.tall <- do.call(rbind, ref.tall.list)
ref.color <- 'red'
gg.bands <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("experiment")+
  geom_line(aes(
    N.data, value, group=paste(ref.name, fun.name)),
    color=ref.color,
    data=ref.tall)+
  geom_text(aes(
    N.data, value, label=fun.name),
    color=ref.color,
    hjust=0,
    data=ref.tall[N.data==max(N.data)])+
  geom_ribbon(aes(
    box.mid, ymin=q25, ymax=q75),
    data=box.models.stats,
    alpha=0.5)+
  geom_line(aes(
    box.mid, median),
    data=box.models.stats)+
  scale_x_log10(
    "N = number of data to segment (log scale)",
    limits=c(NA, 10^7.5),
    labels=paste)+
  scale_y_log10("Peaks in models with min label error\n(log scale)")
gg.points <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("experiment")+
  geom_point(aes(
    bedGraph.lines, mid.peaks),
    data=min.err.ranges)+
  scale_x_log10()+
  scale_y_log10()+
  geom_line(aes(
    N.data, value, group=paste(ref.name, fun.name)),
    color=ref.color,
    data=ref.tall)+
  geom_text(aes(
    N.data, value, label=fun.name),
    color=ref.color,
    hjust=0,
    data=ref.tall[N.data==max(N.data)])

## Make boxes with the median and quartiles of the number of peaks.
log10.range <- log10(range(min.err.ranges$bedGraph.lines))
box.dt <- data.table(
  box.mid=10^seq(log10.range[1], log10.range[2], l=8))
(diff.vec <- diff(log10(box.dt$box.mid)))
box.w <- diff.vec[1]/2
box.dt[, box.min := 10^(log10(box.mid)-box.w) ]
box.dt[, box.max := 10^(log10(box.mid)+box.w) ]
box.models <- box.dt[min.err.ranges, on=list(
  box.min < bedGraph.lines,
  box.max > bedGraph.lines)]
stopifnot(nrow(box.models) == nrow(min.err.ranges))
box.models.stats <- box.models[, list(
  median=median(mid.peaks),
  q75=quantile(mid.peaks, 0.75),
  q25=quantile(mid.peaks, 0.25),
  min=min(mid.peaks),
  max=max(mid.peaks),
  models=.N
), by=list(box.mid)][order(box.mid)]
ref.dt <- data.table(N.data=10^seq(log10.range[1], log10.range[2], l=100))
log10.range <- log10(range(box.models.stats$box.mid))
N.data <- 10^seq(log10.range[1], log10.range[2], l=100)
fun.list <- list(
  "O(N)"=identity,
  "O(log N)"=log,
  "loglog(N)"=function(x)log(log(x)),
  "O(sqrt N)"=sqrt)
ref.line.list <- list(
  ##OP=list(y=9, lines=c("log(N)", "sqrt(N)", "loglog(N)")),
  SN=list(y=13.5, lines=c("O(N)", "O(log N)", "O(sqrt N)")))
ref.tall.list <- list()
for(ref.name in names(ref.line.list)){
  ref.info <- ref.line.list[[ref.name]]
  for(fun.name in ref.info$lines){
    fun <- fun.list[[fun.name]]
    first.y <- fun(min(N.data))
    ref.tall.list[[paste(fun.name, ref.name)]] <- data.table(
      N.data,
      ref.name,
      fun.name,
      value=fun(N.data)/first.y*ref.info$y)
  }
}
ref.tall <- do.call(rbind, ref.tall.list)
ref.color <- 'red'
gg <- ggplot()+
  theme_bw()+
  geom_line(aes(
    N.data, value, group=paste(ref.name, fun.name)),
    color=ref.color,
    data=ref.tall)+
  geom_text(aes(
    N.data, value, label=fun.name),
    color=ref.color,
    hjust=0,
    data=ref.tall[N.data==max(N.data)])+
  geom_ribbon(aes(
    box.mid, ymin=q25, ymax=q75),
    data=box.models.stats,
    alpha=0.5)+
  geom_line(aes(
    box.mid, median),
    data=box.models.stats)+
  scale_x_log10(
    "N = number of data to segment (log scale)",
    limits=c(NA, 10^7.5),
    labels=paste)+
  scale_y_log10("Peaks in models with min label error\n(log scale)")

pdf("jss-figure-data-peaks.pdf", 7, 3)
print(gg)
dev.off()

## Make boxes with the median and quartiles of the number of peaks. (TIKZ)
log10.range <- log10(range(min.err.ranges$bedGraph.lines))
box.dt <- data.table(
  box.mid=10^seq(log10.range[1], log10.range[2], l=8))
(diff.vec <- diff(log10(box.dt$box.mid)))
box.w <- diff.vec[1]/2
box.dt[, box.min := 10^(log10(box.mid)-box.w) ]
box.dt[, box.max := 10^(log10(box.mid)+box.w) ]
box.models <- box.dt[min.err.ranges, on=list(
  box.min < bedGraph.lines,
  box.max > bedGraph.lines)]
stopifnot(nrow(box.models) == nrow(min.err.ranges))
box.models.stats <- box.models[, list(
  median=median(mid.peaks),
  q75=quantile(mid.peaks, 0.75),
  q25=quantile(mid.peaks, 0.25),
  min=min(mid.peaks),
  max=max(mid.peaks),
  models=.N
), by=list(box.mid)][order(box.mid)]
ref.dt <- data.table(N.data=10^seq(log10.range[1], log10.range[2], l=100))
log10.range <- log10(range(box.models.stats$box.mid))
N.data <- 10^seq(log10.range[1], log10.range[2], l=100)
fun.list <- list(
  "$O(N)$"=identity,
  "$O(\\log N)$"=log,
  "$O(\\sqrt N)$"=sqrt)
ref.line.list <- list(
  SN=list(y=13.5, lines=names(fun.list)))
ref.tall.list <- list()
for(ref.name in names(ref.line.list)){
  ref.info <- ref.line.list[[ref.name]]
  for(fun.name in ref.info$lines){
    fun <- fun.list[[fun.name]]
    first.y <- fun(min(N.data))
    ref.tall.list[[paste(fun.name, ref.name)]] <- data.table(
      N.data,
      ref.name,
      fun.name,
      value=fun(N.data)/first.y*ref.info$y)
  }
}
ref.tall <- do.call(rbind, ref.tall.list)
ref.color <- 'red'
gg <- ggplot()+
  theme_bw()+
  geom_line(aes(
    N.data, value, group=paste(ref.name, fun.name)),
    color=ref.color,
    data=ref.tall)+
  geom_text(aes(
    N.data, value, label=fun.name),
    color=ref.color,
    hjust=0,
    data=ref.tall[N.data==max(N.data)])+
  geom_ribbon(aes(
    box.mid, ymin=q25, ymax=q75),
    data=box.models.stats,
    alpha=0.5)+
  geom_line(aes(
    box.mid, median),
    data=box.models.stats)+
  scale_x_log10(
    "N = number of data to segment (log scale)",
    limits=c(NA, 10^7.5),
    labels=paste)+
  scale_y_log10("Peaks in models with min label error\n(log scale)")

tikz("jss-figure-data-peaks.tex", 6, 3)
print(gg)
dev.off()
