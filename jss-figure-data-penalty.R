library(data.table)
library(directlabels)
library(ggplot2)

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)]
bench.models[, minutes := seconds/60]
bench.models[, hours := minutes/60]
bench.models[, gigabytes := megabytes/1024]

min.err.ranges <- bench.models[, .SD[errors==min(errors), list(
  min.penalty=min(penalty),
  max.penalty=max(penalty),
  min.peaks=min(peaks),
  max.peaks=max(peaks)
)], by=list(bedGraph.lines, prob.dir)]

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
  ## geom_segment(aes(
  ##   bedGraph.lines, min.penalty,
  ##   xend=bedGraph.lines, yend=max.penalty),
  ##   color="red",
  ##   data=min.err.ranges)+
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
  scale_x_log10("N = number of weighted data to segment")+
  scale_y_log10()
dl <- direct.label(leg, "last.polygons")
png("jss-figure-data-penalty.png", 7, 3, units="in", res=200)
print(dl)
dev.off()
