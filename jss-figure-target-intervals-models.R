library(data.table)
library(ggplot2)
library(directlabels)

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)]
bench.models[, gigabytes := megabytes/1024]
prob.dir.vec <- c(
  "Most bases"=bench.models[which.max(bases), prob.dir],
  "Most weighted data"=bench.models[which.max(bedGraph.lines), prob.dir],
  "Largest mean intervals"=bench.models[which.max(mean.intervals), prob.dir],
  "Largest max intervals"=bench.models[which.max(max.intervals), prob.dir],
  "Most megabytes stored"=bench.models[which.max(megabytes), prob.dir],
  "Most computation time"=bench.models[which.max(seconds), prob.dir])
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
  measure.vars=c("minutes", "gigabytes"))

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
  max.gigabytes=max(gigabytes)
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
  measure.vars=c("gigabytes", "minutes", "mean.intervals"))
bench.models.tall[, var := ifelse(
  variable=="mean.intervals", "intervals", paste(variable))]
bench.models.tall[, stat := ifelse(
  variable=="mean.intervals", "mean", "total")]
bench.models.tall.ranges <- bench.models.tall[, list(
  min.value=min(value),
  max.value=max(value)
), by=list(bedGraph.lines, prob.dir, variable, var, stat)]
bench.models.max.intervals <- bench.models[, list(
  max.intervals=max(max.intervals),
  stat="max",
  var="intervals"
), by=list(bedGraph.lines, prob.dir)]

## this data point min value stands out -- probably an optimization
## error.
bench.models.tall.ranges[6 < log10(bedGraph.lines) & min.value < 0.1]
show.segments <- bench.models.tall.ranges[!(6 < log10(bedGraph.lines) & min.value < 0.1)]
show.max <- rbind(show.segments[, {
  .SD[which.max(max.value)]
}, by=list(var, stat)][, list(
  var, stat, bedGraph.lines, prob.dir, max.value
)], bench.models.max.intervals[which.max(max.intervals), list(
  var, stat, bedGraph.lines, prob.dir, max.value=max.intervals
  )])
blank.dt <- data.table(
  var=c("gigabytes", "minutes", "intervals"),
  y=c(10^2.9, 10^3, 10^3.1))
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

pdf("jss-figure-target-intervals-models.pdf", 7, 3)
print(gg)
dev.off()

