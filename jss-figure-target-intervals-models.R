library(data.table)
library(ggplot2)
library(directlabels)

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)]
prob.dir.vec <- c(
  "Most bases"=bench.models[which.max(bases), prob.dir],
  "Most compressed data"=bench.models[which.max(bedGraph.lines), prob.dir],
  "Largest mean intervals"=bench.models[which.max(mean.intervals), prob.dir],
  "Largest max intervals"=bench.models[which.max(max.intervals), prob.dir],
  "Most megabytes stored"=bench.models[which.max(megabytes), prob.dir],
  "Most computation time"=bench.models[which.max(seconds), prob.dir])
t(t(prob.dir.vec))
prob.label.dt <- data.table(
  prob.label=names(prob.dir.vec),
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

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(
    log(penalty), value, color=stat),
    shape=1,
    data=some.probs.both)+
  scale_y_log10("")+
  scale_x_continuous("log(penalty)")+
  facet_grid(variable ~ prob.label, scales="free_y")

bench.models[, gigabytes := megabytes/1024]
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
  measure.vars=c("gigabytes", "minutes"))
bench.models.tall.ranges <- bench.models.tall[, list(
  min.value=min(value),
  max.value=max(value)
), by=list(bedGraph.lines, prob.dir, variable)]

## this data point min value stands out -- probably an optimization
## error.
bench.models.tall.ranges[6 < log10(bedGraph.lines) & min.value < 0.1]

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(variable ~ ., scales="free")+
  geom_segment(aes(
    log10(bedGraph.lines), min.value,
    xend=log10(bedGraph.lines), yend=max.value),
    data=bench.models.tall.ranges[!(6 < log10(bedGraph.lines) & min.value < 0.1)])+
  scale_y_log10("")+
  scale_x_continuous("log10(number of data after compression = lines in bedGraph file)")
pdf("jss-figure-target-intervals-models.pdf")
print(gg)
dev.off()
