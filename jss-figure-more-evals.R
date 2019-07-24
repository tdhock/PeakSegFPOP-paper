source("jss-packages.R")

jss.more.evals <- readRDS("jss.more.evals.rds")

evals.dt <- jss.more.evals[, .(
  penalties=.N,
  max.peaks=max(others.peaks)
), by=list(bedGraph.lines, peaks.arg, loss.peaks)]

pen0 <- unique(
  jss.more.evals[others.penalty==0, .(bedGraph.lines, others.peaks)])
ggplot()+
  theme_bw()+
  coord_equal()+
  geom_abline(aes(
    slope=slope, intercept=intercept),
    data=data.table(slope=0.5, intercept=0),
    color="grey")+
  geom_point(aes(
    bedGraph.lines, max.peaks),
    data=evals.dt)+
  geom_point(aes(
    bedGraph.lines, others.peaks),
    color="red",
    shape=1,
    data=pen0)

ggplot()+
  geom_point(aes(
    peaks.arg, penalties),
    data=evals.dt)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap("bedGraph.lines")

ggplot()+
  scale_color_gradient(low="white", high="red")+
  geom_line(aes(
    peaks.arg, penalties, group=max.peaks, color=log10(max.peaks)),
    size=2,
    data=evals.dt)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))

ggplot()+
  scale_color_gradient(low="white", high="red")+
  geom_line(aes(
    max.peaks, penalties, group=peaks.arg, color=log10(peaks.arg)),
    size=2,
    data=evals.dt)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))

ggplot()+
  geom_point(aes(
    bedGraph.lines, penalties),
    data=evals.dt)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap("peaks.arg")

ggplot()+
  geom_point(aes(
    max.peaks, penalties),
    data=evals.dt)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap("peaks.arg")

evals.dt[, desired.peaks := peaks.arg]
problem.counts <- evals.dt[, .(
  problems=.N
), by=list(desired.peaks)]
max.counts <- problem.counts[problems==max(problems)][desired.peaks < 2000]
max.evals <- evals.dt[max.counts, on=list(desired.peaks)]
ggplot()+
  geom_point(aes(
    max.peaks, penalties),
    data=max.evals)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap("peaks.arg")

evals.tall <- melt(
  max.evals,
  measure.vars=c("max.peaks", "desired.peaks"))
stats.dt <- evals.tall[, .(
  mean.penalties=mean(penalties),
  sd.penalties=sd(penalties),
  penalties=.N
  ), by=list(variable, value)]
ggplot()+
  geom_ribbon(aes(
    value,
    ymax=mean.penalties+sd.penalties,
    ymin=mean.penalties-sd.penalties),
    data=stats.dt,
    alpha=0.5)+
  geom_line(aes(
    value, mean.penalties),
    data=stats.dt)+
  scale_y_log10(
  "GFPOP calls required")+
  scale_x_log10("")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ variable, scales="free")

gg <- ggplot()+
  geom_point(aes(
    max.peaks, penalties),
    data=max.evals[desired.peaks %in% c(10, 100, 1000)])+
  scale_y_log10(
    "GFPOP calls required
to compute model with
desired peaks (log scale)")+
  scale_x_log10(
    "Maximum number of peaks for data set")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap("desired.peaks", labeller=label_both)
png("jss-figure-more-evals.png", 5, 2, units="in", res=300)
print(gg)
dev.off()
