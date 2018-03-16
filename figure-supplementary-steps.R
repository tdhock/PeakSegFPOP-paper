library(data.table)
fun.dt <- fread("supplementary.steps.csv")
unc.dt <- fread("supplementary.unconstrained.csv")
fit <- PeakSegDP::cDPA(as.integer(c(3, 9, 18, 15, 20, 2)), maxSegments=5L)
CDPA.dt <- with(fit, data.table(
  total.cost=as.numeric(loss),
  prev_seg_end=as.integer(ends),
  log.mean=log(as.numeric(mean)),
  timestep=as.integer(col(ends)),
  segments=as.integer(row(ends))))[segments <= timestep]
CDPA.dt[, cost := total.cost / timestep]

makeLines <- function(dt, l.grid=100){
  res <- dt[, {
    log.mean <- sort(unique(c(
      seq(min_log_mean, max_log_mean, l=l.grid),
      log(seq(exp(min_log_mean), exp(max_log_mean), l=l.grid)))))
    data.table(log.mean, cost=Log*log.mean+Linear*exp(log.mean)+Constant)
  }, by=list(segments=changes+1, timestep=data_point+1, fun, min_log_mean, max_log_mean, prev_seg_end=data_i+1)]
  res[, Fun := trans.vec[fun] ]
  res
}
lines.dt <- makeLines(fun.dt)
library(ggplot2)
trans.vec <- c(
  "cost model"="cost",
  "min prev cost"="minLess",
  "new cost model"="minEnv",
  "new cost model = min prev cost"="minLess",
  "new cost model + new data point"="new",
  "prev cost model"="prev")
points.dt <- lines.dt[log.mean==min_log_mean & 2 < exp(log.mean)]
text.dt <- lines.dt[, {
  L <- approx(log.mean, cost, log((exp(min_log_mean)+exp(max_log_mean))/2))
  with(L, data.table(
    log.mean=x,
    cost=y
  ))
}, by=list(segments, timestep, min_log_mean, max_log_mean, Fun, prev_seg_end)]

i.colors <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
  "#A6761D", "#666666")

ggplot(lines.dt, aes(exp(log.mean), cost))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ timestep, labeller=label_both)+
  geom_line(aes(
    color=factor(prev_seg_end), group=paste(
      fun, timestep)))+
  geom_point(data=points.dt, shape=1)+
  directlabels::geom_dl(aes(
    label=Fun),
    method="last.qp")+
  xlim(2, 25)

ggplot(lines.dt[Fun=="new"], aes(exp(log.mean), cost))+
  geom_point(data=points.dt[Fun=="new"], size=1)+
  ##geom_text(aes(label=prev_seg_end, color=factor(prev_seg_end)), vjust=-0.5, data=text.dt[Fun=="new"])+
  ##geom_point(data=CDPA.dt)+
  theme_bw()+
  theme(
    panel.margin=grid::unit(0, "lines"),
    legend.position="bottom"
  )+
  facet_grid(segments ~ timestep, labeller=label_both)+
  geom_line(aes(
    color=factor(prev_seg_end), group=paste(
      min_log_mean)),
    size=1)+
  scale_color_manual(values=i.colors, guide=guide_legend(ncol=10))

unc.lines.dt <- makeLines(unc.dt)
ggplot(unc.lines.dt[Fun=="new"], aes(exp(log.mean), cost))+
  theme_bw()+
  theme(
    panel.margin=grid::unit(0, "lines"),
    legend.position="bottom"
  )+
  facet_grid(segments ~ timestep, labeller=label_both)+
  geom_line(aes(
    color=factor(prev_seg_end), group=paste(
      min_log_mean)),
    size=1)+
  scale_color_manual(values=i.colors, guide=guide_legend(ncol=10))


ggplot(, aes(exp(log.mean), cost))+
  theme_bw()+
  theme(
    panel.margin=grid::unit(0, "lines"),
    legend.position="bottom"
  )+
  facet_grid(segments ~ timestep, labeller=label_both)+
  geom_line(aes(
    color=model, group=paste(
      model, min_log_mean)),
    data=rbind(
      data.table(model="constrained", lines.dt[Fun=="new"]),
      data.table(model="unconstrained", unc.lines.dt[Fun=="new"])),
    size=0.5)
  scale_color_manual(values=i.colors, guide=guide_legend(ncol=10))

both.3.4 <- rbind(
  data.table(model="constrained", lines.dt),
  data.table(model="unconstrained", unc.lines.dt)
)[segments==3 & timestep==4]
both.3.4[Fun=="minLess" & model=="unconstrained", Fun := "uncMin"]
ggplot(both.3.4, aes(exp(log.mean), cost))+
  coord_cartesian( xlim=c(3, 25), ylim=c(-14.6, -13.6))+
  theme_bw()+
  theme(
    panel.margin=grid::unit(0, "lines"),
    legend.position="bottom"
  )+
  facet_grid(. ~ model, labeller=label_both)+
  geom_line(aes(
    color=factor(prev_seg_end), group=paste(
      model, Fun, min_log_mean)),
    size=0.5)+
  scale_color_manual(values=i.colors, guide=guide_legend(ncol=10))+
  directlabels::geom_dl(aes(
    label=Fun),
    method="last.qp")


ggplot(lines.dt[segments %in% c(2, 3) & timestep %in% c(3,4)], aes(exp(log.mean), cost))+
  coord_cartesian( xlim=c(3, 25), ylim=c(-14.6, -13.6))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ timestep, labeller=label_both)+
  geom_line(aes(
    color=factor(prev_seg_end), group=paste(
      fun, timestep)),
    size=1)+
  directlabels::geom_dl(aes(
    label=Fun),
    method="last.qp")
