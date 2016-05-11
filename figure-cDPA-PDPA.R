source("packages.R")

load("data/H3K4me3_XJ_immune/2/dp.model.reverse.RData")
load("data/H3K4me3_XJ_immune/2/counts.RData")
dp.model.reverse <- dp.model
load("data/H3K4me3_XJ_immune/2/dp.model.RData")
sid <- "McGill0102"
load("H3K4me3_XJ_immune_chunk2_McGill0102_PDPA_result.RData")
one.sample <- data.table(counts)[sample.id==sid,]

loss.dt <- rbind(with(dp.model.reverse[[sid]]$error, {
  data.table(
    algorithm="cDPA.reverse",
    peaks, poisson.loss=error, constraint="inactive")
}),
with(dp.model[[sid]]$error, {
  data.table(
    algorithm="cDPA.forward",
    peaks, poisson.loss=error, constraint="inactive")
}),
result$model$models[, data.table(
  algorithm="PDPA", peaks, poisson.loss=min.cost, constraint)])

segment.dt <- rbind(with(dp.model[[sid]]$segments, {
  data.table(
    algorithm="cDPA.forward",
    peaks,
    chromStart, chromEnd, mean)
}), with(dp.model.reverse[[sid]]$segments, {
  data.table(
    algorithm="cDPA.reverse",
    peaks,
    chromStart, chromEnd, mean)
}), result$model$segments[, data.table(
  algorithm="PDPA",
  peaks,
  chromStart, chromEnd,
  mean=min.cost.mean)
])

break.dt <- rbind(with(dp.model[[sid]]$breaks, {
  data.table(
    algorithm="cDPA.forward",
    peaks,
    chromEnd)
}), with(dp.model.reverse[[sid]]$breaks, {
  data.table(
    algorithm="cDPA.reverse",
    peaks,
    chromEnd)
}), result$model$segments[segment.end<max(segment.end), data.table(
  algorithm="PDPA",
  peaks,
  chromEnd)
])

ggplot()+
  geom_point(aes(peaks, poisson.loss),
             data=loss.dt[constraint=="active",])+
  geom_line(aes(peaks, poisson.loss, group=algorithm),
            data=loss.dt)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(algorithm ~ peaks)+
  geom_step(aes(chromStart/1e3, coverage),
            color="grey50",
            data=one.sample)+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               color="green",
               data=segment.dt)+
  geom_vline(aes(xintercept=(chromEnd+0.5)/1e3), 
             color="green",
             linetype="dashed",
             data=break.dt)

viz <- list(
  loss=ggplot()+
    geom_tallrect(aes(xmin=peaks-0.5, xmax=peaks+0.5,
                      clickSelects=peaks),
                  alpha=0.5,
                  data=data.table(peaks=0:9))+
    geom_line(aes(peaks, poisson.loss, group=algorithm,
                  clickSelects=algorithm),
              size=4,
              data=loss.dt)+
    geom_point(aes(peaks, poisson.loss, color=constraint),
               data=loss.dt[constraint=="active",]),
  data=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=800)+
    facet_grid(algorithm ~ .)+
    geom_step(aes(chromStart/1e3, coverage),
              color="grey50",
              data=one.sample)+
    geom_segment(aes(chromStart/1e3, mean,
                     showSelected=peaks,
                     clickSelects=algorithm,
                     tooltip=paste("mean =", mean),
                     xend=chromEnd/1e3, yend=mean),
                 color="green",
                 size=4,
                 alpha=0.8,
                 data=segment.dt)+
    geom_vline(aes(xintercept=(chromEnd+0.5)/1e3,
                   clickSelects=algorithm,
                   showSelected=peaks),
               alpha=0.8,
               color="green",
               linetype="dashed",
               data=break.dt)
)
animint2dir(viz, "figure-cDPA-PDPA")
