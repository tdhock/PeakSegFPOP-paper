source("packages.R")

missing.RData.vec <- Sys.glob("data/H3K*/*/PDPA.missing/*")
loss.list <- list()
segs.list <- list()
breaks.list <- list()
counts.list <- list()
for(missing.i in seq_along(missing.RData.vec)){
  missing.RData <- missing.RData.vec[[missing.i]]
  (objs <- load(missing.RData))
  chunk.dir <- dirname(dirname(missing.RData))
  chunk.name <- sub("^data/", "", chunk.dir)
  counts.RData <- file.path(chunk.dir, "counts.RData")
  load(counts.RData)
  sid <- sub("[.]RData$", "", basename(missing.RData))
  one.sample <- data.table(counts)[sample.id==sid,]
  first.chromStart <- one.sample$chromStart[1]
  last.chromEnd <- one.sample[.N, chromEnd]
  bases <- last.chromEnd-first.chromStart
  normalize <- function(pos)(pos-first.chromStart)/bases
  one.sample[, normStart := normalize(chromStart)]
  one.sample[, normEnd := normalize(chromStart)]
  max.coverage <- max(one.sample$coverage)
  one.sample[, norm := coverage/max.coverage]
  meta <- data.table(chunk.name, sample.id=sid, missing.RData)
  counts.list[[missing.RData]] <- data.table(
    meta,
    one.sample)
  load(file.path(chunk.dir, "dp.model.reverse.RData"))
  dp.model.reverse <- dp.model
  load(file.path(chunk.dir, "dp.model.RData"))
  loss.list[[missing.RData]] <- rbind(with(dp.model.reverse[[sid]]$error, {
    data.table(
      meta,
      algorithm="cDPA.reverse",
      peaks, poisson.loss=error, constraint="inactive")
  }),
  with(dp.model[[sid]]$error, {
    data.table(
      meta,
      algorithm="cDPA.forward",
      peaks, poisson.loss=error, constraint="inactive")
  }),
  result$model$models[, data.table(
      meta,
    algorithm="PDPA", peaks, poisson.loss=min.cost, constraint)])
  segs.list[[missing.RData]] <- rbind(with(dp.model[[sid]]$segments, {
    data.table(
      meta,
      algorithm="cDPA.forward",
      peaks,
      chromStart, chromEnd,
      normStart=normalize(chromStart),
      normEnd=normalize(chromEnd),
      mean,
      norm=mean/max.coverage)
  }), with(dp.model.reverse[[sid]]$segments, {
    data.table(
      meta,
      algorithm="cDPA.reverse",
      peaks,
      chromStart, chromEnd,
      normStart=normalize(chromStart),
      normEnd=normalize(chromEnd),
      mean,
      norm=mean/max.coverage)
  }), result$model$segments[, data.table(
      meta,
      algorithm="PDPA",
      peaks,
      chromStart, chromEnd,
      normStart=normalize(chromStart),
      normEnd=normalize(chromEnd),
      mean=min.cost.mean,
      norm=min.cost.mean/max.coverage)
    ])
  breaks.list[[missing.RData]] <- rbind(with(dp.model[[sid]]$breaks, {
    data.table(
      meta,
      algorithm="cDPA.forward",
      peaks,
      normEnd=normalize(chromEnd),
      chromEnd)
  }), with(dp.model.reverse[[sid]]$breaks, {
    data.table(
      meta,
      algorithm="cDPA.reverse",
      peaks,
      normEnd=normalize(chromEnd),
      chromEnd)
  }), result$model$segments[segment.end<max(segment.end), data.table(
      meta,
      algorithm="PDPA",
      peaks,
      normEnd=normalize(chromEnd),
      chromEnd)
    ])
}
segs <- do.call(rbind, segs.list)
breaks <- do.call(rbind, breaks.list)
loss <- do.call(rbind, loss.list)
counts <- do.call(rbind, counts.list)
scatter.dt <- loss[, list(
  cDPA.models=sum(grepl("cDPA", algorithm)),
  PDPA.feasible=sum(constraint=="inactive" & algorithm=="PDPA")),
  by=.(missing.RData)]
  
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chunk.name + sample.id ~ algorithm, scales="free")+
  geom_point(aes(peaks, poisson.loss, color=constraint),
             data=loss)

abline.dt <- data.table(intercept=0, slope=1)
viz <- list(
  scatter=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(height=300)+
    geom_abline(aes(slope=slope, intercept=intercept),
                color="grey",
                data=abline.dt)+
    geom_point(aes(cDPA.models, PDPA.feasible, clickSelects=missing.RData),
               size=5,
               alpha=0.8,
               data=scatter.dt)+
    coord_equal(),
  loss=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(height=300)+
    facet_grid(. ~ algorithm)+
    geom_tallrect(aes(xmin=peaks-0.5, xmax=peaks+0.5,
                      clickSelects=peaks),
                  alpha=0.3,
                  data=data.table(peaks=0:9))+
    geom_point(aes(peaks, poisson.loss, color=constraint,
                   showSelected=missing.RData),
               data=loss),
  data=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=800)+
    facet_grid(algorithm ~ .)+
    geom_step(aes(normStart, norm,
                  showSelected=missing.RData),
              color="grey50",
              data=counts)+
    geom_text(aes(0.5, 0.8, label=sprintf("Poisson loss = %.1f", poisson.loss),
                  showSelected=missing.RData,
                  showSelected2=peaks),
              data=loss)+
    geom_segment(aes(normStart, norm,
                     showSelected=peaks,
                     showSelected2=missing.RData,
                     tooltip=paste("mean =", mean),
                     xend=normEnd, yend=norm),
                 color="green",
                 size=4,
                 data=segs)+
    geom_segment(aes(normEnd, 0,
                     xend=normEnd, yend=0.8,
                     showSelected=missing.RData,
                     showSelected2=peaks),
               color="green",
               linetype="dashed",
               data=breaks)
)
animint2dir(viz, "figure-cDPA-PDPA-all")
