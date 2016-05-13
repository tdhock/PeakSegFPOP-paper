source("packages.R")

missing.RData.vec <- Sys.glob("data/H3K*/*/PDPA.missing/*")
algo <- function(x){
  factor(x, c("PDPA intervals", "PDPA", "cDPA.forward", "cDPA.reverse"))
}
loss.list <- list()
segs.list <- list()
breaks.list <- list()
counts.list <- list()
intervals.list <- list()
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
  sample.intervals <- data.table(
    algorithm=algo("PDPA intervals"), missing.RData, result$model$intervals)
  sample.intervals[, chromStart := one.sample$chromStart[timestep]]
  sample.intervals[, normStart := normalize(chromStart)]
  sample.intervals[, normIntervals := intervals/max(intervals)]
  sample.intervals[, peaks := (total.segments-1)/2]
  sample.intervals$data <- nrow(one.sample)
  intervals.list[[missing.RData]] <- sample.intervals
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
  sample.loss <- rbind(
    result$model$models[, data.table(
      meta,
      algorithm=algo("PDPA"), peaks, poisson.loss=min.cost, constraint)],    
    with(dp.model.reverse[[sid]]$error, {
      data.table(
        meta,
        algorithm=algo("cDPA.reverse"),
        peaks, poisson.loss=error, constraint="inactive")
    }),
    with(dp.model[[sid]]$error, {
      data.table(
        meta,
        algorithm=algo("cDPA.forward"),
        peaks, poisson.loss=error, constraint="inactive")
    }))
  min.loss <- min(sample.loss$poisson.loss)
  sample.loss[, normLoss := (poisson.loss-min.loss)/(max(poisson.loss)-min.loss)]
  loss.list[[missing.RData]] <- sample.loss
  segs.list[[missing.RData]] <- rbind(
    result$model$segments[, data.table(
      meta,
      algorithm=algo("PDPA"),
      peaks,
      chromStart, chromEnd,
      normStart=normalize(chromStart),
      normEnd=normalize(chromEnd),
      mean=min.cost.mean,
      norm=min.cost.mean/max.coverage)
      ],
    with(dp.model[[sid]]$segments, {
    data.table(
      meta,
      algorithm=algo("cDPA.forward"),
      peaks,
      chromStart, chromEnd,
      normStart=normalize(chromStart),
      normEnd=normalize(chromEnd),
      mean,
      norm=mean/max.coverage)
  }), with(dp.model.reverse[[sid]]$segments, {
    data.table(
      meta,
      algorithm=algo("cDPA.reverse"),
      peaks,
      chromStart, chromEnd,
      normStart=normalize(chromStart),
      normEnd=normalize(chromEnd),
      mean,
      norm=mean/max.coverage)
  }))
  breaks.fwd <- dp.model[[sid]]$breaks
  breaks.rev <- dp.model.reverse[[sid]]$breaks
  breaks.list[[missing.RData]] <- rbind(
    result$model$segments[segment.end<max(segment.end), data.table(
      meta,
      algorithm=algo("PDPA"),
      peaks,
      normEnd=normalize(chromEnd),
      chromEnd)
    ],
    if(is.data.frame(breaks.fwd))with(breaks.fwd, {
    data.table(
      meta,
      algorithm=algo("cDPA.forward"),
      peaks,
      normEnd=normalize(chromEnd),
      chromEnd)
  }), if(is.data.frame(breaks.rev))with(breaks.rev, {
    data.table(
      meta,
      algorithm=algo("cDPA.reverse"),
      peaks,
      normEnd=normalize(chromEnd),
      chromEnd)
  }))
}
segs <- do.call(rbind, segs.list)
breaks <- do.call(rbind, breaks.list)
loss <- do.call(rbind, loss.list)
counts <- do.call(rbind, counts.list)
intervals <- do.call(rbind, intervals.list)
interval.means <-
  intervals[, list(mean.intervals=mean(intervals)), by=missing.RData]
scatter.some <- loss[, list(
  cDPA.models=sum(grepl("cDPA", algorithm)),
  PDPA.feasible=sum(constraint=="inactive" & algorithm=="PDPA")),
  by=.(missing.RData)]
setkey(interval.means, missing.RData)
setkey(scatter.some, missing.RData)
scatter.dt <- scatter.some[interval.means]
interval.max <-
  intervals[, .SD[which.max(intervals),], by=.(missing.RData, peaks)]
interval.stats <-
  intervals[, list(
    max=max(intervals),
    q75=quantile(intervals, 0.75),
    median=as.numeric(median(intervals)),
    q25=quantile(intervals, 0.25),
    min=min(intervals)
  ), by=.(missing.RData, peaks)]
  
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chunk.name + sample.id ~ algorithm, scales="free")+
  geom_point(aes(peaks, poisson.loss, color=constraint),
             data=loss)

abline.dt <- data.table(intercept=0, slope=1)
loss[, feasible := ifelse(constraint=="active", "infeasible","feasible")]
loss[, optimal := ifelse(
  poisson.loss-.SD[algorithm=="PDPA", poisson.loss] < 1e-2,
  "optimal",
  "sub-optimal"), by=.(missing.RData, peaks)]
viz <- list(
  title="Algorithms for computing PeakSeg model",
  scatter=ggplot()+
    ggtitle("select data set")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(height=300, width=200)+
    geom_abline(aes(slope=slope, intercept=intercept),
                color="grey",
                data=abline.dt)+
    xlab("all cDPA models")+
    ylab("PDPA feasible models")+
    geom_rect(aes(xmin=cDPA.models-0.5, xmax=cDPA.models+0.5,
                  ymin=PDPA.feasible-0.5, ymax=PDPA.feasible+0.5,
                  fill=mean.intervals,
                  clickSelects=missing.RData),
               data=scatter.dt),
  intervals=ggplot()+
    ggtitle("quartiles of intervals")+
    theme_animint(height=300, width=250)+
    scale_x_continuous(breaks=0:9)+
    ylab("intervals")+
    geom_segment(aes(peaks, min,
                     xend=peaks, yend=q25,
                     showSelected=missing.RData),
                 data=interval.stats)+
    geom_segment(aes(peaks, max,
                     xend=peaks, yend=q75,
                     showSelected=missing.RData),
                 data=interval.stats)+
    geom_point(aes(peaks, median,
                   showSelected=missing.RData),
               data=interval.stats)+
    geom_tallrect(aes(xmin=peaks-0.5, xmax=peaks+0.5,
                      clickSelects=peaks),
                  alpha=0.2,
                  data=data.table(peaks=0:9)),
  loss=ggplot()+
    ggtitle("PeakSeg objective and feasibility")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(height=300, width=350)+
    facet_grid(. ~ algorithm)+
    scale_x_continuous(breaks=seq(1, 9,by=2))+
    scale_color_manual(values=c("sub-optimal"="black", optimal="white"))+
    geom_point(aes(peaks, normLoss,
                   color=optimal,
                   fill=feasible,
                   showSelected=missing.RData),
               size=4,
               data=loss)+
    geom_tallrect(aes(xmin=peaks-0.5, xmax=peaks+0.5,
                      clickSelects=peaks),
                  alpha=0.2,
                  data=data.table(peaks=0:9))+
    scale_y_continuous("relative Poisson loss", breaks=c()),
  data=ggplot()+
    ggtitle("segmentation models and intervals")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=1000)+
    facet_grid(algorithm ~ .)+
    geom_line(aes(normStart, normIntervals, group=peaks,
                  showSelected2=peaks,
                  showSelected=missing.RData),
              size=3,
              data=intervals)+
    geom_text(aes(normStart, normIntervals, label=paste(
      intervals, "intervals /", data, "data points"),
                  showSelected2=peaks,
                  showSelected=missing.RData),
              data=interval.max)+
    geom_step(aes(normStart, norm,
                  showSelected=missing.RData),
              color="grey50",
              data=counts)+
    geom_text(aes(0.5, 0.8, label=sprintf("Poisson loss = %.1f", poisson.loss),
                  showSelected=missing.RData,
                  showSelected2=peaks),
              color="green",
              data=loss)+
    geom_segment(aes(normStart, norm,
                     showSelected=peaks,
                     showSelected2=missing.RData,
                     tooltip=paste("mean =", mean),
                     xend=normEnd, yend=norm),
                 color="green",
                 alpha=0.5,
                 size=4,
                 data=segs)+
    geom_segment(aes(normEnd, 0,
                     xend=normEnd, yend=0.8,
                     showSelected=missing.RData,
                     showSelected2=peaks),
                 color="green",
                 alpha=0.5,
                 linetype="dashed",
                 data=breaks)+
    scale_x_continuous("relative position on chromosome", breaks=c())+
    scale_y_continuous("relative aligned read counts / intervals", breaks=c())
)
animint2dir(viz, "figure-cDPA-PDPA-all")
