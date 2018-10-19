source("packages.R")

## TODO: check problematic sample H3K4me3_XJ_immune/1 McGill0024 where
## PDPA loss is bigger than cDPA loss.

PDPA.RData.vec <- Sys.glob("../chip-seq-paper/chunks/H3K*/*/PDPA.model/*")
"H3K4me3_PGP_immune/10/"
"H3K4me3_PGP_immune/2/"
"H3K4me3_XJ_immune/1/"
##PDPA.RData.vec <- grep("H3K4me3_XJ_immune/1/|H3K4me3_PGP_immune/2/|H3K4me3_PGP_immune/10/", PDPA.RData.vec, value=TRUE)

algo <- function(x){
  factor(x, c("PDPA intervals", "PDPA", "cDPA.forward", "cDPA.reverse"))
}
loss.list <- list()
segs.list <- list()
breaks.list <- list()
counts.list <- list()
intervals.list <- list()
for(PDPA.i in seq_along(PDPA.RData.vec)){
  PDPA.RData <- PDPA.RData.vec[[PDPA.i]]
  cat(sprintf("%4d / %4d %s\n", PDPA.i, length(PDPA.RData.vec), PDPA.RData))
  (objs <- load(PDPA.RData))
  chunk.dir <- dirname(dirname(PDPA.RData))
  chunk.name <- sub("^chunks/", "", chunk.dir)
  counts.RData <- file.path(chunk.dir, "counts.RData")
  load(counts.RData)
  sid <- sub("[.]RData$", "", basename(PDPA.RData))
  one.sample <- data.table(counts)[sample.id==sid,]
  ##fit <- cDPA(one.sample$coverage, one.sample[,chromEnd-chromStart], 19)
  first.chromStart <- one.sample$chromStart[1]
  last.chromEnd <- one.sample[.N, chromEnd]
  bases <- last.chromEnd-first.chromStart
  normalize <- function(pos)(pos-first.chromStart)/bases
  sample.intervals <- data.table(
    algorithm=algo("PDPA intervals"), PDPA.RData, result$model$intervals)
  sample.intervals[, chromStart := one.sample$chromStart[timestep]]
  sample.intervals[, normStart := normalize(chromStart)]
  sample.intervals[, normIntervals := intervals/max(intervals)]
  sample.intervals[, peaks := (total.segments-1)/2]
  sample.intervals$data <- nrow(one.sample)
  intervals.list[[PDPA.RData]] <- sample.intervals
  one.sample[, normStart := normalize(chromStart)]
  one.sample[, normEnd := normalize(chromStart)]
  max.coverage <- max(one.sample$coverage)
  one.sample[, norm := coverage/max.coverage]
  meta <- data.table(chunk.name, sample.id=sid, PDPA.RData)
  norm.pos <- seq(0, 1, l=200)
  counts.list[[PDPA.RData]] <- one.sample[, data.table(
    meta,
    norm=approx(normStart, norm, norm.pos)$y,
    normStart=norm.pos)]
  counts.list[[PDPA.RData]] <- data.table(meta, one.sample)
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
  loss.list[[PDPA.RData]] <- sample.loss
  segs.list[[PDPA.RData]] <- rbind(
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
  breaks.list[[PDPA.RData]] <- rbind(
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
  intervals[, list(mean.intervals=mean(intervals)), by=PDPA.RData]
loss.wide <- dcast(
  loss,
  PDPA.RData + chunk.name + sample.id + peaks ~ algorithm,
  value.var="poisson.loss")
loss.wide[, poisson.loss.diff.fwd := cDPA.forward-PDPA]
loss.wide[, poisson.loss.diff.rev := cDPA.reverse-PDPA]
loss.wide[order(poisson.loss.diff.fwd),] # TODO: check these!
loss.diff.min.max <- loss.wide[, list(
  min=min(poisson.loss.diff.fwd, poisson.loss.diff.rev),
  max=max(poisson.loss.diff.fwd, poisson.loss.diff.rev)),
  by=.(PDPA.RData, chunk.name, sample.id)]
scatter.some <- loss[, list(
  cDPA.fwd.models=sum(grepl("cDPA.forward", algorithm)),
  cDPA.rev.models=sum(grepl("cDPA.rev", algorithm)),
  PDPA.feasible=sum(constraint=="inactive" & algorithm=="PDPA")),
  by=.(PDPA.RData)]
scatter.some[, cDPA.mean.models := (cDPA.fwd.models + cDPA.rev.models)/2]
setkey(interval.means, PDPA.RData)
setkey(scatter.some, PDPA.RData)
setkey(loss.diff.min.max, PDPA.RData)
scatter.dt <- scatter.some[interval.means][loss.diff.min.max]
interval.max <-
  intervals[, .SD[which.max(intervals),], by=.(PDPA.RData, peaks)]
interval.stats <-
  intervals[, list(
    max=max(intervals),
    q75=quantile(intervals, 0.75),
    median=as.numeric(median(intervals)),
    q25=quantile(intervals, 0.25),
    min=min(intervals)
  ), by=.(PDPA.RData, peaks)]
  
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chunk.name + sample.id ~ algorithm, scales="free")+
  geom_point(aes(peaks, poisson.loss, color=constraint),
             data=loss)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("PDPA.RData")+
  coord_equal()+
  geom_point(aes(poisson.loss.diff.rev, poisson.loss.diff.fwd),
             shape=1,
             data=loss.wide)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("PDPA.RData")+
  geom_line(aes(peaks, poisson.loss.diff.fwd),
            data=loss.wide)+
  geom_line(aes(peaks, poisson.loss.diff.rev),
            data=loss.wide)

abline.dt <- data.table(intercept=0, slope=1)
loss[, feasible := ifelse(constraint=="active", "infeasible","feasible")]
loss[, optimal := ifelse(
  poisson.loss-.SD[algorithm=="PDPA", poisson.loss] < 1e-2,
  "optimal",
  "sub-optimal"), by=.(PDPA.RData, peaks)]
viz <- list(
  title="Algorithms for computing PeakSeg model",
  feasibilty=ggplot()+
    ggtitle("select data set")+
    theme_grey()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(height=300, width=200)+
    xlab("mean PDPA intervals")+
    ylab("PDPA feasible - cDPA models")+
    scale_fill_gradient("max loss diff", low="black", high="white")+
    scale_color_gradient("min loss diff", low="white", high="blue")+
    geom_point(aes(mean.intervals, PDPA.feasible-cDPA.mean.models,
                   color=min, fill=max,
                   clickSelects=PDPA.RData),
               alpha=0.8,
               shape=21,
               size=4,
               data=scatter.dt),
  optimality=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(height=300, width=200)+
    geom_line(aes(peaks, poisson.loss.diff.fwd,
                  showSelected=PDPA.RData),
              data=loss.wide)+
    geom_line(aes(peaks, poisson.loss.diff.rev,
                  showSelected=PDPA.RData),
              data=loss.wide)+
    geom_point(aes(peaks, poisson.loss.diff.rev,
                   showSelected=PDPA.RData),
               data=loss.wide[poisson.loss.diff.rev < -1,])+
    geom_tallrect(aes(xmin=peaks-0.5, xmax=peaks+0.5,
                      clickSelects=peaks),
                  alpha=0.2,
                  data=data.table(peaks=0:9)),
  intervals=ggplot()+
    ggtitle("quartiles of intervals")+
    theme_animint(height=300, width=250)+
    scale_x_continuous(breaks=0:9)+
    ylab("intervals")+
    geom_segment(aes(peaks, min,
                     xend=peaks, yend=q25,
                     showSelected=PDPA.RData),
                 data=interval.stats)+
    geom_segment(aes(peaks, max,
                     xend=peaks, yend=q75,
                     showSelected=PDPA.RData),
                 data=interval.stats)+
    geom_point(aes(peaks, median,
                   showSelected=PDPA.RData),
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
                   showSelected=PDPA.RData),
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
                  showSelected=PDPA.RData),
              size=3,
              data=intervals)+
    geom_text(aes(normStart, normIntervals, label=paste(
      intervals, "intervals /", data, "data points"),
                  showSelected2=peaks,
                  showSelected=PDPA.RData),
              data=interval.max)+
    geom_step(aes(normStart, norm,
                  showSelected=PDPA.RData),
              color="grey50",
              data=counts)+
    geom_text(aes(0.5, 0.8, label=sprintf("Poisson loss = %.1f", poisson.loss),
                  showSelected=PDPA.RData,
                  showSelected2=peaks),
              color="green",
              data=loss)+
    geom_segment(aes(normStart, norm,
                     showSelected=peaks,
                     showSelected2=PDPA.RData,
                     tooltip=paste("mean =", mean),
                     xend=normEnd, yend=norm),
                 color="green",
                 alpha=0.5,
                 size=4,
                 data=segs)+
    geom_segment(aes(normEnd, 0,
                     xend=normEnd, yend=0.8,
                     showSelected=PDPA.RData,
                     showSelected2=peaks),
                 color="green",
                 alpha=0.5,
                 linetype="dashed",
                 data=breaks)+
    scale_x_continuous("relative position on chromosome", breaks=c())+
    scale_y_continuous("relative aligned read counts / intervals", breaks=c())
)
viz$feasibilty

animint2dir(viz, "figure-cDPA-PDPA-all")
