source("packages.R")

chunk.path <- "data/H3K4me3_TDH_immune/5"
peaks.path <- file.path(chunk.path, "peaks", "macs.trained.RData")
load(peaks.path)
counts.path <- file.path(chunk.path, "counts.RData")
load(counts.path)
counts.by.sample <- split(counts, counts$sample.id)

sid <- "McGill0104"
sample.counts <- data.table(counts.by.sample[[sid]])
sample.counts[, chromStart1 := chromStart+1]
setkey(sample.counts, chromStart1, chromEnd)
sample.peak.list <- lapply(peaks, subset, sample.id==sid)
n.peaks <- sapply(sample.peak.list, nrow)
n.peaks[order(as.numeric(names(n.peaks)))]

png.h <- 4
png.w <- 6
param.vec <- c(
  "1.30103"=-5,
  "7.5"=-10,
  "15"=-15)
param.labels <- data.table(
  param.name=names(param.vec),
  y=param.vec,
  PoissonLoss=NA_real_)
setkey(param.labels, param.name)
param.peaks.list <- list()
param.segs.list <- list()
y.val <- -5
x.val <- 118097
for(param.name in names(param.vec)){
  one.param <- data.table(sample.peak.list[[param.name]])
  change.vec <- one.param[, sort(c(chromStart, chromEnd))]
  segStart <- c(min(sample.counts$chromStart), change.vec)
  seg.dt <- data.table(
    segStart,
    segStart1=segStart+1,
    segEnd=c(change.vec, max(sample.counts$chromEnd)))
  setkey(seg.dt, segStart1, segEnd)
  over.dt <- foverlaps(sample.counts, seg.dt, nomatch=0L)
  over.dt[, bases := chromEnd-chromStart]
  stopifnot(
    over.dt[, sum(chromEnd-chromStart)]==
      sample.counts[, sum(chromEnd-chromStart)])
  seg.means <- over.dt[, list(
    mean=sum(coverage*bases)/sum(bases)
    ), by=.(segStart1, segEnd)]
  over.means <- foverlaps(sample.counts, seg.means, nomatch=0L)
  stopifnot(
    over.means[, sum(chromEnd-chromStart)]
    ==
      sample.counts[, sum(chromEnd-chromStart)]
    )
  ploss <- over.means[, PoissonLoss(coverage, mean, chromEnd-chromStart)]
  param.labels[param.name, PoissonLoss := ploss]
  param.peaks.list[[param.name]] <- data.table(
    param.name, one.param, y=param.vec[[param.name]])
  param.segs.list[[param.name]] <- data.table(
    param.name, seg.means)
  
  gg <- ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_segment(aes(chromStart/1e3, y.val,
                     xend=chromEnd/1e3, yend=y.val),
                 size=2,
                 data=one.param)+
    coord_cartesian(
      ylim=range(param.vec, sample.counts$coverage))+
    geom_step(aes(chromStart/1e3, coverage),
              color="grey50",
              data=sample.counts)+
    geom_vline(aes(xintercept=(segStart1-1)/1e3),
               data=seg.means[-1,],
               color="green",
               linetype="dashed")+
    geom_segment(aes((segStart1-1)/1e3, mean,
                     xend=segEnd/1e3, yend=mean),
                 size=1,
                 color="green",
                 data=seg.means)+
    geom_text(aes(
      x.val, y.val*2,
      label=sprintf(
        "macs log(qvalue)=%s, Poisson loss = %.1f",
        sub("0103", "", param.name), ploss)),
              hjust=0,
              ##size=3,
              data=param.labels[param.name,])+
    scale_y_continuous("aligned read counts")+
    scale_x_continuous("position on chromosome (kb = kilo bases)")
  print(gg)

  png.name <- sprintf("figure-macs-problem-%s.png", sub("[.]", "-", param.name))
  print(png.name)
  png(png.name,
      png.w, png.h, units="in", res=200)
  print(gg)
  dev.off()
  
}
param.peaks <- do.call(rbind, param.peaks.list)
param.segs <- do.call(rbind, param.segs.list)

param.breaks <- param.segs[segStart1!=min(segStart1),]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(param.name ~ .)+
  geom_vline(aes(xintercept=(segStart1-1)/1e3),
             data=param.breaks,
             color="green",
             linetype="dashed")+
  geom_segment(aes(chromStart/1e3, y.val,
                   xend=chromEnd/1e3, yend=y.val),
               size=2,
               data=param.peaks)+
  coord_cartesian(
    ylim=range(param.vec, sample.counts$coverage))+
  geom_text(aes(
    x.val, y.val,
    label=sprintf(
      "macs log(qvalue)=%s, Poisson loss = %.1f",
      sub("0103", "", param.name), PoissonLoss)),
            hjust=0,
            data=param.labels)+
  geom_step(aes(chromStart/1e3, coverage),
            color="grey50",
            data=sample.counts)+
  geom_segment(aes((segStart1-1)/1e3, mean,
                   xend=segEnd/1e3, yend=mean),
               size=1,
               color="green",
               data=param.segs)

sample.counts[, count := coverage]
fit <- PeakSegPDPAchrom(sample.counts, 2L)

peakseg.means <- subset(fit$segments, peaks==2)
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_segment(aes(chromStart/1e3, y.val,
                   xend=chromEnd/1e3, yend=y.val),
               size=2,
               data=subset(peakseg.means, status=="peak"))+
  geom_step(aes(chromStart/1e3, coverage),
            color="grey50",
            data=sample.counts)+
  geom_vline(aes(xintercept=chromStart/1e3),
             data=peakseg.means[-1,],
             color="green",
             linetype="dashed")+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               size=1,
               color="green",
               data=peakseg.means)+
  geom_text(aes(
    x.val, y.val*2,
    label=sprintf(
      "PeakSeg(maxPeaks=2), Poisson loss = %.1f",
      PoissonLoss)),
            hjust=0,
            ##size=3,
            data=fit$loss[3,])+
  coord_cartesian(
    ylim=range(param.vec, sample.counts$coverage))+
  scale_y_continuous("aligned read counts")+
  scale_x_continuous("position on chromosome (kb = kilo bases)")

png("figure-macs-problem-PeakSeg.png", png.w, png.h, units="in", res=200)
print(gg)
dev.off()

gg <- ggplot()+
  theme_bw()+
  geom_segment(aes(chromStart/1e3, y,
                   xend=chromEnd/1e3, yend=y),
               size=2,
               data=param.peaks)+
  geom_text(aes(x.val, y, label=paste0("macs log(qvalue)=",
                            sub("0103", "", param.name))),
            hjust=0,
            data=param.labels)+
  geom_step(aes(chromStart/1e3, coverage),
            color="grey50",
            data=sample.counts)+
  coord_cartesian(xlim=sample.counts[, range(chromStart, chromEnd)/1e3])+
  scale_y_continuous("aligned read counts")+
  scale_x_continuous("position on chromosome (kb = kilo bases)")
png("figure-macs-problem.png", png.w, png.h, units="in", res=200)
print(gg)
dev.off()

