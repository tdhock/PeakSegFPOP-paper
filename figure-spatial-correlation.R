library(data.table)
library(ggplot2)

files.dt <- rbind(
  data.table(sample.id="McGill0004", experiment="H3K36me3", chrom="chr9", chunk="H3K36me3_AM_immune/8"),
  data.table(sample.id="McGill0002", experiment="H3K4me3", chrom="chr2", chunk="H3K4me3_PGP_immune/7"))

peaks.dt.list <- list()
cov.dt.list <- list()
for(file.i in 1:nrow(files.dt)){
  f <- files.dt[file.i]
  bed.gz <- f[, paste0(sample.id, "_", experiment, ".bed.gz")]
  if(!file.exists(bed.gz)){
    suffix <- f[, paste0(sample.id, "/", experiment)]
    u <- paste0("http://hubs.hpc.mcgill.ca/~thocking/bed/", suffix)
    download.file(u, bed.gz)
  }
  reads.dt <- fread(paste0("zcat ", bed.gz, "|grep ^", f$chrom), select=2:3)
  setnames(reads.dt, c("chromStart", "chromEnd"))
  load(paste0("data/", f$chunk, "/counts.RData"))
  sample.dt <- data.table(counts)[sample.id==f$sample.id]
  first <- sample.dt$chromStart[1]
  last <- sample.dt[, chromEnd[.N] ]
  some.reads <- reads.dt[!(chromEnd < first | last < chromStart)]
  end.counts <- some.reads[, list(count=.N), by=list(chromEnd)]
  end.counts[, chromStart := chromEnd-1L]

  ggplot()+
    theme_bw()+
    geom_rect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3,
      ymin=0, ymax=coverage),
              data=sample.dt)+
    geom_point(aes(
      chromEnd/1e3, count),
               shape=1,
               data=end.counts)

  ## Create rle/compressed data for PeakSegPDPA.
  u.pos <- end.counts[, sort(unique(c(chromStart, chromEnd, first, last)))]
  zero.cov <- data.table(
    chromStart = u.pos[-length(u.pos)], 
    chromEnd = u.pos[-1], count = 0L)
  setkey(zero.cov, chromEnd)
  zero.cov[J(end.counts$chromEnd), `:=`(count, end.counts$count)]
  dup.cov <- zero.cov[first <= chromStart & chromEnd <= last]
  out.cov <- dup.cov[c(diff(count), Inf)!=0]
  out.cov[, chromStart := c(first, chromEnd[-.N])]
  out.cov[, stopifnot(chromEnd[-.N] == chromStart[-1])]
  out.cov[, stopifnot(all(diff(count)!=0))]

  count.list <- list(
    coverage=sample.dt[, list(chromStart, chromEnd, count=coverage)],
    last=out.cov)

  for(count.method in names(count.list)){
    count.dt <- count.list[[count.method]]
    count.dt[, stopifnot(
      chromEnd[-.N] == chromStart[-1],
      all(diff(count)!=0),
      first==chromStart[1],
      last==chromEnd[.N])]
    fit <- PeakSegOptimal::PeakSegPDPAchrom(count.dt, 2L)
    cov.dt.list[[paste(count.method, file.i)]] <- data.table(count.method, f, count.dt)
    peaks.dt.list[[paste(count.method, file.i)]] <- data.table(count.method, f, subset(fit$segments, peaks==2))
  }
}

cov.dt <- do.call(rbind, cov.dt.list)
peaks.dt <- do.call(rbind, peaks.dt.list)[status=="peak"]

cov.dt[, list(max=max(count)), by=list(experiment)]
scale.dt <- data.table(
  y=15, count.method="last",
  yy=-c(54, 377)/10,
  start=c(111750, 175440),
  end=c(111800, 175490),
  experiment=c("H3K36me3", "H3K4me3"))
just.diff <- 0.3
max.dt <- cov.dt[, list(
  max.count=max(count)
  ), by=count.method]
blank.dt <- cov.dt[, list(min.chromStart=min(chromStart)), by=list(count.method, experiment)][max.dt, on=list(count.method)]
size <- 3
cov.peaks <- do.call(rbind, peaks.dt.list)[count.method=="coverage"]
gg <- ggplot()+
  theme_bw()+
  geom_segment(aes(
    start, yy,
    xend=end, yend=yy),
               color="grey",
               size=4,
               data=scale.dt)+
  geom_text(aes(
    (start+end)/2, yy, label=paste(end-start, "kb")),
               color="grey",
             vjust=-0.7,
               data=scale.dt)+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("experiment", scales="free", labeller=function(df){
    if("count.method" %in% names(df)){
      df$count.method <- c(
        ## coverage="each read counted at 100 positions, one for each aligned base: spatial correlation present",
        ## last="each read counted at one position, the last aligned base: spatial correlation absent"
        coverage="Some spatial correlation",
        last="No spatial correlation"
        )[df$count.method]
    }
    if("experiment" %in% names(df)){
      df$experiment <- c(
        H3K36me3="Broad histone mark H3K36me3",
        H3K4me3="Sharp histone mark H3K4me3"
        )[df$experiment]
    }
    df
  })+
  geom_step(aes(
    chromStart/1e3, count),
            data=cov.dt[count.method=="coverage"])+
  geom_segment(aes(
    chromStart/1e3, mean,
    xend=chromEnd/1e3, yend=mean),
               color="deepskyblue",
               size=1,
               data=cov.peaks)+
  xlab("position on chromosome (kb = kilo bases)")+
  scale_y_continuous("aligned read coverage")+
  coord_cartesian(expand=FALSE)
png("figure-spatial-correlation-mean.png", 1800, 600, res=200)
print(gg)
dev.off()


gg <- ggplot()+
  theme_bw()+
  geom_segment(aes(
    start, y,
    xend=end, yend=y),
               color="grey",
               size=2,
               data=scale.dt)+
  geom_text(aes(
    (start+end)/2, y, label=paste(end-start, "kb")),
               color="grey",
             vjust=-0.5, 
               data=scale.dt)+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(count.method ~ experiment, scales="free", labeller=function(df){
    if("count.method" %in% names(df)){
      df$count.method <- c(
        ## coverage="each read counted at 100 positions, one for each aligned base: spatial correlation present",
        ## last="each read counted at one position, the last aligned base: spatial correlation absent"
        coverage="Some spatial correlation",
        last="No spatial correlation"
        )[df$count.method]
    }
    if("experiment" %in% names(df)){
      df$experiment <- c(
        H3K36me3="Broad histone mark H3K36me3",
        H3K4me3="Sharp histone mark H3K4me3"
        )[df$experiment]
    }
    df
  })+
  geom_step(aes(
    chromStart/1e3, count),
            data=cov.dt)+
  geom_segment(aes(
    chromStart/1e3, mean,
    xend=chromEnd/1e3, yend=mean),
               color="deepskyblue",
               size=1,
               data=peaks.dt)+
  geom_text(aes(
    chromStart/1e3, 0, label=as.integer(chromStart/1e3)),
            color="deepskyblue",
            vjust=1.5,
            size=size,
            hjust=1-just.diff,
               data=peaks.dt)+
  geom_text(aes(
    chromEnd/1e3, 0, label=as.integer(chromEnd/1e3)),
            color="deepskyblue",
            vjust=2.7,
            size=size,
            hjust=just.diff,
            data=peaks.dt)+
  geom_blank(aes(
    min.chromStart/1e3, -max.count/10),
             data=blank.dt)+
  xlab("position on chromosome (kb = kilo bases)")+
  scale_y_continuous("aligned read counts")
print(gg)

png("figure-spatial-correlation.png", 1800, 1200, res=200)
print(gg)
dev.off()
