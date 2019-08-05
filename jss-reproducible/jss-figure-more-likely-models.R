source("jss-packages.R")

chunk.dir <- "jss-figure-more-likely-models"
counts.RData <- file.path(chunk.dir, "counts.RData")
load(counts.RData)
counts.dt <- data.table(counts)

sample.num.vec <- c(#101,
  ##322,
  ##4, 91,
  104
)
sample.id.vec <- sprintf("McGill%04d", sample.num.vec)
some.counts <- counts.dt[sample.id %in% sample.id.vec]

peaks.RData <- file.path(chunk.dir, "macs.default.RData")
load(peaks.RData)
peaks.dt <- data.table(peaks[[1]])
some.peaks <- peaks.dt[sample.id %in% sample.id.vec]

PoissonLogLik <- function(data.dt, seg.dt.in){
  dataStart <- data.dt$chromStart[1]
  dataEnd <- data.dt[.N, chromEnd]
  seg.dt <- data.table(seg.dt.in)
  seg.dt[, segStart1 := segStart+1]
  data.dt[, chromStart1 := chromStart+1]
  setkey(seg.dt, segStart1, segEnd)
  setkey(data.dt, chromStart1, chromEnd)
  over.dt <- foverlaps(data.dt, seg.dt, nomatch=0L)
  over.dt[chromStart < segStart, chromStart := segStart]
  over.dt[segEnd < chromEnd, chromEnd := segEnd]
  over.dt[, weight := chromEnd-chromStart]
  stopifnot(dataEnd-dataStart == sum(over.dt$weight))
  seg.means <- over.dt[, list(
    mean=sum(coverage*weight)/sum(weight)
  ), by=list(segStart, segEnd)]
  seg.means[, status := rep(c("background", "peak"), l=.N)]
  data.means <- seg.means[over.dt, on=list(segStart, segEnd)]
  data.means[, list(
    total.logLik=sum(dpois(coverage, mean, log=TRUE)*weight),
    total.weight=sum(weight)
    )]
}

segs.list <- list()
loss.list <- list()
logLik.list <- list()
chrom <- "chr11"
n.peaks.vec <- 5:1
n.peaks.vec <- as.integer(c(5, 3))
for(sample.i in seq_along(sample.id.vec)){
  sample.id <- sample.id.vec[[sample.i]]
  cat(sprintf("%4d / %4d %s\n", sample.i, length(sample.id.vec), sample.id))
  s.dt <- data.table(sample.id)
  sample.counts <- counts.dt[s.dt, on=list(sample.id)]
  sample.peaks <- peaks.dt[s.dt, on=list(sample.id)]
  problemStart <- sample.counts$chromStart[1]
  problemEnd <- sample.counts[.N, chromEnd]
  problem <- data.table(chrom, problemStart, problemEnd)
  problem.dir <- file.path(
    "jss-data", sample.id,
    sprintf("%s:%d-%d", chrom, problemStart, problemEnd))
  dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
  problem.bed <- file.path(problem.dir, "problem.bed")
  fwrite(problem, problem.bed, sep="\t", col.names=FALSE)
  coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  coverage.dt <- sample.counts[, data.table(
    chrom, chromStart, chromEnd, coverage)]
  fwrite(coverage.dt, coverage.bedGraph, sep="\t", col.names=FALSE)
  seg.limit.vec <- c(
    problemStart,
    sample.peaks[, rbind(chromStart, chromEnd)],
    problemEnd)
  n.limits <- length(seg.limit.vec)
  seg.dt <- data.table(
    segStart=seg.limit.vec[-n.limits],
    segEnd=seg.limit.vec[-1])
  seg.dt[, segStart1 := segStart+1]
  sample.counts[, chromStart1 := chromStart+1]
  setkey(seg.dt, segStart1, segEnd)
  setkey(sample.counts, chromStart1, chromEnd)
  over.dt <- foverlaps(sample.counts, seg.dt, nomatch=0L)
  over.dt[chromStart < segStart, chromStart := segStart]
  over.dt[segEnd < chromEnd, chromEnd := segEnd]
  over.dt[, weight := chromEnd-chromStart]
  stopifnot(problemEnd-problemStart == sum(over.dt$weight))
  seg.means <- over.dt[, list(
    mean=sum(coverage*weight)/sum(weight)
  ), by=list(segStart, segEnd)]
  seg.means[, status := rep(c("background", "peak"), l=.N)]
  data.means <- seg.means[over.dt, on=list(segStart, segEnd)]
  logLik.macs <- data.means[, sum(dpois(coverage, mean, log=TRUE)*weight)]
  equality.constraints <- seg.means[, sum(diff(mean) == 0) ]
  total.loss.macs <- data.means[, PeakSegOptimal::PoissonLoss(
    coverage, mean, weight)]
  loss.list[[paste(sample.id, "macs2")]] <- data.table(
    sample.id, 
    model="macs2",
    equality.constraints,
    total.loss=total.loss.macs,
    peaks=nrow(sample.peaks))
  segs.list[[paste(sample.id, "macs2")]] <- data.table(
    sample.id, model="macs2", seg.means)
  for(n.peaks in n.peaks.vec){
    ## better.list <- NULL
    ## while(is.null(better.list)){
    ##   better.list <- tryCatch({
    ##     PeakSegDisk::sequentialSearch_dir(problem.dir, n.peaks)
    ##   }, error=function(e){
    ##     print("trying again")
    ##     NULL
    ##   })
    ## }
    better.list <- PeakSegDisk::sequentialSearch_dir(problem.dir, n.peaks)
    loss.list[[paste(sample.id, n.peaks)]] <- data.table(
      sample.id, 
      model=n.peaks,
      better.list$loss[, .(equality.constraints, total.loss, peaks)])
    segs.list[[paste(sample.id, n.peaks)]] <- data.table(
      sample.id, 
      model=n.peaks,
      better.list$segments[, .(
        segStart=chromStart,
        segEnd=chromEnd,
        mean,
        status
      )]
    )
  }
}
##unlink("jss-figure-more-likely-models", recursive=TRUE)
##unlink("jss-data", recursive=TRUE)

segs <- do.call(rbind, segs.list)
loss <- do.call(rbind, loss.list)
logLik <- segs[, {
  data.table(
    PoissonLogLik(sample.counts, .SD),
    peaks=(.N-1)/2
  )
}, by=list(sample.id, model)]
changes <- segs[order(segStart), data.table(
  position=segStart[-1],
  constraint=ifelse(diff(mean)==0, "equality", "inequality")
), by=list(sample.id, model)]

max.dt <- some.counts[, list(
  max=max(coverage)
), by=list(sample.id)]
max.logLik <- logLik[max.dt, on=list(sample.id)]

possible.vec <- c("macs2", n.peaks.vec)
lab.vec <- gsub(" ", "\n", ifelse(
  possible.vec=="macs2",
  "Default MACS2 model",
  paste0("Optimal ", possible.vec, "-peak model")))
mfactor <- function(val){
  factor(
    val,
    possible.vec,
    lab.vec)
}
max.logLik[, model.fac := mfactor(model)]
segs[, model.fac := mfactor(model)]
changes[, model.fac := mfactor(model)]
gg <- ggplot()+
  geom_text(aes(
    118120, max, label=sprintf(
      "logLik=%.1f\n%d peak%s",
      total.logLik, peaks, ifelse(peaks==1, "", "s"))),
    hjust=1,
    vjust=1,
    color="deepskyblue",
    data=max.logLik)+
  geom_step(aes(
    chromEnd/1e3, coverage),
    color="grey",
    data=some.counts)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sample.id ~ model.fac, scales="free", labeller=function(df){
    if("sample.id" %in% names(df)){
      df$sample.id <- sub("McGill0", "", df$sample.id)
    }else{
      df$model.fac <- paste(df$model.fac)
    }
    df
  })+
  geom_segment(aes(
    segStart/1e3, 0,
    xend=segEnd/1e3, yend=0),
    data=segs[status=="peak"],
    color="deepskyblue",
    size=1)+
  geom_point(aes(
    segStart/1e3, 0),
    data=segs[status=="peak"],
    shape=1,
    color="deepskyblue")+
  geom_segment(aes(
    segStart/1e3, mean,
    xend=segEnd/1e3, yend=mean),
    color="green",
    size=0.5,
    data=segs)+
  xlab("position on chromosome (kb)")+
  ylab("aligned read coverage")

one <- function(dt){
  dt[sample.id=="McGill0322"]
}
gg <- ggplot()+
  geom_text(aes(
    118120000, max, label=sprintf(
      "logLik=%.1f\n%d peak%s",
      total.logLik, peaks, ifelse(peaks==1, "", "s"))),
    hjust=1,
    vjust=1,
    color="deepskyblue",
    data=one(max.logLik))+
  geom_step(aes(
    chromEnd, coverage),
    color="grey",
    data=one(some.counts))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(model.fac ~ .)+
  geom_segment(aes(
    segStart, 0,
    xend=segEnd, yend=0),
    data=one(segs[status=="peak"]),
    color="deepskyblue",
    size=1)+
  geom_point(aes(
    segStart, 0),
    data=one(segs[status=="peak"]),
    shape=1,
    color="deepskyblue")+
  geom_segment(aes(
    segStart, mean,
    xend=segEnd, yend=mean),
    color="green",
    size=1,
    data=one(segs))+
  geom_text(aes(
    segStart, mean, label=format(segStart, big.mark=",")),
    color="green",
    hjust=1,
    data=one(segs[status=="peak"]))+
  geom_text(aes(
    segEnd, mean, label=format(segEnd, big.mark=",")),
    color="green",
    hjust=0,
    data=one(segs[status=="peak"]))+
  xlab("position on chromosome")+
  ylab("aligned read coverage")

one <- function(dt){
  dt[sample.id=="McGill0104"]
}
xmin <- 118122000
xmax <- 118124000
gg <- ggplot()+
  geom_text(aes(
    118120000, max, label=sprintf(
      "logLik=%.1f\n%d peak%s",
      total.logLik, peaks, ifelse(peaks==1, "", "s"))),
    hjust=1,
    vjust=1,
    color="deepskyblue",
    data=one(max.logLik))+
  geom_step(aes(
    chromEnd, coverage),
    color="grey",
    data=one(some.counts))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(model.fac ~ .)+
  geom_segment(aes(
    segStart, 0,
    xend=segEnd, yend=0),
    data=one(segs[status=="peak"]),
    color="deepskyblue",
    size=1)+
  geom_point(aes(
    segStart, 0),
    data=one(segs[status=="peak"]),
    shape=1,
    color="deepskyblue")+
  geom_segment(aes(
    segStart, mean,
    xend=segEnd, yend=mean),
    color="green",
    alpha=0.5,
    size=1,
    data=one(segs))+
  xlab("position on chromosome")+
  scale_linetype_manual(values=c(
    equality="solid",
    inequality="dotted"))+
  geom_vline(aes(
    xintercept=position,
    linetype=constraint),
    data=changes,
    alpha=0.5,
    size=0.4,
    color="green")+
  ylab("aligned read coverage")
png("jss-figure-more-likely-models-three-peaks.png",
    units="in", res=300, width=4, height=3)
print(gg+guides(linetype="none"))
dev.off() 

gg.zoom <- gg+coord_cartesian(xlim=c(xmin, xmax))+
  ggtitle("Zoom to right peak")+
  scale_x_continuous(breaks=NULL)+
  theme(
    legend.position="bottom",
    text=element_text(size=8))
png("jss-figure-more-likely-models-three-peaks-zoom.png",
    units="in", res=300, width=2.5, height=3)
print(gg.zoom)
dev.off() 

