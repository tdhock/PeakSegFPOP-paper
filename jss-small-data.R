library(data.table)
library(ggplot2)
chunk.name <- "H3K4me3_TDH_immune/5"
chunk.dir <- file.path("data", chunk.name)
counts.RData <- file.path(chunk.dir, "counts.RData")
load(counts.RData)
counts.dt <- data.table(counts)

sample.num.vec <- c(101, 104, 4, 91, 322)
sample.id.vec <- sprintf("McGill%04d", sample.num.vec)
some.counts <- counts.dt[sample.id %in% sample.id.vec]

peaks.RData <- file.path(chunk.dir, "peaks", "macs.default.RData")
load(peaks.RData)
peaks.dt <- data.table(peaks[[1]])
some.peaks <- peaks.dt[sample.id %in% sample.id.vec]

problem.findPeaks <- function
### Compute a model with better (more likely) peaks, with at most the
### number of peaks given by peaks.int. This function repeated calls
### problem.PeakSegFPOP with different penalty values, until either
### (1) it finds the peaks.int model, or (2) it concludes that there
### is no peaks.int model, in which case it returns the next simplest
### model (with fewer peaks than peaks.int).
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  getCandidateIndex,
### function: how to choose next/optimal model?
  verbose=0
### Print messages?
 ){
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  stopifnot(is.function(getCandidateIndex))
  problem.coverage(problem.dir)
  model.list <- list()
  next.pen <- c(0, Inf)
  iteration <- 0
  while(length(next.pen)){
    if(verbose)cat(
      "Next =", paste(next.pen, collapse=", "),
      "mc.cores=", getOption("mc.cores"),
      "\n")
    next.str <- paste(next.pen)
    iteration <- iteration+1
    model.list[next.str] <- mclapply.or.stop(next.str, function(penalty.str){
      problem.PeakSegFPOP(problem.dir, penalty.str)
    })
    loss.list <- lapply(model.list, "[[", "loss")
    loss.dt <- do.call(rbind, loss.list)[order(penalty)]
    if(!is.numeric(loss.dt$penalty)){
      stop("penalty column is not numeric -- check loss in _loss.tsv files")
    }
    unique.peaks <- loss.dt[, data.table(
      .SD[which.min(iteration)],
      penalties=.N
    ), by=list(peaks)]
    path.dt <- data.table(penaltyLearning::modelSelection(
      unique.peaks, "total.cost", "peaks"))
    path.dt[, next.pen := max.lambda]
    path.dt[, already.computed := next.pen %in% names(loss.list)]
    path.dt[, no.next := c(diff(peaks) == -1, NA)]
    path.dt[, done := already.computed | no.next]
    candidate.i <- getCandidateIndex(path.dt)
    if(!(
      is.integer(candidate.i) &&
      length(candidate.i)==1 &&
      1 <= candidate.i && candidate.i <= nrow(path.dt)-1
    )){
      str(candidate.i)
      stop(
        "getCandidateIndex must return integer scalar ",
        "between 1 and N-1 ",
        "where N is the number of selectable models")
    }
    candidate <- path.dt[candidate.i]
    next.pen <- if(!candidate$done)candidate$next.pen
  }#while(!is.null(pen))
  out <- model.list[[paste(candidate$penalty)]]
  out$others <- path.dt
  out
### TODO
}

problem.betterPeaks <- function
### TODO
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  peaks.int,
### integer scalar: number of peaks to find.
  verbose=0
### Print messages?
){
  stopifnot(
    is.integer(peaks.int) &&
    length(peaks.int)==1 &&
    0 <= peaks.int)
  problem.findPeaks(problem.dir, function(dt){
    dt[, max(which(peaks.int <= peaks))]
  }, verbose)
}

problem.fewerPeaks <- function
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  total.cost.upper.bound,
### numeric scalar: upper bound on loss (we want to find the model
### with fewest peaks, and loss which is below this quantity).
  verbose=0
### Print messages?
){
  stopifnot(
    is.numeric(total.cost.upper.bound) &&
    length(total.cost.upper.bound)==1)
  problem.findPeaks(problem.dir, function(dt){
    dt[, max(which(total.cost <= total.cost.upper.bound))]
  }, verbose)
}

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
for(sample.i in seq_along(sample.id.vec)){
  sample.id <- sample.id.vec[[sample.i]]
  cat(sprintf("%4d / %4d %s\n", sample.i, length(sample.id.vec), sample.id))
  s.dt <- data.table(sample.id)
  sample.counts <- counts.dt[s.dt, on=list(sample.id)]
  sample.peaks <- peaks.dt[s.dt, on=list(sample.id)]
  problemStart <- sample.counts$chromStart[1]
  problemEnd <- sample.counts[.N, chromEnd]
  problem <- data.table(chrom, problemStart, problemEnd)
  problem.dir <- file.path("jss-data", sample.id, sprintf("%s:%d-%d", chrom, problemStart, problemEnd))
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
  infeasible.changes <- seg.means[, sum(sign(diff(mean)) != c(1, -1))]
  total.cost.macs <- data.means[, PeakSegOptimal::PoissonLoss(coverage, mean, weight)]
  ##better.list <- problem.betterPeaks(problem.dir, nrow(sample.peaks))
  better.list <- problem.fewerPeaks(problem.dir, total.cost.macs)
  loss.list[[paste(sample.id)]] <- data.table(
    sample.id, rbind(
      data.table(
        model="macs2",
        status=ifelse(infeasible.changes==0, "feasible", "infeasible"),
        total.cost=total.cost.macs,
        peaks=nrow(sample.peaks)),
      data.table(
        model="better",
        better.list$loss[, .(status, total.cost, peaks)])))
  segs.list[[paste(sample.id)]] <- sample.segs <- data.table(
    sample.id, rbind(
      data.table(model="macs2", seg.means),
      data.table(
        model="better",
        better.list$segments[, .(segStart=chromStart, segEnd=chromEnd, mean, status)])))
  logLik.list[[paste(sample.id)]] <- sample.segs[, {
    data.table(
      PoissonLogLik(sample.counts, .SD),
      peaks=(.N-1)/2
      )
  }, by=list(sample.id, model)]
}
logLik <- do.call(rbind, logLik.list)
segs <- do.call(rbind, segs.list)
loss <- do.call(rbind, loss.list)

max.dt <- some.counts[, list(
  max=max(coverage)
), by=list(sample.id)]
max.logLik <- logLik[max.dt, on=list(sample.id)]

ggplot()+
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
  facet_grid(sample.id ~ model, scales="free")+
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
    size=1,
    data=segs)


