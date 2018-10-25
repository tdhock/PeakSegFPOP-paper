library(data.table)
library(ggplot2)

load("Segmentor.peaks.error.RData")
load("Segmentor.infeasible.error.RData")
load("PDPA.infeasible.error.RData")
load("PDPA.peaks.error.RData")
load("dp.peaks.error.RData")

CDPA.error.list <- list()
for(chunk.name in names(dp.peaks.error)){
  chunk.error <- dp.peaks.error[[chunk.name]]
  for(cell.type in names(chunk.error)){
    type.dt <- data.table(
      chunk.name,
      chunk.error[[cell.type]])
    names(type.dt)[3] <- "peaks"
    CDPA.error.list[[paste(chunk.name, cell.type)]] <- type.dt
  }
}
CDPA.error <- do.call(rbind, CDPA.error.list)

CDPA.error[, segments := as.integer(paste(peaks))*2+1]
PDPA.peaks.error[, segments := as.integer(paste(peaks))*2+1]
Segmentor.peaks.error[, segments := as.integer(paste(peaks))*2+1]
col.name.vec <- names(PDPA.infeasible.error)
DT <- function(algo, dt){
  data.table(algo, dt[, ..col.name.vec])
}
all.error <- rbind(
  DT("CDPA", CDPA.error),
  DT("G.ig", PDPA.peaks.error),
  DT("G.jo", PDPA.infeasible.error),
  DT("S.ig", Segmentor.peaks.error),
  DT("S.jo", Segmentor.infeasible.error))

all.error[, table(chunk.name, algo)]

all.totals <- all.error[, list(
  total.fp=sum(fp),
  total.fn=sum(fn),
  total.errors=sum(fp+fn),
  possible.fp=sum(possible.fp),
  possible.fn=sum(possible.tp),
  labels=.N
  ), by=list(algo, chunk.name, sample.id, peaks, segments)]

all.totals[, table(chunk.name, algo)]

all.min <- all.totals[, list(
  min.errors=min(total.errors)
  ), by=list(algo, chunk.name, sample.id)]

all.min[, table(chunk.name, algo)]

min.wide <- dcast(all.min, chunk.name+sample.id~algo)

## as expected there are some (81) problems where the peak joining
## achieves fewer errors than ignoring the infeasible models with
## equality constraints.
min.wide[G.jo < G.ig]
min.wide[G.ig < G.jo]

## Figure for illustrating infeasible.
some.diff <- min.wide[G.jo==0 & 0<G.ig & 0<CDPA]
big.diff <- min.wide[G.jo==0 & 2==G.ig & 2==CDPA]
small.peaks <- dcast(
  min.models[some.diff, on=list(chunk.name, sample.id)],
  chunk.name+sample.id~algo,
  value.var="range")[order(G.jo)][6:10]
big.diff.tall <- all.totals[small.peaks, on=list(chunk.name, sample.id)]
(big.diff.wide <- dcast(
  big.diff.tall,
  chunk.name+sample.id+peaks~algo,
  value.var="total.errors"))
##24: H3K4me3_TDH_immune/11 McGill0009     3    2          NA         0
chunk.name <- "H3K4me3_TDH_immune/11"
sample.id <- "McGill0009"
peaks <- 3

## Somewhat strangely there are 35 problems when CDPA gets fewer label
## errors, and also 35 problems when GPDPA with join when infeasible
## gets fewer label errors.
min.wide[CDPA < G.jo]
min.wide[G.jo < CDPA]

## Even more strangely the distribution of differences is symmetrical!
min.wide[, table(CDPA-G.jo)]
min.wide[, set.name := sub("/.*", "", chunk.name)]
min.wide[, table(set.name, CDPA-G.jo)]
## > min.wide[, table(diff)]
## diff
##   -2   -1    0    1    2 
##    3   32 2682   32    3 
## >

min.wide[, list(
  count=.N),
  by=list(GPDPA.jo.ig=G.jo-G.ig)]

## ignore vs join
g <- function(x, Comparison){
  dt <- data.table(table(x), Comparison)
  names(dt)[1] <- "diff"
  dt
}
count.tall <- min.wide[, rbind(
  g(G.jo-G.ig, "Join-ignore GPDPA"),
  g(S.jo-S.ig, "Join-ignore PDPA"),
  g(G.jo-S.jo, "GPDPA-PDPA Join"),
  g(G.ig-S.ig, "GPDPA-PDPA Ignore"),
  g(CDPA-G.jo, "CDPA-GPDPAjoin"))]
count.tall[, diff.fac := factor(diff, (-10):2)]
(count.wide <- dcast(count.tall, Comparison ~diff.fac, value.var="N"))
library(xtable)
xt <- xtable(count.wide, align="lrrrrrrrrrr", digits=0)
print(
  xt,
  file="PDPA-infeasible-error-compare.tex",
  include.rownames=FALSE, floating=FALSE)

## compare total min error.
[min.by.sample <- all.totals[, {
  .SD[which.min(total.errors)]
}, by=list(algo, chunk.name, sample.id)]
min.by.sample[, set.name := sub("/.*", "", chunk.name)]

min.by.sample[, list(
  error=sum(total.errors)/sum(labels)*100,
  fp=sum(total.fp)/sum(possible.fp)*100,
  fn=sum(total.fn)/sum(possible.fn)*100
  ), by=list(algo)][order(error)]

min.by.sample[, list(
  Errors=sum(total.errors),
  FP=sum(total.fp),
  FN=sum(total.fn)
  ), by=list(algo)][order(Errors)]

algo.map <- c(
  CDPA="Constrained, approximate (CDPA)",
  G.jo="Constrained, optimal (GPDPA), join equality constraints",
  G.ig="Constrained, optimal (GPDPA), ignore equality constraints",
  S.jo="Unconstrained, optimal (PDPA), join",
  S.ig="Unconstrained, optimal (PDPA), ignore")
min.totals <- min.by.sample[, list(
  `Error Rate (%)`=sum(total.errors)/sum(labels)*100,
  `False Positive Rate (%)`=sum(total.fp)/sum(possible.fp)*100,
  `False Negative Rate (%)`=sum(total.fn)/sum(possible.fn)*100
  ), by=list(Algorithm=algo.map[algo])][order(`Error Rate (%)`)]
set.totals <- min.by.sample[, list(
  error=sum(total.errors)/sum(labels)*100,
  fp=sum(total.fp)/sum(possible.fp)*100,
  fn=sum(total.fn)/sum(possible.fn)*100
), by=list(set.name, algo)][order(set.name, error)]

set.totals[, list(
  mean.error=mean(error),
  sd.error=sd(error)
  ), by=list(algo)][order(mean.error)]


## How many of these problems have possible.fp > 0?
nonzero.fp <- all.totals[algo=="G.jo" & peaks==0 & 0 < possible.fp, list(
  chunk.name, sample.id)]

## For hundreds of problems, the model with 9 peaks has the fewest
## label errors -- this suggests increasing the max number of
## peaks. 
min.models <- all.totals[, {
  .SD[total.errors==min(total.errors), {
    p <- as.integer(paste(peaks))
    list(
      m=min(p),
      M=max(p)
    )
  }]
}, by=list(algo, chunk.name, sample.id)]
min.models[, a := substr(algo, 1, 1)]
min.models[, range := ifelse(
  m==M, m, paste0(m, "-", M))]
models.wide <- dcast(
  min.models[algo %in% c("G.jo", "CDPA")],
  chunk.name+sample.id ~ algo,
  value.var="range")
models.wide[grepl("9", G.jo)]
models.wide[grepl("9", CDPA)]

## take a look at the models where CDPA is better.
load("PDPA.infeasible.RData")
load("dp.peaks.RData")
CDPA.diff <- min.wide[diff != 0][order(diff)]
models.wide[CDPA.diff[, list(chunk.name, sample.id, diff)], on=list(chunk.name, sample.id)]
some.error <- all.error[algo %in% c("CDPA", "G.jo")]
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
some.totals <- all.totals[algo %in% c("CDPA", "G.jo")]
for(problem.i in 1:nrow(CDPA.diff)){
  problem <- CDPA.diff[problem.i]
  problem.totals <- some.totals[problem, on=list(chunk.name, sample.id)]
  problem.wide <- dcast(
    problem.totals,
    peaks ~ algo,
    value.var=c("total.fp", "total.fn"))
  print(problem)
  print(problem.wide)
  problem.error <- some.error[problem, on=list(chunk.name, sample.id)]
  problem.error[, segments := as.integer(paste(peaks))*2+1]
  problem.totals[, segments := as.integer(paste(peaks))*2+1]
  problem.totals[, is.min := total.errors==min(total.errors), by=list(algo)]
  pre <- paste0(
    "../chip-seq-paper/chunks/",
    problem$chunk.name)
  counts.RData <- paste0(
    pre,
    "/counts.RData")
  regions.RData <- paste0(
    pre,
    "/regions.RData")
  load(counts.RData)
  load(regions.RData)
  counts.dt <- data.table(counts)[sample.id==problem$sample.id]
  regions.dt <- data.table(regions)[sample.id==problem$sample.id]
  chunk.PDPA <- PDPA.infeasible[[problem$chunk.name]]
  chunk.CDPA <- dp.peaks[[problem$chunk.name]]
  sample.PDPA <- chunk.PDPA[[paste(problem$sample.id)]]
  sample.CDPA <- chunk.CDPA[[paste(problem$sample.id)]]
  sample.peaks <- rbind(
    do.call(rbind, sample.PDPA)[, data.table(
      algo="G.jo", segments, chromStart, chromEnd)],
    data.table(do.call(rbind, sample.CDPA))[, data.table(
      algo="CDPA", segments, chromStart, chromEnd)])
  max.coverage <- max(counts.dt$coverage)
  peak.y <- c(
    CDPA=3,
    G.jo=1)
  getY <- function(a){
    peak.y[a]*-0.1*max.coverage
  }
  h <- 0.03*max.coverage
  sample.peaks[, y := getY(algo)]
  problem.error[, y := getY(algo)]
  problem.totals[, y := getY(algo)]
  text.size <- 3
  gg <- ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(segments ~ .)+
    penaltyLearning::geom_tallrect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3,
      fill=annotation),
      color="grey",
      data=regions.dt)+
    geom_step(aes(
      chromStart/1e3, coverage),
      color="grey50",
      data=counts.dt)+
    scale_fill_manual(values=ann.colors)+
    geom_segment(aes(
      chromStart/1e3, y,
      xend=chromEnd/1e3, yend=y),
      data=sample.peaks,
      size=2,
      color="deepskyblue")+
    scale_linetype_manual(
      "error type",
      values=c(
        correct=0,
        "false negative"=3,
        "false positive"=1))+
    geom_rect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3,
      linetype=status,
      ymin=y-h, ymax=y+h),
      color="black",
      fill=NA,
      data=problem.error)+
    geom_text(aes(
      min(regions.dt$chromStart)/1e3, y,
      color=is.min,
      label=sprintf(
        "%s errors=%d ", algo, total.errors)),
      hjust=1,
      size=text.size,
      data=problem.totals)+
    scale_color_manual(values=c(
      "TRUE"="black",
      "FALSE"="grey"))+
    geom_text(aes(
      max(regions.dt$chromEnd)/1e3, y,
      color=is.min,
      label=sprintf(
        " fp=%d fn=%d", total.fp, total.fn)),
      hjust=0,
      size=text.size,
      data=problem.totals)
  dir.create("PDPA.infeasible.error.compare", showWarnings=FALSE)
  png.name <- problem[, file.path(
    "PDPA.infeasible.error.compare",
    paste0(
      G.jo-CDPA,
      "_",
      sub("/", "_", chunk.name),
      "_",
      sample.id,
      ".png"
    ))]
  png(png.name, 19, 9, units="in", res=100)
  print(gg)
  dev.off()
  system(paste("falkon", png.name))
}
