library(data.table)
library(ggplot2)

load("PDPA.infeasible.error.RData")
load("PDPA.peaks.error.RData")
load("dp.peaks.error.RData")

load("PDPA.infeasible.RData")
load("dp.peaks.RData")

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

all.error <- rbind(
  data.table(algo="CDPA", CDPA.error),
  data.table(algo="PDPA.ignore", PDPA.peaks.error),
  data.table(algo="PDPA.join", PDPA.infeasible.error))

all.error[, table(chunk.name, algo)]

all.totals <- all.error[, list(
  total.fp=sum(fp),
  total.fn=sum(fn),
  total.errors=sum(fp+fn),
  possible.fp=sum(possible.fp),
  possible.fn=sum(possible.tp),
  labels=.N
  ), by=list(algo, chunk.name, sample.id, peaks)]

all.totals[, table(chunk.name, algo)]

all.min <- all.totals[, list(
  min.errors=min(total.errors)
  ), by=list(algo, chunk.name, sample.id)]

all.min[, table(chunk.name, algo)]

min.wide <- dcast(all.min, chunk.name+sample.id~algo)

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
  min.models[algo %in% c("PDPA.join", "CDPA")],
  chunk.name+sample.id ~ algo,
  value.var="range")

## as expected there are some (81) problems where the peak joining
## achieves fewer errors than ignoring the infeasible models with
## equality constraints.
min.wide[PDPA.join < PDPA.ignore]
min.wide[PDPA.ignore < PDPA.join]

## Somewhat strangely there are 35 problems when CDPA gets fewer label
## errors, and also 35 problems when PDPA with join when infeasible
## gets fewer label errors.
min.wide[CDPA < PDPA.join]
min.wide[PDPA.join < CDPA]

## Even more strangely the distribution of differences is symmetrical!
min.wide[, diff := CDPA-PDPA.join]
min.wide[, table(diff)]
## > min.wide[, table(diff)]
## diff
##   -2   -1    0    1    2 
##    3   32 2682   32    3 
## > 

## Not so when ignoring infeasible models -- CDPA is better.
min.wide[, table(CDPA-PDPA.ignore)]

## For hundreds of problems, the model with 9 peaks has the fewest
## label errors -- this suggests increasing the max number of
## peaks. 
models.wide[grepl("9", PDPA.join)]
models.wide[grepl("9", CDPA)]

## How many of these problems have possible.fp > 0?
nonzero.fp <- all.totals[algo=="PDPA.join" & peaks==0 & 0 < possible.fp, list(
  chunk.name, sample.id)]

## take a look at the models where CDPA is better.
CDPA.diff <- min.wide[diff != 0][order(diff)]
models.wide[CDPA.diff[, list(chunk.name, sample.id, diff)], on=list(chunk.name, sample.id)]
some.error <- all.error[algo %in% c("CDPA", "PDPA.join")]
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
some.totals <- all.totals[algo %in% c("CDPA", "PDPA.join")]
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
      algo="PDPA.join", segments, chromStart, chromEnd)],
    data.table(do.call(rbind, sample.CDPA))[, data.table(
      algo="CDPA", segments, chromStart, chromEnd)])
  max.coverage <- max(counts.dt$coverage)
  peak.y <- c(
    CDPA=3,
    PDPA.join=1)
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
      PDPA.join-CDPA,
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
