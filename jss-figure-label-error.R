library(data.table)
library(ggplot2)
library(directlabels)
library(PeakSegPipeline)
library(penaltyLearning)

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)][log(bedGraph.lines) < penalty & penalty < bedGraph.lines & 1000 < bedGraph.lines]
bench.models[prob.dir=="H3K4me3_TDH_immune/samples/tcell/McGill0008/problems/chr11:96437584-134946516"]

some.models <- bench.models[prob.dir=="H3K36me3_TDH_immune/samples/bcell/McGill0322/problems/chr12:7239876-34856694"]
model.i.vec <- c(3, 8, 13, 18)

chunk.name <- "H3K36me3_TDH_immune/4"
chunk.dir <- file.path("data", chunk.name)
counts.RData <- file.path(chunk.dir, "counts.RData")
load(counts.RData)
counts.dt <- data.table(counts)
some.counts <- counts.dt[sample.id == "McGill0104"]

load(file.path(chunk.dir, "regions.RData"))
regions.dt <- data.table(regions)
some.regions <- regions.dt[sample.id=="McGill0104"]

some.counts[, plot(chromStart, coverage)]

small.counts <- some.counts[, approx(chromStart, coverage, seq(min(chromStart), max(chromStart), l=1000))]

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
ggplot()+
  theme_bw()+
  scale_fill_manual(values=ann.colors)+
  geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd, fill=annotation),
    alpha=0.5,
    color="grey",
    data=some.regions)+
  geom_step(aes(
    x, y),
    color="grey50",
    data=small.counts)

(some.models <- bench.models[prob.dir=="H3K36me3_TDH_immune/samples/monocyte/McGill0104/problems/chr12:7239876-34856694"])
model.i.vec <- c(
  "too many peaks"=5,
  "min incorrect labels"=15,
  "too few peaks"=19)

error.dt.list <- list()
peak.dt.list <- list()
for(i in seq_along(model.i.vec)){
  model.name <- names(model.i.vec)[[i]]
  model.i <- model.i.vec[[i]]
  cat(sprintf("%4d / %4d %s\n", model.i, length(model.i.vec), model.name))
  model <- some.models[model.i]
  prob.dir <- file.path(
    "~/projects/feature-learning-benchmark/data",
    model$prob.dir)
  pen.str <- paste(model$penalty)
  system(paste("gunzip", file.path(prob.dir, "coverage.bedGraph.gz")))
  fit <- PeakSegPipeline::problem.PeakSegFPOP(prob.dir, pen.str)
  fit.peaks <- fit$segments[status=="peak"]
  fit.errors <- PeakError::PeakErrorChrom(fit.peaks, some.regions)
  meta <- data.table(model.name, model.y=-i*3, penalty=model$penalty)
  peak.dt.list[[model.name]] <- data.table(meta, fit.peaks)
  error.dt.list[[model.name]] <- data.table(meta, fit.errors)
  ##system(paste("gzip", file.path(prob.dir, coverage.bedGraph.gz)))
}
peak.dt <- do.call(rbind, peak.dt.list)
error.dt <- do.call(rbind, error.dt.list)

ggplot()+
  theme_bw()+
  scale_fill_manual(values=ann.colors)+
  geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd, fill=annotation),
    alpha=0.5,
    color="grey",
    data=some.regions)+
  geom_line(aes(
    x, y),
    color="grey50",
    data=small.counts)+
  coord_cartesian(xlim=range(small.counts$x), expand=FALSE)+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    linetype=status,
    ymin=model.y-1, ymax=model.y+1),
    fill=NA,
    color="black",
    size=1,
    data=error.dt)+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  geom_segment(aes(
    chromStart, model.y,
    xend=chromEnd, yend=model.y),
    color="deepskyblue",
    size=2,
    data=peak.dt)+
  geom_point(aes(
    chromStart, model.y),
    color="deepskyblue",
    shape=1,
    data=peak.dt)
