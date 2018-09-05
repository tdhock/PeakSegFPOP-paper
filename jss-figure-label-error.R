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
small.counts <- some.counts[, approx(chromStart, coverage, seq(min(chromStart), max(chromStart), l=1000))]

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
gg <- ggplot()+
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
  "no errors, max peaks"=12,
  ##"min incorrect labels"=15,
  "no errors, min peaks"=18,
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
  meta <- data.table(
    pen.str=paste(round(model$penalty)),
    model.name, model.y=-i*3, penalty=model$penalty, peaks=model$peaks)
  peak.dt.list[[model.name]] <- data.table(meta, fit.peaks)
  error.dt.list[[model.name]] <- data.table(meta, fit.errors)
  ##system(paste("gzip", file.path(prob.dir, coverage.bedGraph.gz)))
}
peak.dt <- do.call(rbind, peak.dt.list)
error.dt <- do.call(rbind, error.dt.list)

coverage.dt <- fread(file.path(prob.dir, "coverage.bedGraph"), drop=c(1,3))
setnames(coverage.dt, c("chromStart", "coverage"))
xlim.vec <- c(1.52e7, 1.73e7)
small.counts <- coverage.dt[, approx(chromStart, coverage, seq(xlim.vec[1], xlim.vec[2], l=1000))]
coverage.dt[xlim.vec[1] < chromStart & chromStart < xlim.vec[2]]

show.peaks <- peak.dt[min(error.dt$chromStart) < chromStart & chromStart < max(error.dt$chromEnd)]
gg.data <- ggplot()+
  theme_bw()+
  scale_fill_manual("label", values=ann.colors)+
  geom_line(aes(
    x, y),
    color="grey50",
    data=small.counts)+
  coord_cartesian(
    xlim=xlim.vec, expand=FALSE,
    ylim=c(-14, 31))+
  xlab("position on chromosome")+
  scale_y_continuous(
    "aligned read counts",
    breaks=seq(0, 30, by=10))
pdf("jss-figure-label-error-data.pdf", 8, 3)
print(gg.data)
dev.off()

gg.labels <- ggplot()+
  theme_bw()+
  scale_fill_manual("label", values=ann.colors)+
  geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd, fill=annotation),
    alpha=0.5,
    size=0.4,
    color="grey",
    data=some.regions[annotation!="noPeaks"])+
  geom_line(aes(
    x, y),
    color="grey50",
    data=small.counts)+
  coord_cartesian(
    xlim=xlim.vec, expand=FALSE,
    ylim=c(-14, 31))+
  xlab("position on chromosome")+
  scale_y_continuous(
    "aligned read counts",
    breaks=seq(0, 30, by=10))
pdf("jss-figure-label-error-data-labels.pdf", 8, 3)
print(gg.labels)
dev.off()

some <- function(dt){
  dt[pen.str=="278653"]
}
gg.few <- ggplot()+
  theme_bw()+
  scale_fill_manual("label", values=ann.colors)+
  geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd, fill=annotation),
    alpha=0.5,
    size=0.4,
    color="grey",
    data=some.regions[annotation!="noPeaks"])+
  geom_line(aes(
    x, y),
    color="grey50",
    data=small.counts)+
  coord_cartesian(
    xlim=xlim.vec, expand=FALSE,
    ylim=c(-14, 31))+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    linetype=status,
    ymin=model.y-1, ymax=model.y+1),
    fill=NA,
    color="black",
    size=0.7,
    data=some(error.dt))+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  geom_segment(aes(
    chromStart, model.y,
    xend=chromEnd, yend=model.y),
    color="deepskyblue",
    size=2,
    data=some(show.peaks))+
  geom_point(aes(
    chromStart, model.y),
    color="deepskyblue",
    shape=1,
    data=some(show.peaks))+
  geom_text(aes(
    chromStart, model.y,
    label=sprintf("penalty=%s ", pen.str)),
    hjust=1,
    data=some(error.dt[chromStart==min(chromStart)]))+
  geom_text(aes(
    chromEnd, model.y,
    label=sprintf(" %s (%d)", model.name, peaks)),
    hjust=0,
    data=some(error.dt[chromStart==max(chromStart)]))+
  xlab("position on chromosome")+
  scale_y_continuous(
    "aligned read counts",
    breaks=seq(0, 30, by=10))
pdf("jss-figure-label-error-too-few.pdf", 8, 3)
print(gg.few)
dev.off()

some <- function(dt){
  dt[pen.str %in% c("278653", "6682") ]
}
gg.many <- ggplot()+
  theme_bw()+
  scale_fill_manual("label", values=ann.colors)+
  geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd, fill=annotation),
    alpha=0.5,
    size=0.4,
    color="grey",
    data=some.regions[annotation!="noPeaks"])+
  geom_line(aes(
    x, y),
    color="grey50",
    data=small.counts)+
  coord_cartesian(
    xlim=xlim.vec, expand=FALSE,
    ylim=c(-14, 31))+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    linetype=status,
    ymin=model.y-1, ymax=model.y+1),
    fill=NA,
    color="black",
    size=0.7,
    data=some(error.dt))+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  geom_segment(aes(
    chromStart, model.y,
    xend=chromEnd, yend=model.y),
    color="deepskyblue",
    size=2,
    data=some(show.peaks))+
  geom_point(aes(
    chromStart, model.y),
    color="deepskyblue",
    shape=1,
    data=some(show.peaks))+
  geom_text(aes(
    chromStart, model.y,
    label=sprintf("penalty=%s ", pen.str)),
    hjust=1,
    data=some(error.dt[chromStart==min(chromStart)]))+
  geom_text(aes(
    chromEnd, model.y,
    label=sprintf(" %s (%d)", model.name, peaks)),
    hjust=0,
    data=some(error.dt[chromStart==max(chromStart)]))+
  xlab("position on chromosome")+
  scale_y_continuous(
    "aligned read counts",
    breaks=seq(0, 30, by=10))
pdf("jss-figure-label-error-too-many.pdf", 8, 3)
print(gg.many)
dev.off()

gg.data <- ggplot()+
  theme_bw()+
  scale_fill_manual("label", values=ann.colors)+
  geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd, fill=annotation),
    alpha=0.5,
    size=0.4,
    color="grey",
    data=some.regions[annotation!="noPeaks"])+
  geom_line(aes(
    x, y),
    color="grey50",
    data=small.counts)+
  coord_cartesian(
    xlim=xlim.vec, expand=FALSE,
    ylim=c(-14, 31))+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    linetype=status,
    ymin=model.y-1, ymax=model.y+1),
    fill=NA,
    color="black",
    size=0.7,
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
    data=show.peaks)+
  geom_point(aes(
    chromStart, model.y),
    color="deepskyblue",
    shape=1,
    data=show.peaks)+
  geom_text(aes(
    chromStart, model.y,
    label=sprintf("penalty=%s ", pen.str)),
    hjust=1,
    data=error.dt[chromStart==min(chromStart)])+
  geom_text(aes(
    chromEnd, model.y,
    label=sprintf(" %s (%d)", model.name, peaks)),
    hjust=0,
    data=error.dt[chromStart==max(chromStart)])+
  xlab("position on chromosome")+
  scale_y_continuous(
    "aligned read counts",
    breaks=seq(0, 30, by=10))
pdf("jss-figure-label-error.pdf", 8, 3)
print(gg)
dev.off()
