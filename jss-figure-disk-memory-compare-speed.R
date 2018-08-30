library(data.table)
library(PeakSegOptimal)
library(PeakSegPipeline)
library(ggplot2)
library(microbenchmark)
library(directlabels)

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)][log(bedGraph.lines) < penalty & penalty < bedGraph.lines & 1000 < bedGraph.lines]
bench.models[, minutes := seconds/60]
bench.models[, hours := minutes/60]
bench.models[, gigabytes := megabytes/1024]
gigabyte.ranges <- bench.models[0 < gigabytes, list(
  min.gigabytes=min(gigabytes),
  max.gigabytes=max(gigabytes),
  models=.N
  ), by=list(bedGraph.lines, prob.dir)]
small.probs <- gigabyte.ranges[max.gigabytes < 1 & 10 < models][order(bedGraph.lines)]
log10.range <- small.probs[, log10(range(bedGraph.lines))]
some.probs <- data.table(target.N=10^seq(log10.range[1], log10.range[2], l=10))[, {
  small.probs[which.min(abs(target.N-bedGraph.lines))]
}, by=list(target.N)]
##some.probs <- small.probs[round(seq(1, .N, l=10))]
small.models <- bench.models[some.probs, on=list(prob.dir, bedGraph.lines)]
small.models[, sum(hours)]

small.models[, plot( log10(penalty) ~ log10(bedGraph.lines) )]

## ggplot()+
##   geom_point(aes(
##     bedGraph.lines, penalty),
##     data=small.models)+
##   scale_x_log10()+
##   scale_y_log10()

bench.dir <- "~/projects/feature-learning-benchmark"
data.dir <- file.path(bench.dir, "data")
if(!file.exists(data.dir)){
  download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00439/peak-detection-data.tar.xz", "data.tar.xz")
  system(paste("cd", bench.dir, "&& tar xvf data.tar.xz && rm data.tar.xz"))
}

fit.dt.list <- list()
zip.dt.list <- list()
model.i.vec <- 1:nrow(small.models)
##model.i.vec <- 1:151
for(model.i in model.i.vec){
  model <- small.models[model.i]
  pen.str <- paste(model$penalty)
  times.RData <- file.path(
    data.dir, model$prob.dir, paste0("times_pen=", pen.str, ".RData"))
  if(file.exists(times.RData)){
    load(times.RData)
  }else{
    cat(sprintf("%4d / %4d models\n", model.i, nrow(small.models)))
    coverage.bedGraph.gz <- file.path(
      data.dir, model$prob.dir, "coverage.bedGraph.gz")
    gunzip.seconds <- system.time({
      system(paste("gunzip", coverage.bedGraph.gz))
    })[["elapsed"]]
    coverage.bedGraph <- file.path(
      data.dir, model$prob.dir, "coverage.bedGraph")
    fread.seconds <- system.time({
      coverage.dt <- fread(coverage.bedGraph)
      setnames(coverage.dt, c("chrom", "chromStart", "chromEnd", "count"))
    })[["elapsed"]]
    time.df <- microbenchmark(
      memory=PeakSegOptimal::PeakSegFPOPchrom(coverage.dt, model$penalty),
      disk=fit <- PeakSegPipeline::PeakSegFPOP_disk(
        coverage.bedGraph, pen.str),
      times=2)
    unlink(fit$db)
    gzip.seconds <- system.time({
      system(paste("gzip", coverage.bedGraph))
    })[["elapsed"]]
    zip.times <- data.table(gunzip.seconds, gzip.seconds, fread.seconds)
    fit.times <- data.table(time.df)
    save(zip.times, fit.times, file=times.RData)
  }
  zip.dt.list[[model.i]] <- data.table(model, zip.times)
  fit.dt.list[[model.i]] <- data.table(model, fit.times)
}
zip.dt <- do.call(rbind, zip.dt.list)
fit.dt <- do.call(rbind, fit.dt.list)

fit.dt[, bench.seconds := time/1e9 ]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("bedGraph.lines")+
  geom_point(aes(
    log(penalty), bench.seconds, color=expr),
    shape=1,
    data=fit.dt)+
  scale_y_log10("seconds")

ggplot()+
  geom_point(aes(
    log10(bedGraph.lines), bench.seconds, color=expr),
    shape=1,
    data=fit.dt)+
  scale_y_log10("seconds")

bench.stats <- fit.dt[, list(
  median=median(bench.seconds),
  q25=quantile(bench.seconds, 0.25),
  q75=quantile(bench.seconds, 0.75)
), by=list(bedGraph.lines, storage=expr)]

wide <- dcast(bench.stats, bedGraph.lines ~ storage, value.var="median")

ggplot()+
  geom_point(aes(
    log10(bedGraph.lines), disk/memory),
    data=wide)

details.dt <- fit.dt[bedGraph.lines==106569]
leg <- ggplot()+
  geom_point(aes(
    log(penalty), bench.seconds, color=expr),
    shape=1,
    data=details.dt)+
  scale_y_log10("seconds")+
  xlim(NA, 13)
dl <- direct.label(leg, list("last.qp", dl.trans(x=x+0.1)))
pdf("jss-figure-disk-memory-compare-speed-penalty.pdf", 3.3, 3)
print(dl)
dev.off()

range.dt <- details.dt[, list(
  min.seconds=min(bench.seconds),
  max.seconds=max(bench.seconds),
  bedGraph.lines=bedGraph.lines[1]
  )]
leg <- ggplot()+
  geom_segment(aes(
    log10(bedGraph.lines), min.seconds,
    xend=log10(bedGraph.lines), yend=max.seconds),
    data=range.dt)+
  geom_ribbon(aes(
    log10(bedGraph.lines), ymin=q25, ymax=q75, fill=storage),
    alpha=0.5,
    data=bench.stats)+
  geom_line(aes(
    log10(bedGraph.lines), median, color=storage),
    data=bench.stats)+
  scale_y_log10("seconds")+
  scale_x_continuous(
    "log10(N = number of data to segment)",
    limits=c(NA, 6.2))
dl <- direct.label(leg, "last.polygons")
pdf("jss-figure-disk-memory-compare-speed.pdf", 3.3, 3)
print(dl)
dev.off()
