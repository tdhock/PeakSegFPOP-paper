source("findPeaks.R")
library(data.table)
library(PeakSegPipeline)

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)][log(bedGraph.lines) < penalty & penalty < bedGraph.lines & 1000 < bedGraph.lines]
bench.models[, gigabytes := megabytes/1024]

min.err <- bench.models[, list(
  n.feasible.0=sum(errors==0 & status=="feasible"),
  min.errors=min(errors),
  max.fp=max(fp),
  max.fn=max(fn),
  max.gigabytes=max(gigabytes)
  ), by=list(bedGraph.lines, prob.dir)]

zero.err <- min.err[min.errors==0 & max.gigabytes < 30 & 0 < max.fp & 0 < max.fn][order(bedGraph.lines)]
picked.models <- bench.models[zero.err, {
  .SD[errors==0, list(peaks=as.integer((min(peaks)+max(peaks))/2))]
}, by=list(prob.dir, bedGraph.lines), on=list(bedGraph.lines, prob.dir)]

picked.models[, log10.peaks := log10(peaks)]
picked.models[, log10.lines := log10(bedGraph.lines)]
fit <- lm(log10.peaks ~ log10.lines, picked.models)
log10.range <- picked.models[, log10(range(bedGraph.lines))]
line.dt <- data.table(log10.lines=seq(log10.range[1], log10.range[2], l=100))
line.dt[, log10.peaks := predict(fit, line.dt)]

picked.models[, residual := log10.peaks-predict(fit, picked.models)]
close.models <- picked.models[abs(residual) < 0.1]

some.probs <- data.table(target.N=10^seq(log10.range[1], log10.range[2], l=10))[, {
  close.models[order(abs(target.N-bedGraph.lines))][1:2]
}, by=list(target.N)][order(bedGraph.lines)]

prob.i.vec <- c(1, 2, 7, 8, 13, 14)
jss.variable.peaks.list <- list()
for(prob.i in prob.i.vec){
  prob <- some.probs[prob.i]
  pdir <- file.path("~/projects/feature-learning-benchmark/data", prob$prob.dir)
  system(paste("gunzip", file.path(pdir, "coverage.bedGraph.gz")))
  match.df <- namedCapture::str_match_named(pdir, paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<problemStart>[0-9]+)",
    "-",
    "(?<problemEnd>[0-9]+)"), list(
      problemStart=as.integer,
      problemEnd=as.integer))
  fwrite(
    match.df, file.path(pdir, "problem.bed"),
    quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
  trivial <- problem.betterPeaks(pdir, 0L, verbose=1)
  most.peaks <- trivial$others$peaks[1]
  most.peaks <- prob$target.N/2
  peak.vec <- as.integer(c(
    1:9,
    10^seq(1, log10(most.peaks), by=0.5)))
  for(peak.i in seq_along(peak.vec)){
    peaks <- peak.vec[[peak.i]]
    cat(sprintf(
      "%4d / %4d problems %4d / %4d peaks=%d\n",
      prob.i, length(prob.i.vec),
      peak.i, length(peak.vec),
      peaks
    ))
    fit.list <- problem.betterPeaks(pdir, peaks, verbose=1)
    jss.variable.peaks.list[[paste(pdir, peaks)]] <- data.table(
      prob,
      loss=fit.list$loss,
      others=fit.list$others)
  }
}
jss.variable.peaks <- do.call(rbind, jss.variable.peaks.list)

saveRDS(jss.variable.peaks, "jss.variable.peaks.rds")
