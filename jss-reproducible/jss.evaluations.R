source("jss-packages.R")

bench.models <- fread("jss.bench.models.csv")
bench.models[, gigabytes := megabytes/1024]

min.err <- bench.models[, list(
  min.errors=min(errors),
  max.fp=max(fp),
  max.fn=max(fn),
  max.gigabytes=max(gigabytes)
  ), by=list(bedGraph.lines, prob.dir)]

zero.err <- min.err[min.errors==0 & max.gigabytes < 30 & 0 < max.fp & 0 < max.fn][order(bedGraph.lines)]
## picking the number of peaks in the middle of the zero-error
## interval -- this makes the sequential search iterations plot looks
## O(log N) like it should be.
picked.models <- bench.models[zero.err, {
  .SD[errors==0, list(peaks=as.integer((min(peaks)+max(peaks))/2))]
}, by=list(prob.dir, bedGraph.lines), on=list(bedGraph.lines, prob.dir)]

picked.models[, log10.peaks := log10(peaks)]
picked.models[, log10.lines := log10(bedGraph.lines)]
fit <- lm(log10.peaks ~ log10.lines, picked.models)
log10.range <- picked.models[, log10(range(bedGraph.lines))]
line.dt <- data.table(log10.lines=seq(log10.range[1], log10.range[2], l=100))
line.dt[, log10.peaks := predict(fit, line.dt)]

gg <- ggplot()+
  geom_point(aes(
    log10.lines, log10.peaks),
    shape=1,
    data=picked.models)+
  geom_line(aes(
    log10.lines, log10.peaks),
    color="red",
    data=line.dt)

picked.models[, residual := log10.peaks-predict(fit, picked.models)]
close.models <- picked.models[abs(residual) < 0.1]

some.probs <- data.table(target.N=10^seq(log10.range[1], log10.range[2], l=10))[, {
  close.models[order(abs(target.N-bedGraph.lines))][1:2]
}, by=list(target.N)][order(bedGraph.lines)]

gg <- ggplot()+
  geom_point(aes(
    log10.lines, log10.peaks),
    shape=1,
    data=picked.models)+
  geom_point(aes(
    log10.lines, log10.peaks),
    color="red",
    data=some.probs)

gg <- ggplot()+
  geom_point(aes(
    bedGraph.lines, peaks),
    data=some.probs)+
  scale_x_log10()+
  scale_y_log10()

data.dir <- "data"
if(!dir.exists(data.dir)){
  download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00439/peak-detection-data.tar.xz", "peak-detection-data.tar.xz")
  system("tar xvf peak-detection-data.tar.xz")
}

prob.i.vec <- 1:nrow(some.probs)
##prob.i.vec <- 1:16
jss.evaluations.list <- list()
for(prob.i in prob.i.vec){
  prob <- some.probs[prob.i]
  cat(sprintf("%4d / %4d problems\n", prob.i, length(prob.i.vec)))
  pdir <- file.path(data.dir, prob$prob.dir)
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
  fit.list <- problem.sequentialSearch(pdir, prob$peaks, verbose=1)
  jss.evaluations.list[[prob.i]] <- data.table(
    prob,
    loss=fit.list$loss,
    others=fit.list$others)
}
jss.evaluations <- do.call(rbind, jss.evaluations.list)

saveRDS(jss.evaluations, "jss.evaluations.rds")

