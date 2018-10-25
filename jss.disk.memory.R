source("jss-packages.R")

bench.models <- fread("jss.bench.models.csv")
bench.models[, minutes := seconds/60]
bench.models[, hours := minutes/60]
bench.models[, gigabytes := megabytes/1024]
bench.models[, sum(hours)/24/30]#months of computation time

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

data.dir <- "data"
if(!dir.exists(data.dir)){
  download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00439/peak-detection-data.tar.xz", "peak-detection-data.tar.xz")
  system("tar xvf peak-detection-data.tar.xz")
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
      disk=fit <- PeakSegDisk::PeakSegFPOP_disk(
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

saveRDS(fit.dt, "jss.disk.memory.rds")
