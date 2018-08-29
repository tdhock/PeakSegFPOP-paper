library(data.table)
library(PeakSegOptimal)
library(PeakSegPipeline)

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)]
bench.models[, minutes := seconds/60]
bench.models[, hours := minutes/60]
bench.models[, gigabytes := megabytes/1024]
gigabyte.ranges <- bench.models[0 < gigabytes, list(
  min.gigabytes=min(gigabytes),
  max.gigabytes=max(gigabytes)
  ), by=list(bedGraph.lines, prob.dir)]
small.probs <- gigabyte.ranges[max.gigabytes < 0.5][order(max.gigabytes)]
small.models <- bench.models[small.probs, on=list(prob.dir)]
small.models[, sum(hours)]

bench.dir <- "~/projects/feature-learning-benchmark"
data.dir <- file.path(bench.dir, "data")
if(!file.exists(data.dir)){
  download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00439/peak-detection-data.tar.xz", "data.tar.xz")
  system(paste("cd", bench.dir, "&& tar xvf data.tar.xz && rm data.tar.xz"))
}


