if(!file.exists("labeled_problems_features.csv")){
  download.file(
    "https://raw.githubusercontent.com/tdhock/feature-learning-benchmark/master/labeled_problems_features.csv",
    "labeled_problems_features.csv")
}
if(!file.exists("target.intervals.models.csv")){
  download.file(
    "http://members.cbio.ensmp.fr/~thocking/data/target.intervals.models.csv",
    "target.intervals.models.csv")
}
system("touch target.intervals.models.csv labeled_problems_features.csv")
source("jss-packages.R")

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
bench.models <-
  target.intervals.models[labeled_problems_features$prob.dir, on=list(
    prob.dir)][log(bedGraph.lines) < penalty &
               penalty < bedGraph.lines &
               1000 < bedGraph.lines]
fwrite(bench.models[, list(
  prob.dir, bedGraph.lines, penalty, megabytes, seconds, peaks,
  bases, mean.pen.cost, total.cost, mean.intervals, max.intervals,
  fn, fp, errors)], "jss.bench.models.csv")
