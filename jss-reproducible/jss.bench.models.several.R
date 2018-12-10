source("jss-packages.R")

max.days <- 2

bench.models <- fread("jss.bench.models.csv")[order(bedGraph.lines, penalty)]
bench.models[, row := 1:.N]
bench.models[, minutes := seconds/60]
bench.models[, hours := minutes/60]
bench.models[, days := hours/24]
out.csv.vec <- dir("jss.bench.models.rules")
bench.models[, done := paste0(row, ".csv") %in% out.csv.vec]
(problems <- bench.models[done==FALSE, list(
  total.days=sum(days)
), by=list(prob.dir)][order(-total.days)])
problems[, job := rep(1:min(.N, 500), l=.N)]

## set.seed(1)
## (rand.models <- problems[order(sample(1:.N))])
## rand.models[, cum.days := cumsum(total.days)]
## rand.models[, job := (cum.days %/% max.days)+1]
## rand.models[, remainder := cum.days %% max.days]
## rand.models[, rank := rank(remainder), by=job]
## rand.models[rank==1, job := job-1]
## rand.models[, list(
##   probs=.N,
##   min.days=min(total.days),
##   max.days=max(total.days),
##   total.days=sum(total.days)
## ), by=list(job)][order(total.days)]

problems[, list(
  probs=.N,
  min.days=min(total.days),
  max.days=max(total.days),
  total.days=sum(total.days)
), by=list(job)][order(-total.days)]

join.dt <- bench.models[problems, on=list(prob.dir)]

fwrite(join.dt[order(-job, row), list(row, job)], "jss.bench.models.several.csv")
