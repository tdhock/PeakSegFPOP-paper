library(data.table)
data.dir.list <- list(
  "jss.bench.models.rules"=function(dt){
    identical(names(dt), c("rule", "fn", "fp", "errors", "peaks")) &&
      identical(dt$rule, c("use", "remove", "join", "join2"))
  },
  "jss.bench.models.one"=function(dt){
    nrow(dt)==1 && identical(
      names(dt),
      c("prob.dir", "bedGraph.lines", "penalty", "megabytes", "seconds", 
        "peaks", "bases", "mean.pen.cost", "total.loss", "mean.intervals", 
        "max.intervals", "fn", "fp", "errors"))
  })

for(data.dir in names(data.dir.list)){
  check.fun <- data.dir.list[[data.dir]]
  base.vec <- dir(data.dir)
  id.vec <- sub(".csv$", "", base.vec)
  tall.dt <- data.table(id=id.vec)[, {
    file.csv <- file.path(data.dir, paste0(id, ".csv"))
    cat(sprintf("%5d / %5d %s\n", .I, length(id.vec), file.csv))
    dt <- fread(
      file.csv,
      colClasses=list(
        character=c("prob.dir"),
        integer=c(
          "bedGraph.lines", "peaks", "bases",
          "max.intervals", "fn", "fp", "errors"),
        numeric=c(
          "penalty", "megabytes", "seconds", 
          "mean.pen.cost", "total.loss", "mean.intervals")))
    if(check.fun(dt)){
      dt
    }else{
      str(dt);print(file.csv)
      system(paste("Rscript jss.bench.models.one.R", id))
      fread(file.csv)
    }
  }, by=list(id)]
  print(tall.dt)
  fwrite(tall.dt, paste0(data.dir, ".csv"))
  cmd <- paste0(
    "tar czvf ",
    data.dir,
    ".tgz ",
    data.dir,
    " && rm -rf ",
    data.dir)
  print(cmd)
  system(cmd)
}

jss.bench.models <- fread("jss.bench.models.csv")[order(bedGraph.lines, penalty)]
jss.bench.models.rules <- fread("jss.bench.models.rules.csv")
jss.bench.models.one <- fread("jss.bench.models.one.csv")
nrow(jss.bench.models.one)
nrow(jss.bench.models)

wide.rules <- dcast(
  jss.bench.models.rules,
  id ~ rule,
  value.var=c("peaks", "errors", "fp", "fn"))
wide.rules[peaks_join != peaks_join2]
wide.rules[errors_join != errors_join2]
wide.rules[, table(errors_join-errors_use)]
wide.rules[, mean(errors_join-errors_use)]
wide.rules[, table(better=ifelse(errors_join<errors_use, "join", ifelse(errors_use<errors_join, "use", "equal")))]
wide.rules[, table(better=ifelse(errors_remove<errors_use, "remove", ifelse(errors_use<errors_remove, "use", "equal")))]
wide.rules[, table(better=ifelse(errors_remove<errors_join, "remove", ifelse(errors_join<errors_remove, "join", "equal")))]
## errors: remove < join = join2 < use
join.dt <- jss.bench.models.rules[jss.bench.models.one, on=list(id)]
min.err <- join.dt[, list(min.errors=min(errors)), by=list(prob.dir, rule)]
min.err[, list(
  best=paste(rule[min.errors==min(min.errors)], collapse=",")
  ), by=list(prob.dir)][, sort(table(best))]
best.err <- min.err[, list(
  best.err=min(min.errors)
  ), by=list(prob.dir)][min.err, on=list(prob.dir)]
best.counts <- best.err[, list(
  n.best=sum(best.err==min.errors),
  total=.N
  ), by=list(rule)][order(n.best)]
best.counts[, percent := n.best/total]
