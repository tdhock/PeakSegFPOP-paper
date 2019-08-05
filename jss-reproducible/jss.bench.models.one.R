source("jss-packages.R")
model.i <- 2

model.i <- commandArgs(trailingOnly=TRUE)

## Compute one model.
bench.models <- fread("jss.bench.models.csv")[order(bedGraph.lines, penalty)]
one.model <- bench.models[as.integer(model.i)]
prob.dir <- file.path("data", one.model$prob.dir)
coverage.bedGraph <- file.path(prob.dir, "coverage.bedGraph")
coverage.bedGraph.gz <- paste0(coverage.bedGraph, ".gz")
system(paste("gunzip", coverage.bedGraph.gz))
fit <- PeakSegDisk::PeakSegFPOP_dir(prob.dir, paste(one.model$penalty))
system(paste("gzip", coverage.bedGraph))

## Read labels.
labels.bed <- file.path(prob.dir, "labels.bed")
label.dt <- fread(labels.bed)
setnames(label.dt, c("chrom", "chromStart", "chromEnd", "annotation"))

## Compute different kinds of peaks based on those models.
diff.vec <- diff(fit$segments$mean)
fit$segments[, diff.after  := c(diff.vec, Inf)]
fit$segments[, diff.before := c(-Inf, diff.vec)]
seg.dt <- data.table(fit$segments)
seg.dt[status=="background" & (diff.after==0|diff.before==0), status := "peak"]
seg.dt[, peak.i := cumsum(status=="background")]
seg.dt[, st := ifelse(
  0<diff.after & diff.before<0,
  "background", "peak")]
seg.dt[, p.i := cumsum(st=="background")]
join.or.norows <- function(DT, by.var){
  if(nrow(DT)==0)return(DT)
  DT[, list(
    chromStart=min(chromStart),
    chromEnd=max(chromEnd)
  ), by=by.var]
}
peaks.list <- list(
  use=fit$segments[status=="peak"],
  remove=fit$segments[status=="peak" & diff.before != 0 & diff.after != 0],
  join=join.or.norows(seg.dt[status=="peak"], "peak.i"),
  join2=join.or.norows(seg.dt[st=="peak"], "p.i"))
rule.errors.list <- list()
for(rule in names(peaks.list)){
  peak.dt <- peaks.list[[rule]]
  error.df <- PeakError::PeakErrorChrom(peak.dt, label.dt)
  error.dt <- data.table(error.df)[, {
    data.table(
      fn=sum(fn),
      fp=sum(fp),
      errors=sum(fp+fn),
      peaks=nrow(peak.dt))
    }]
  rule.errors.list[[rule]] <- data.table(rule, error.dt)
}
rule.errors <- do.call(rbind, rule.errors.list)
col.names <- names(one.model)
one.computed <- data.table(
  prob.dir=one.model$prob.dir,
  fit$loss,
  rule.errors[rule=="use"])[, ..col.names]
rbind(
  data.table(data="stored", one.model),
  data.table(data="computed", one.computed))
dir.create("jss.bench.models.one", showWarnings=FALSE)
out.csv <- file.path("jss.bench.models.one", paste0(model.i, ".csv"))
cat("writing", out.csv, "\n")
fwrite(one.computed, out.csv)

dir.create("jss.bench.models.rules", showWarnings=FALSE)
out.csv <- file.path("jss.bench.models.rules", paste0(model.i, ".csv"))
fwrite(rule.errors, out.csv)
cat("wrote", out.csv, "\n")
