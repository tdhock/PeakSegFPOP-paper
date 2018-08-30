problem.findPeaks <- function
### Compute a model with better (more likely) peaks, with at most the
### number of peaks given by peaks.int. This function repeated calls
### problem.PeakSegFPOP with different penalty values, until either
### (1) it finds the peaks.int model, or (2) it concludes that there
### is no peaks.int model, in which case it returns the next simplest
### model (with fewer peaks than peaks.int).
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  getCandidateIndex,
### function: how to choose next/optimal model?
  verbose=0
### Print messages?
 ){
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  stopifnot(is.function(getCandidateIndex))
  PeakSegPipeline::problem.coverage(problem.dir)
  model.list <- list()
  next.pen <- c(0, Inf)
  iteration <- 0
  while(length(next.pen)){
    if(verbose)cat(
      "Next =", paste(next.pen, collapse=", "),
      "mc.cores=", getOption("mc.cores"),
      "\n")
    next.str <- paste(next.pen)
    iteration <- iteration+1
    model.list[next.str] <- PeakSegPipeline::mclapply.or.stop(
      next.str, function(penalty.str){
        L <- PeakSegPipeline::problem.PeakSegFPOP(problem.dir, penalty.str)
        L$loss$seconds <- L$timing$seconds
        L$loss$megabytes <- L$timing$megabytes
        L
      }
    )
    loss.list <- lapply(model.list, "[[", "loss")
    loss.dt <- do.call(rbind, loss.list)[order(penalty)]
    if(2 <= verbose){
      print(loss.dt[, .(penalty, segments, peaks, bases, mean.pen.cost, total.cost, status)])
    }
    if(!is.numeric(loss.dt$penalty)){
      stop("penalty column is not numeric -- check loss in _loss.tsv files")
    }
    unique.peaks <- loss.dt[, data.table(
      .SD[which.min(iteration)],
      penalties=.N
    ), by=list(peaks)]
    path.dt <- data.table(penaltyLearning::modelSelection(
      unique.peaks, "total.cost", "peaks"))
    path.dt[, next.pen := max.lambda]
    path.dt[, already.computed := next.pen %in% names(loss.list)]
    path.dt[, no.next := c(diff(peaks) == -1, NA)]
    path.dt[, done := already.computed | no.next]
    candidate.i <- getCandidateIndex(path.dt)
    if(!(
      is.integer(candidate.i) &&
      length(candidate.i)==1 &&
      1 <= candidate.i && candidate.i <= nrow(path.dt)-1
    )){
      str(candidate.i)
      stop(
        "getCandidateIndex must return integer scalar ",
        "between 1 and N-1 ",
        "where N is the number of selectable models")
    }
    candidate <- path.dt[candidate.i]
    next.pen <- if(!candidate$done)candidate$next.pen
  }#while(!is.null(pen))
  out <- model.list[[paste(candidate$penalty)]]
  out$others <- path.dt
  out
### TODO
}

problem.betterPeaks <- function
### TODO
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  peaks.int,
### integer scalar: number of peaks to find.
  verbose=0
### Print messages?
){
  stopifnot(
    is.integer(peaks.int) &&
    length(peaks.int)==1 &&
    0 <= peaks.int)
  problem.findPeaks(problem.dir, function(dt){
    dt[, max(which(peaks.int <= peaks))]
  }, verbose)
}

problem.fewestBetterPeaks <- function
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  total.cost.upper.bound,
### numeric scalar: upper bound on loss (we want to find the model
### with fewest peaks, and loss which is below this quantity).
  verbose=0
### Print messages?
){
  stopifnot(
    is.numeric(total.cost.upper.bound) &&
    length(total.cost.upper.bound)==1)
  problem.findPeaks(problem.dir, function(dt){
    dt[, max(which(total.cost <= total.cost.upper.bound))]
  }, verbose)
}

problem.mostFeasibleBetterPeaks <- function
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  total.cost.upper.bound,
### numeric scalar: upper bound on loss (we want to find the model
### with fewest peaks, and loss which is below this quantity).
  verbose=0
### Print messages?
){
  stopifnot(
    is.numeric(total.cost.upper.bound) &&
    length(total.cost.upper.bound)==1)
  problem.findPeaks(problem.dir, function(dt){
    i <- dt[, min(which(total.cost <= total.cost.upper.bound & status=="feasible"))]
    if(i==dt[, .N]){
      i-1L
    }else{
      i
    }
  }, verbose)
}

