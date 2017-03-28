source("packages.R")

library(penaltyLearning)

selection.list <- list()
model.files <- Sys.glob("data/H*/*/PDPA.model.RData")
i.vec <- seq_along(model.files)
##i.vec <- 1:2
for(model.file.i in i.vec){
  model.file <- model.files[[model.file.i]]
  chunk.name <- sub("data/", "", dirname(model.file))
  cat(sprintf("%4d / %4d %s\n", model.file.i, length(model.files), chunk.name))
  load(model.file)
  segments.list <- list()
  oseg.list <- list()
  for(sample.id in names(PDPA.model)){
    fit <- PDPA.model[[sample.id]]
    lik.vec <- rep(NA, 19L)
    force.na <- rep(FALSE, 19L)
    for(n.segments in seq_along(lik.vec)){
      seg.mean.vec <- fit$mean.mat[n.segments, 1:n.segments]
      ends.vec <- c(fit$ends.mat[n.segments, 1:n.segments], fit$n.data)
      is.feasible <- all(diff(seg.mean.vec) != 0)
      lik.vec[[n.segments]] <- if(is.feasible && n.segments %% 2){
        data.mean.vec <- rep(seg.mean.vec, diff(ends.vec))
        data.lik.vec <- dpois(fit$count.vec, data.mean.vec, log=TRUE)
        last.lik <- -sum(data.lik.vec)
      }else{
        force.na[[n.segments]] <- TRUE
        last.lik
      }
    }
    n.bases <- sum(fit$weight.vec)
    sample.oracle <- alice.oracle(n.bases, lik.vec)
    oseg.list[[sample.id]] <- my.oracle(n.bases, lik.vec)
    segSeq <- seq(1, 19, by=2)
    write.na <- force.na[segSeq]
    penalty.mat <-  #peaks x penalties
      cbind(oracle=sample.oracle$crit[segSeq],
            none=lik.vec[segSeq])
    rownames(penalty.mat) <- segSeq
    penalty.mat[write.na, ] <- NA
    selected.rows <- apply(penalty.mat, 2, which.min)
    selected.segs <- segSeq[selected.rows]
    names(selected.segs) <- names(selected.rows)
    segments.list[[sample.id]] <- selected.segs
  }
  unsupervised.pdpa[[chunk.name]] <- do.call(rbind, segments.list)
  oracle.pdpa[[chunk.name]] <- do.call(rbind, oseg.list)
}
test <- do.call(rbind, unsupervised.pdpa)
stopifnot(!is.null(test))


load("PDPA.peaks.error.RData")
load("PDPA.peaks.RData")

error.dt <- data.table(PDPA.peaks.error)[, list(
  possible.fn=sum(possible.tp),
  possible.fp=sum(possible.fp),
  fp=sum(fp),
  fn=sum(fn),
  errors=sum(fn+fp),
  labels=.N
  ), by=list(chunk.name, sample.id, peaks)]
