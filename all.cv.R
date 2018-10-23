source("packages.R")

load("all.features.RData")
load("dp.peaks.sets.RData")
load("all.modelSelection.RData")
all.modelSelection[, row.id := paste(chunk.name, sample.id)]
all.modelSelection[, errors := total.errors]
all.modelSelection[, fp := total.fp]
all.modelSelection[, fn := total.fn]

load("all.targets.RData")
all.targets[, row.id := paste(chunk.name, sample.id)]
setkey(all.targets, algo, row.id)
algo.vec <- unique(all.targets$algo)

##future::plan(multi

roc.list <- list()
roc.thresh.list <- list()
for(set.name in names(dp.peaks.sets)){
  fold.list <- dp.peaks.sets[[set.name]]
  for(fold.i in seq_along(fold.list)){
    testSet <- paste(set.name, "split", fold.i)
    train.chunk.list <- fold.list[-fold.i]
    chunk.fold.list <- list()
    for(chunk.fold in seq_along(train.chunk.list)){
      some.train.chunks <- train.chunk.list[[chunk.fold]]
      chunk.fold.list[[chunk.fold]] <- structure(
        rep(chunk.fold, length(some.train.chunks)),
        names=some.train.chunks)
    }
    chunk.fold.vec <- unlist(chunk.fold.list)
    set.chunk.list <- list(
      test=fold.list[[fold.i]],
      train=unlist(train.chunk.list))
    set.features.list <- list()
    for(train.or.test in names(set.chunk.list)){
      set.chunk.vec <- set.chunk.list[[train.or.test]]
      L <- list()
      for(chunk.name in set.chunk.vec){
        f.mat <- all.features[[chunk.name]]
        rownames(f.mat) <- paste(chunk.name, rownames(f.mat))
        L[[chunk.name]] <- f.mat
      }
      set.features.list[[train.or.test]] <- do.call(rbind, L)
    }
    for(algo in algo.vec){
      select.dt <- data.table(algo, row.id=rownames(set.features.list$train))
      train.target.dt <- all.targets[select.dt]
      train.fold.vec <- chunk.fold.vec[train.target.dt$chunk.name]
      train.target.mat <- train.target.dt[, cbind(
        min.log.lambda, max.log.lambda)]
      any.finite <- apply(is.finite(train.target.mat), 1, any)
      make.pred.mat <- function(x){
        matrix(x, nrow(set.features.list$test), dimnames=list(
          rownames(set.features.list$test)))
      }
      pred.mat <- if(all(train.target.mat[,1] == -Inf)){
        make.pred.mat(min(train.target.mat[,2])-0.1)
      }else if(all(train.target.mat[,2] == Inf)){
        make.pred.mat(max(train.target.mat[,1])+0.1)
      }else{
        fit <- IntervalRegressionCV(
          set.features.list$train[any.finite,],
          train.target.mat[any.finite,],
          fold.vec=train.fold.vec[any.finite],
          verbose=0)
        fit$predict(set.features.list$test)
      }
      pred.dt <- data.table(
        algo,
        row.id=rownames(pred.mat))
      models.test <- all.modelSelection[pred.dt, on=list(
        algo, row.id)]
      pred.dt[, pred.log.lambda := as.numeric(pred.mat)]
      L <- ROChange(models.test, pred.dt, "row.id")
      roc.list[[paste(set.name, fold.i, algo)]] <- data.table(
        set.name, fold.i, algo, L$roc)
      roc.thresh.list[[paste(set.name, fold.i, algo)]] <- with(L, data.table(
        set.name, fold.i, algo, thresholds, auc))
    }
  }
}
roc <- do.call(rbind, roc.list)
roc.thresh <- do.call(rbind, roc.thresh.list)

roc.thresh[, accuracy.percent := 100-error.percent]
dcast(roc.thresh[threshold=="predicted"], set.name+fold.i~algo, value.var=c("auc", "accuracy.percent"))

roc.tall <- melt(
  roc.thresh[threshold=="predicted"],
  measure.vars=c("auc", "accuracy.percent"),
  id.vars=c("set.name", "fold.i", "algo"))
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(variable ~ set.name, scales="free")+
  geom_point(aes(
    algo, value),
    data=roc.tall)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(set.name ~ variable, scales="free")+
  geom_point(aes(
    value, algo),
    data=roc.tall)

save(roc, roc.thresh, file="all.cv.RData")
