load("dp.peaks.matrices.RData")

sapply(dp.peaks.matrices, length)

## How many folds in cross-validation?
n.folds <- 4

set.seed(1)

dp.peaks.sets <- list()
for(set.name in names(dp.peaks.matrices)){
  set <- dp.peaks.matrices[[set.name]]
  chunk.names <- names(set)
  fold.vec <- sample(rep(1:n.folds, l=length(set)))
  test.sets <- list()
  for(fold.i in 1:n.folds){
    test.sets[[fold.i]] <- chunk.names[fold.vec == fold.i]
  }
  dp.peaks.sets[[set.name]] <- test.sets
}

save(dp.peaks.sets, file="dp.peaks.sets.RData")
