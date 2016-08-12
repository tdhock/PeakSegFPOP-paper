source("packages.R")
source("pick.best.index.R")

load("unsupervised.RData")
load("unsupervised.pdpa.RData")
load("dp.peaks.matrices.RData")
load("dp.peaks.sets.RData")

default.params <-
  c(macs.trained="1.30103",
    hmcan.broad.trained="2.30258509299405")

seg.mat.list <- list(
  coseg=oracle.pdpa,
  PeakSegDP=oracle.segments)

elist <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  chunk.list <- dp.peaks.matrices[[set.name]]
  for(set.i in seq_along(train.sets)){
    testSet <- paste(set.name, "split", set.i)
    test.chunks <- train.sets[[set.i]]
    train.chunks <- names(chunk.list)[! names(chunk.list) %in% test.chunks]

    baseline.list <- list()
    for(algorithm in c("hmcan.broad.trained", "macs.trained")){
      train.mat.list <- list()
      for(train.chunk in train.chunks){
        train.mat.list[[train.chunk]] <-
          chunk.list[[train.chunk]][[algorithm]]
      }
      train.mat <- do.call(rbind, train.mat.list)
      error.curve <- colSums(train.mat)
      error.sorted <- error.curve[order(as.numeric(names(error.curve)))]
      picked <- pick.best.index(error.sorted)
      plot(error.sorted)
      abline(v=picked)
      baseline.list[[algorithm]] <- names(error.sorted)[[picked]]
    }
    
    best.list <- c()
    for(algorithm in names(seg.mat.list)){
      train.list <- list()
      for(train.chunk in train.chunks){
        seg.mat <- seg.mat.list[[algorithm]][[train.chunk]]
        err.mat <- chunk.list[[train.chunk]][[algorithm]]
        err.big <- apply(seg.mat, 2, function(segs){
          err.mat[cbind(seq_along(segs), (segs-1)/2+1)]
        })
        train.list[[train.chunk]] <-
          colSums(err.big)
      }
      train.mat <- do.call(rbind, train.list)
      train.err <- colSums(train.mat)
      plot(train.err)
      picked <- pick.best.index(train.err)
      abline(v=picked)
      best.list[[algorithm]] <- picked
    }
    
    ##print(c(best.list, baseline.list))
    cat(sprintf("%d / %d %s\n", set.i, length(train.sets), set.name))
    for(test.chunk in test.chunks){
      test.info <- chunk.list[[test.chunk]]
      pred.seg.list <- list(unsupervised=list(
        PeakSegDP=unsupervised[[test.chunk]][, "oracle"],
        coseg=unsupervised.pdpa[[test.chunk]][, "oracle"]))
      for(algorithm in names(best.list)){
        best.beta <- best.list[[algorithm]]
        pred.seg.list[["supervised"]][[algorithm]] <-
          seg.mat.list[[algorithm]][[test.chunk]][, best.beta]
      }
      for(train.type in names(pred.seg.list)){
        pred.seg.by.algo <- pred.seg.list[[train.type]]
        for(algorithm in names(pred.seg.by.algo)){
          pred.seg.vec <- pred.seg.by.algo[[algorithm]]
          err.mat <- chunk.list[[test.chunk]][[algorithm]]
          pred.peaks.vec <- (pred.seg.vec-1)/2
          param.name <- as.character(pred.peaks.vec)
          sample.id <- names(pred.seg.vec)
          i.mat <- cbind(sample.id, param.name)
          errors <- err.mat[i.mat]
          stopifnot(!is.na(errors))
          elist[[paste(set.name, set.i, test.chunk, algorithm, train.type)]] <- 
            data.table(set.name, set.i, testSet, test.chunk,
                       algorithm, train.type,
                       param.name, sample.id, errors,
                       regions=test.info$regions)
        }
      }
      for(algorithm in names(default.params)){
        param.list <- list(
          unsupervised=default.params[[algorithm]],
          supervised=baseline.list[[algorithm]])
        err.mat <- test.info[[algorithm]]
        for(train.type in names(param.list)){
          param.name <- param.list[[train.type]]
          errors <- err.mat[, param.name]
          sample.id <- names(errors)
          elist[[paste(set.name, set.i, test.chunk, algorithm, train.type)]] <- 
            data.table(set.name, set.i, testSet, test.chunk,
                       algorithm=ifelse(
                         algorithm=="macs.trained", "MACS", "HMCanBroad"),
                       train.type,
                       param.name, sample.id, errors,
                       regions=test.info$regions)
        }#train.type
      }#algorithm
    }#test.chunk
  }#set.i
}#set

test.error <- do.call(rbind, elist)

save(test.error, file="test.error.RData")
