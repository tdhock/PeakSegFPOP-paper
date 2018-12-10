source("packages.R")
source("pick.best.index.R")

load("unsupervised.RData")
load("unsupervised.inf.RData")
load("unsupervised.Segmentor.RData")
load("dp.peaks.matrices.RData")
load("dp.peaks.sets.RData")

default.params <- c(
  macs.trained="1.30103",
  hmcan.broad.trained="2.30258509299405")
seg.mat.list <- list(
  coseg.inf=oracle.inf,
  Segmentor=oracle.Segmentor,
  PeakSegDP=oracle.segments)

elist <- list()
roc.list <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  chunk.list <- dp.peaks.matrices[[set.name]]
  tp.list <- dp.peaks.matrices.tp[[set.name]]
  fp.list <- dp.peaks.matrices.fp[[set.name]]
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
        if(any(is.na(err.mat)) && algorithm != "PeakSegDP"){
          stop("missing values in error matrix")
        }
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
      test.tp <- tp.list[[test.chunk]]
      test.fp <- fp.list[[test.chunk]]
      pred.seg.list <- list(
        unsupervised=list(
          PeakSegDP=unsupervised[[test.chunk]][, "oracle"],
          Segmentor=unsupervised.Segmentor[[test.chunk]][, "oracle"],
          coseg.inf=unsupervised.inf[[test.chunk]][, "oracle"]))
      for(algorithm in names(best.list)){
        best.beta <- best.list[[algorithm]]
        pred.seg.list[["supervised"]][[algorithm]] <-
          seg.mat.list[[algorithm]][[test.chunk]][, best.beta]
      }
      ## optimal segmentation based methods.
      for(train.type in names(pred.seg.list)){
        pred.seg.by.algo <- pred.seg.list[[train.type]]
        for(algorithm in names(pred.seg.by.algo)){
          pred.seg.vec <- pred.seg.by.algo[[algorithm]]
          err.mat <- test.info[[algorithm]]
          tp.mat <- test.tp[[algorithm]]
          fp.mat <- test.fp[[algorithm]]
          sample.id <- names(pred.seg.vec)
          param.name <- as.character(pred.seg.vec)
          i.mat <- cbind(sample.id, param.name)
          if(train.type=="supervised"){
            seg.mat <- seg.mat.list[[algorithm]][[test.chunk]]
            ## seg.mat[sample.id, param.name]
            tp.big <- apply(seg.mat, 2, function(segs){
              tp.mat[cbind(seq_along(segs), (segs-1)/2+1)]
            })
            fp.big <- apply(seg.mat, 2, function(segs){
              fp.mat[cbind(seq_along(segs), (segs-1)/2+1)]
            })
            roc.list[[paste(set.name, set.i, test.chunk, algorithm)]] <-
              data.table(set.name, set.i, test.chunk, algorithm,
                         param.i=1:ncol(tp.big),
                         tp=colSums(tp.big),
                         possible.tp=sum(test.tp$possible.tp),
                         fp=colSums(fp.big),
                         possible.fp=sum(test.fp$possible.fp))
          }
          elist[[paste(set.name, set.i, test.chunk, algorithm, train.type)]] <- 
            data.table(set.name, set.i, testSet, test.chunk,
                       algorithm, train.type,
                       ## param.name is the number of peaks here.
                       param.name, sample.id,
                       errors=err.mat[i.mat],
                       tp=tp.mat[i.mat],
                       possible.tp=test.tp$possible.tp,
                       fp=fp.mat[i.mat],
                       possible.fp=test.fp$possible.fp,
                       regions=test.info$regions)
        }#for(algorithm
      }#for(train.type
      ## other baseline methods:
      for(algorithm in names(default.params)){
        param.list <- list(
          unsupervised=default.params[[algorithm]],
          supervised=baseline.list[[algorithm]])
        err.mat <- test.info[[algorithm]]
        tp.mat <- test.tp[[algorithm]]
        fp.mat <- test.fp[[algorithm]]
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
                       tp=tp.mat[, param.name],
                       possible.tp=test.tp$possible.tp,
                       fp=fp.mat[, param.name],
                       possible.fp=test.fp$possible.fp,
                       regions=test.info$regions)
        }
        ord <- order(as.numeric(colnames(err.mat)))
        roc.dt <- data.table(
          set.name, set.i, test.chunk, 
          algorithm=ifelse(
            algorithm=="macs.trained", "MACS", "HMCanBroad"),
          param.i=NA_integer_,
          tp=colSums(tp.mat),
          possible.tp=sum(test.tp$possible.tp),
          fp=colSums(fp.mat),
          possible.fp=sum(test.fp$possible.fp))[ord,]
        roc.dt[, param.i := 1:.N]
        roc.list[[paste(set.name, set.i, test.chunk, algorithm)]] <-
          roc.dt
      }#algorithm
    }#test.chunk
  }#set.i
}#set

test.error <- do.call(rbind, elist)
roc <- do.call(rbind, roc.list)

save(test.error, roc, file="test.error.RData")
