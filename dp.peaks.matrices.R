source("packages.R")

load("dp.peaks.error.RData")
load("PDPA.peaks.error.RData")
load("Segmentor.peaks.error.RData")
setkey(PDPA.peaks.error, chunk.name)
setkey(Segmentor.peaks.error, chunk.name)

set.names <- unique(sub("/.*", "", names(dp.peaks.error)))

dp.peaks.matrices <- list()
dp.peaks.matrices.fp <- list()
dp.peaks.matrices.tp <- list()
for(set.name in set.names){
  chunk.names <- grep(set.name, names(dp.peaks.error), value=TRUE)
  for(chunk.name in chunk.names){
    chunk.list <- dp.peaks.error[[chunk.name]]
    chunk.df <- do.call(rbind, chunk.list)
    if(is.data.frame(chunk.df)){
      print(chunk.name)
      chunk.dt <- data.table(chunk.df)
      chunk.dt[, param.num := as.numeric(as.character(param.name))]
      long <- chunk.dt[, list(
        errors=sum(fp+fn),
        fp=sum(fp),
        tp=sum(tp),
        possible.fp=sum(possible.fp),
        possible.tp=sum(possible.tp),
        regions=.N), by=.(sample.id, param.num)]
      pdpa <- PDPA.peaks.error[chunk.name, list(
        errors=sum(fp+fn),
        fp=sum(fp),
        tp=sum(tp),
        possible.fp=sum(possible.fp),
        possible.tp=sum(possible.tp),
        regions=.N), by=.(sample.id, peaks)]
      Seg <- Segmentor.peaks.error[chunk.name, list(
        errors=sum(fp+fn),
        fp=sum(fp),
        tp=sum(tp),
        possible.fp=sum(possible.fp),
        possible.tp=sum(possible.tp),
        regions=.N), by=.(sample.id, peaks)]
      long.list <- split(long, long$sample.id, drop=TRUE)
      err.mat <- fp.mat <- tp.mat <-
        Seg.mat <- Seg.fp <- Seg.tp <- 
          pdpa.mat <- pdpa.fp <- pdpa.tp <- 
        matrix(NA, length(long.list), 10,
               dimnames=list(sample.id=names(long.list),
                 param.name=0:9))
      for(row.i in seq_along(long.list)){
        sample.df <- long.list[[row.i]]
        param.name <- as.character(sample.df$param.num)
        err.mat[row.i, param.name] <- sample.df$errors
        fp.mat[row.i, param.name] <- sample.df$fp
        tp.mat[row.i, param.name] <- sample.df$tp
      }
      pdpa.by.peaks <- split(pdpa, pdpa$peaks)
      for(peaks.str in names(pdpa.by.peaks)){
        peaks.dt <- pdpa.by.peaks[[peaks.str]]
        pdpa.mat[paste(peaks.dt$sample.id), peaks.str] <- peaks.dt$errors
        pdpa.fp[paste(peaks.dt$sample.id), peaks.str] <- peaks.dt$fp
        pdpa.tp[paste(peaks.dt$sample.id), peaks.str] <- peaks.dt$tp
      }
      Seg.by.peaks <- split(Seg, Seg$peaks)
      for(peaks.str in names(Seg.by.peaks)){
        peaks.dt <- Seg.by.peaks[[peaks.str]] 
        Seg.mat[paste(peaks.dt$sample.id), peaks.str] <- peaks.dt$errors
        Seg.fp[paste(peaks.dt$sample.id), peaks.str] <- peaks.dt$fp
        Seg.tp[paste(peaks.dt$sample.id), peaks.str] <- peaks.dt$tp
      }
      err.list <-
        list(PeakSegDP=err.mat,
             coseg=pdpa.mat,
             Segmentor=Seg.mat,
             regions=sapply(long.list, function(x)x$regions[[1]]))
      fp.list <-
        list(PeakSegDP=fp.mat,
             coseg=pdpa.fp,
             Segmentor=Seg.fp
             possible.fp=sapply(long.list, function(x)x$possible.fp[[1]]))
      tp.list <-
        list(PeakSegDP=tp.mat,
             coseg=pdpa.tp,
             Segmentor=Seg.tp
             possible.tp=sapply(long.list, function(x)x$possible.tp[[1]]))
      for(algorithm in c("macs.trained", "hmcan.broad.trained")){
        load(sprintf("data/%s/error/%s.RData", chunk.name, algorithm))
        error.subset <- data.table(error)[sample.id %in% rownames(err.mat),]
        error.subset[, param.num := as.numeric(as.character(param.name))]
        a.df <- error.subset[, list(
          errors=sum(fp+fn),
          fp=sum(fp),
          tp=sum(tp),
          possible.fp=sum(possible.fp),
          possible.tp=sum(possible.tp),
          regions=.N), by=.(sample.id, param.num)]
        alist <- split(a.df, a.df$sample.id, drop=TRUE)
        get.mat <- function(var.name){
          alg.mat <- t(sapply(alist, "[[", var.name))
          colnames(alg.mat) <- alist[[1]]$param.num
          alg.mat
        }
        err.list[[algorithm]] <- get.mat("errors")
        fp.list[[algorithm]] <- get.mat("fp")
        tp.list[[algorithm]] <- get.mat("tp")
      }
      dp.peaks.matrices[[set.name]][[chunk.name]] <- err.list
      dp.peaks.matrices.fp[[set.name]][[chunk.name]] <- fp.list
      dp.peaks.matrices.tp[[set.name]][[chunk.name]] <- tp.list
    }#if(is.data.frame
  }#for(chunk.name
}#for(set.name

save(dp.peaks.matrices,
     dp.peaks.matrices.fp, dp.peaks.matrices.tp,
     file="dp.peaks.matrices.RData")
