source("packages.R")

load("PDPA.infeasible.error.RData")
load("dp.peaks.error.RData")
load("Segmentor.infeasible.error.RData")

set.names <- grep("^H", unique(sub("/.*", "", names(dp.peaks.error))), value=TRUE)
dp.peaks.matrices <- list()
dp.peaks.matrices.fp <- list()
dp.peaks.matrices.tp <- list()
set.i <- 1
chunk.name <- "H3K36me3_AM_immune/13"
for(set.i in seq_along(set.names)){
  set.name <- set.names[[set.i]]
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
        regions=.N), by=.(sample.id, segments=param.num*2+1)]
      PDPA.select.dt <- data.table(chunk.name, rule="remove")
      pdpa.inf <- PDPA.infeasible.error[PDPA.select.dt, list(
        errors=sum(fp+fn),
        fp=sum(fp),
        tp=sum(tp),
        possible.fp=sum(possible.fp),
        possible.tp=sum(possible.tp),
        regions=.N), by=.(
          sample.id, segments
        ), on=list(
          chunk.name, rule)]
      Seg.select.dt <- data.table(chunk.name, rule="rm")
      Seg <- Segmentor.infeasible.error[Seg.select.dt, list(
        errors=sum(fp+fn),
        fp=sum(fp),
        tp=sum(tp),
        possible.fp=sum(possible.fp),
        possible.tp=sum(possible.tp),
        regions=.N), by=.(
          sample.id, segments
        ), on=list(
          rule,
          chunk.name)]
      long.list <- split(long, long$sample.id, drop=TRUE)
      err.mat <- fp.mat <- tp.mat <-
        Seg.mat <- Seg.fp <- Seg.tp <- 
          pdpa.mat <- pdpa.fp <- pdpa.tp <- 
            inf.mat <- inf.fp <- inf.tp <- 
        matrix(NA, length(long.list), 10,
               dimnames=list(sample.id=names(long.list),
                 segments=seq(1, 19, by=2)))
      for(row.i in seq_along(long.list)){
        sample.df <- long.list[[row.i]]
        param.name <- as.character(sample.df$segments)
        stopifnot(0<length(param.name))
        err.mat[row.i, param.name] <- sample.df$errors
        fp.mat[row.i, param.name] <- sample.df$fp
        tp.mat[row.i, param.name] <- sample.df$tp
      }
      inf.by.segments <- split(pdpa.inf, pdpa.inf$segments)
      for(segments.str in names(inf.by.segments)){
        peaks.dt <- inf.by.segments[[segments.str]]
        inf.mat[paste(peaks.dt$sample.id), segments.str] <- peaks.dt$errors
        inf.fp[paste(peaks.dt$sample.id), segments.str] <- peaks.dt$fp
        inf.tp[paste(peaks.dt$sample.id), segments.str] <- peaks.dt$tp
      }
      Seg.by.segments <- split(Seg, Seg$segments)
      for(segments.str in names(Seg.by.segments)){
        peaks.dt <- Seg.by.segments[[segments.str]] 
        Seg.mat[paste(peaks.dt$sample.id), segments.str] <- peaks.dt$errors
        Seg.fp[paste(peaks.dt$sample.id), segments.str] <- peaks.dt$fp
        Seg.tp[paste(peaks.dt$sample.id), segments.str] <- peaks.dt$tp
      }
      err.list <-
        list(PeakSegDP=err.mat,
             coseg.inf=inf.mat,
             Segmentor=Seg.mat,
             regions=sapply(long.list, function(x)x$regions[[1]]))
      fp.list <-
        list(PeakSegDP=fp.mat,
             coseg.inf=inf.fp,
             Segmentor=Seg.fp,
             possible.fp=sapply(long.list, function(x)x$possible.fp[[1]]))
      tp.list <-
        list(PeakSegDP=tp.mat,
             coseg.inf=inf.tp,
             Segmentor=Seg.tp,
             possible.tp=sapply(long.list, function(x)x$possible.tp[[1]]))
      for(algorithm in c("macs.trained", "hmcan.broad.trained")){
        load(sprintf("../chip-seq-paper/chunks/%s/error/%s.RData", chunk.name, algorithm))
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
