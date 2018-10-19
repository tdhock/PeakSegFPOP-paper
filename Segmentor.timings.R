source("packages.R")

count.files <- Sys.glob("../chip-seq-paper/chunks/H*/*/counts.RData")
max.segments <- 19L
Segmentor.timings.list <- list()
file.i <- 42
sample.i <- 24
file.i.vec <- seq_along(count.files)
for(file.i in file.i.vec){
  local.f <- count.files[[file.i]]
  chunk.id.dir <- dirname(paste(local.f))
  chunk.id <- basename(chunk.id.dir)
  set.dir <- dirname(chunk.id.dir)
  set.name <- basename(set.dir)
  print(local.f)
  time.f <- sub("counts", "Segmentor.timing", local.f)
  model.f <- sub("counts", "Segmentor.model", local.f)
  ## run DP based on if we have a model computed or not.
  if(file.exists(time.f)){
    load(time.f)
  }else{
    load(local.f) #counts
    PDPA.file <- sub("counts", "PDPA.model", local.f)
    load(PDPA.file)
    counts$bases <- with(counts, chromEnd-chromStart)
    counts$count <- counts$coverage
    sample.list <- split(counts, counts$sample.id, drop=TRUE)
    sample.ids <- names(sample.list)
    sample.results <- lapply(seq_along(sample.ids), function(sample.i){
      sample.id <- sample.ids[[sample.i]]
      compressed <- data.table(sample.list[[sample.id]])
      bases <- sum(compressed$bases)
      n.data <- nrow(compressed)
      cat(sprintf("%4d / %4d chunks %4d / %4d sample %s %d bases %d data\n",
                  file.i, length(count.files),
                  sample.i, length(sample.list),
                  sample.id,
                  bases, n.data))
      bases.vec <- compressed$bases
      data.vec <- compressed$count
      break.mat <- matrix(0, max.segments, max.segments)
      param.mat <- matrix(0, max.segments, max.segments)
      lik.mat <- rep(0, max.segments)
      seconds <- system.time({
        model.list <- .C(
          "SegmentPoisson",
          Size = as.integer(n.data),
          KMax = as.integer(max.segments),
          Data = as.integer(data.vec),
          DataComp = as.integer(bases.vec),
          Breakpoints = as.integer(break.mat),
          Parameters = as.double(param.mat),
          Likelihood = as.double(lik.mat),
          PACKAGE="Segmentor3IsBack")
      })[["elapsed"]]
      end.mat <- matrix(model.list$Breakpoints, max.segments, max.segments)
      mean.mat <- matrix(model.list$Parameters, max.segments, max.segments)
      mean1 <- sum(bases.vec*data.vec)/sum(bases.vec)
      lik1 <- PoissonLoss(data.vec, mean1, bases.vec)
      stopifnot(abs(lik1-model.list$Likelihood[1]) < 1e-8)
      constrained.lik <- PDPA.model[[sample.id]]$cost.mat[, n.data]
      should.be.positive <- constrained.lik-model.list$Likelihood
      compare.mat <- rbind(
        Segmentor=model.list$Likelihood,
        coseg=constrained.lik)
      if(any(should.be.positive < -1e-7)){
        print(should.be.positive)
        print(compare.mat)
        stop("coseg model more likely than Segmentorx model")
      }
      seg.list <- list()
      peak.list <- list()
      for(n.segments in 1:ncol(end.mat)){
        last <- end.mat[1:n.segments, n.segments]
        first <- c(1, last[-length(last)]+1)
        status.str <- rep(c("background", "peak"), l=n.segments)
        peaks.str <- paste((n.segments-1)/2)
        segs <- data.frame(
          mean=mean.mat[1:n.segments, n.segments],
          first,
          last,
          chromStart=compressed$chromStart[first],
          chromEnd=compressed$chromEnd[last],
          status=factor(status.str, c("background", "peak")),
          peaks=peaks.str,
          segments=n.segments)
        seg.list[[paste(n.segments)]] <- segs
        peakseg.feasible <- all(cumsum(sign(diff(segs$mean))) %in% c(0, 1))
        if(n.segments %% 2 && peakseg.feasible){
          peak.list[[peaks.str]] <- subset(segs, status=="peak")
        }
      }
      model.list$peaks <- peak.list
      model.list$segments <- seg.list
      list(
        model=model.list,
        timing=data.frame(
          set.name, chunk.id, sample.id,
          seconds, data=n.data, bases)
        )
    })
    names(sample.results) <- sample.ids
    Segmentor.model <- lapply(sample.results, "[[", "model")
    timing.list <- lapply(sample.results, "[[", "timing")
    Segmentor.timing <- do.call(rbind, timing.list)
    ## ggplot()+
    ##   geom_point(aes(data, seconds),
    ##              data=dp.timing, pch=1)+
    ##   scale_x_log10()+
    ##   scale_y_log10()
    save(Segmentor.model, file=model.f)
    save(Segmentor.timing, file=time.f)
  }
  Segmentor.timings.list[[file.i]] <- Segmentor.timing
}
Segmentor.timings <- do.call(rbind, Segmentor.timings.list)

save(Segmentor.timings, file="Segmentor.timings.RData")
