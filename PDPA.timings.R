source("packages.R")

library(data.table)
library(PeakSegOptimal)

count.files <- Sys.glob("../chip-seq-paper/chunks/H*/*/counts.RData")
##max.segments <- 199L
max.segments <- 19L
PDPA.timings.list <- list()
file.i <- 63
sample.i <- 20
file.i.vec <- seq_along(count.files)
for(file.i in file.i.vec){
  local.f <- count.files[[file.i]]
  chunk.id.dir <- dirname(paste(local.f))
  chunk.id <- basename(chunk.id.dir)
  set.dir <- dirname(chunk.id.dir)
  set.name <- basename(set.dir)
  print(local.f)
  time.f <- sub("counts", "PDPA.timing", local.f)
  model.f <- sub("counts", "PDPA.model", local.f)
  ## run DP based on if we have a model computed or not.
  if(file.exists(time.f)){
    load(time.f)
  }else{
    load(local.f) #counts
    if(FALSE){
      load(sub("counts", "regions", local.f))
      regions.dt <- data.table(regions)
    }
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
      seconds <- system.time({
        model.list <- PeakSegPDPA(data.vec, bases.vec, max.segments)
      })[["elapsed"]]
      if(FALSE){
        fit <- model.list
        n.segments <- max.segments
        mean.vec <- fit$mean.mat[n.segments, 1:n.segments]
        diff.vec <- diff(mean.vec)
        break.vec <- if(n.segments==1){
          c()
        }else{
          fit$ends.mat[n.segments, 2:n.segments]
        }
        first <- c(1, break.vec+1)
        last <- c(break.vec, nrow(compressed))
        status.str <- rep(c("background", "peak"), l=n.segments)
        peak.dt <- data.table(
          mean=mean.vec,
          first,
          last,
          is.peak=(seq_along(mean.vec)-1) %% 2,
          diff.before=c(Inf, diff.vec),
          diff.after=c(diff.vec, Inf))
        peak.dt[is.peak==0 & (diff.after==0|diff.before==0), is.peak := 1]
        peak.dt[, peak.i := cumsum(is.peak==0)]
        only.peaks <- peak.dt[is.peak==1, list(
          chromStart=compressed$chromStart[min(first)],
          chromEnd=compressed$chromEnd[max(last)]),
          by=list(peak.i)]
        ggplot()+
          geom_segment(aes(
            chromStart, peak.i,
            xend=chromEnd, yend=peak.i),
            data=only.peaks)
        sample.regions <- regions.dt[sample.id, on=list(sample.id)]
        PeakError::PeakErrorChrom(only.peaks, sample.regions)
      }
      dp.file <- sub("counts", "dp.model", local.f)
      if(file.exists(dp.file)){
        load(dp.file)
        dp <- dp.model[[sample.id]]$error
        if(is.numeric(dp$error)){
          compare.mat <- rbind(
            dp=dp$error,
            pdpa=model.list$cost.mat[dp$segments, n.data])
          should.be.positive <- compare.mat["dp",]-compare.mat["pdpa",]
          if(any(should.be.positive < -1e-6)){
            print(should.be.positive)
            print(compare.mat)
            stop("dp model more likely than pdpa model")
            fit <- cDPA(data.vec, bases.vec, 19L)
            all.loss <- data.table(
              dp=as.numeric(fit$loss),
              pdpa=as.numeric(model.list$cost.mat),
              segments=as.integer(row(fit$loss)),
              data=as.integer(col(fit$loss)))
            all.loss[, should.be.positive := dp - pdpa]
            all.loss[should.be.positive < -1e-8,]
          }
        }
      }
      list(
        model=model.list,
        timing=data.frame(
          set.name, chunk.id, sample.id,
          seconds, data=n.data, bases)
        )
    })
    names(sample.results) <- sample.ids
    PDPA.model <- lapply(sample.results, "[[", "model")
    timing.list <- lapply(sample.results, "[[", "timing")
    PDPA.timing <- do.call(rbind, timing.list)
    ## ggplot()+
    ##   geom_point(aes(data, seconds),
    ##              data=dp.timing, pch=1)+
    ##   scale_x_log10()+
    ##   scale_y_log10()
    save(PDPA.model, file=model.f)
    save(PDPA.timing, file=time.f)
  }
  PDPA.timings.list[[file.i]] <- PDPA.timing
}
PDPA.timings <- do.call(rbind, PDPA.timings.list)

save(PDPA.timings, file="PDPA.timings.RData")
