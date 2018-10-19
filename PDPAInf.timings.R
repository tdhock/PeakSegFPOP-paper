source("packages.R")

count.files <- Sys.glob("../chip-seq-paper/chunks/H*/*/counts.RData")
max.segments <- 19L
PDPAInf.timings.list <- list()
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
  time.f <- sub("counts", "PDPAInf.timing", local.f)
  model.f <- sub("counts", "PDPAInf.model", local.f)
  ## run DP based on if we have a model computed or not.
  if(file.exists(time.f)){
    load(time.f)
  }else{
    load(local.f) #counts
    counts$bases <- with(counts, chromEnd-chromStart)
    counts$count <- counts$coverage
    sample.list <- split(counts, counts$sample.id, drop=TRUE)
    sample.ids <- names(sample.list)
    PDPA.file <- sub("counts", "PDPA.model", local.f)
    objs <- load(PDPA.file)
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
        model.list <- PeakSegPDPAInf(data.vec, bases.vec, max.segments)
      })[["elapsed"]]
      PDPA.loss.vec <- PDPA.model[[sample.id]]$loss.vec
      PDPAInf.loss.vec <- with(model.list, cost.mat[, n.data])
      stopifnot(all.equal(PDPA.loss.vec, PDPAInf.loss.vec))
      small.list <- with(model.list, {
        intervals.vec <- intervals.mat[0 < intervals.mat]
        list(
          loss.vec=PDPA.loss.vec,
          mean.mat=mean.mat,
          ends.mat=ends.mat,
          mean.intervals=mean(intervals.vec),
          max.intervals=max(intervals.vec))
      })
      list(
        model=small.list,
        timing=data.frame(
          set.name, chunk.id, sample.id,
          seconds, data=n.data, bases)
        )
    })
    names(sample.results) <- sample.ids
    PDPAInf.model <- lapply(sample.results, "[[", "model")
    timing.list <- lapply(sample.results, "[[", "timing")
    PDPAInf.timing <- do.call(rbind, timing.list)
    ## ggplot()+
    ##   geom_point(aes(data, seconds),
    ##              data=dp.timing, pch=1)+
    ##   scale_x_log10()+
    ##   scale_y_log10()
    save(PDPAInf.model, file=model.f)
    save(PDPAInf.timing, file=time.f)
  }
  PDPAInf.timings.list[[file.i]] <- PDPAInf.timing
}
PDPAInf.timings <- do.call(rbind, PDPAInf.timings.list)

save(PDPAInf.timings, file="PDPAInf.timings.RData")
