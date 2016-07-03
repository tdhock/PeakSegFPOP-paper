source("packages.R")

library(parallel)

prefix <- "http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db/"
u <- paste0(prefix, "RData.count.signal.txt")
count.files <- read.table(u, header=TRUE)
max.segments <- 19L
PDPA.timings.list <- list()
for(file.i in 1:nrow(count.files)){
  f <- count.files$file[[file.i]]
  chunk.id.dir <- dirname(paste(f))
  chunk.id <- basename(chunk.id.dir)
  set.dir <- dirname(chunk.id.dir)
  set.name <- basename(set.dir)
  u <- paste0(prefix, f)
  local.f <- file.path("data", f)
  print(local.f)
  if(!file.exists(local.f)){
    local.dir <- dirname(local.f)
    dir.create(local.dir, showWarnings=FALSE, recursive=TRUE)
    download.file(u, local.f)
  }
  time.f <- sub("counts", "PDPA.timing", local.f)
  model.f <- sub("counts", "PDPA.model", local.f)
  ## run DP based on if we have a model computed or not.
  if(file.exists(time.f)){
    load(time.f)
  }else{
    load(local.f) #counts
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
                  file.i, nrow(count.files),
                  sample.i, length(sample.list),
                  sample.id,
                  bases, n.data))
      bases.vec <- compressed$bases
      data.vec <- compressed$count
      seconds <- system.time({
        model.list <- PeakSegPDPA(data.vec, bases.vec, max.segments)
      })[["elapsed"]]
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
