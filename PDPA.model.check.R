source("packages.R")

PDPA.RData.vec <- Sys.glob("data/*/*/PDPA.model.RData")

for(file.i in seq_along(PDPA.RData.vec)){
  PDPA.RData <- PDPA.RData.vec[[file.i]]
  load(PDPA.RData)
  chunk.dir <- dirname(PDPA.RData)
  for(compare.base in c("dp.model.RData", "dp.model.reverse.RData")){
    dp.RData <- file.path(chunk.dir, compare.base)
    objs <- load(dp.RData)
    for(sample.id in names(dp.model)){
      dp.sample <- dp.model[[sample.id]]
      pdpa.sample <- PDPA.model[[sample.id]]
      n.data <- ncol(pdpa.sample$cost.mat)
      loss.dt <- data.table(dp.sample$error)
      loss.dt[, dp := error]
      loss.dt[, pdpa := pdpa.sample$cost.mat[segments, n.data] ]
      loss.dt[, should.be.positive := dp - pdpa]
      bad.loss <- loss.dt[should.be.positive < -1e-9, ]
      if(nrow(bad.loss)){
        cat(sprintf("%s %s %s\n", PDPA.RData, dp.RData, sample.id))
        print(bad.loss)
        load(file.path(chunk.dir, "dp.timing.RData"))
        print(dp.timing[sample.id,])
        load(file.path(chunk.dir, "counts.RData"))
        counts.by.sample <- split(counts, counts$sample.id)
        one.sample <- counts.by.sample[[sample.id]]
        dp.segments <-
          data.table(dp.sample$segments)[segments %in% bad.loss$segments,]
        pdpa.segments.list <- list()
        for(n.segments in bad.loss$segments){
          break.vec <- pdpa.sample$ends.mat[n.segments, 2:n.segments]
          first <- c(1, break.vec+1)
          last <- c(break.vec, n.data)
          pdpa.segments.list[[paste(n.segments)]] <- data.table(
            mean=pdpa.sample$mean.mat[n.segments, 1:n.segments],
            first,
            last,
            chromStart=one.sample$chromStart[first],
            chromEnd=one.sample$chromEnd[last],
            status=rep(c("background", "peak"), l=n.segments),
            peaks=(n.segments-1)/2,
            segments=n.segments)
        }
        pdpa.segments <- do.call(rbind, pdpa.segments.list)
        rbind(
          pdpa.last=pdpa.segments$last,
          dp.last=dp.segments$last)
        data.vec <- one.sample$coverage
        weight.vec <- with(one.sample, chromEnd-chromStart)
        fit <- cDPA(data.vec, weight.vec, 19L)
        all.loss <- data.table(
          dp=as.numeric(fit$loss),
          pdpa=as.numeric(pdpa.sample$cost.mat),
          segments=as.integer(row(fit$loss)),
          data=as.integer(col(fit$loss)))
        all.loss[, should.be.positive := dp - pdpa]
        all.loss[should.be.positive < -1e-8,]
        stop("dp model more likely than pdpa model")
      }
    }
  }
}
