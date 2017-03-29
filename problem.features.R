source("packages.R")

files <- Sys.glob("data/H*/*/counts.RData")
dp.peaks.features <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  regions.str <- paste0(set.name, "/", chunk.id)
  count.f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), count.f))
  load(count.f)
  count.list <- split(counts, counts$sample.id)
  sample.ids <- names(count.list)
  chunk.mat <-
    matrix(NA, length(sample.ids), 14*4,
           dimnames=list(sample.id=sample.ids, feature=NULL))
  for(sample.id in sample.ids){
    bases <- with(count.df, chromEnd-chromStart)
    long <- rep(count.df$coverage, bases)
    n.bases <- sum(bases)
    n.data <- nrow(count.df)
    under.sqrt <- 1.1 + log(n.data/n.segments)
    in.square <- 1 + 4 * sqrt(under.sqrt)
    cleynen.factor <- in.square * in.square
    cleynen <- n.segments * cleynen.factor
    names(under.sqrt) <- names(in.square) <- names(cleynen) <- peaks.str
    count.df <- count.list[[sample.id]]
    feature.vec <-
      c(unweighted.quartile=quantile(count.df$coverage),
        weighted.quartile=quantile(long),
        unweighted.mean=mean(count.df$coverage),
        weighted.mean=mean(long),
        bases=n.bases,
        ## sqrt=under.sqrt,
        ## square=in.square,
        ## cleynen=cleynen,
        ## p.err,
        data=n.data)        
    log.features <-
      c(feature.vec,
        `log+1`=log(feature.vec+1),
        log=log(feature.vec),
        log.log=log(log(feature.vec)))
    chunk.mat[sample.id, ] <- log.features
  }
  colnames(chunk.mat) <- names(log.features)
  colnames(chunk.mat)[colMeans(is.finite(chunk.mat)) == 1]
  dp.peaks.features[[regions.str]] <- chunk.mat
}

save(problem.features, file="problem.features.RData")

