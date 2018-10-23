source("packages.R")

files <- Sys.glob("../chip-seq-paper/chunks/H*/*/PDPA.model.RData")
pattern <-
  paste0("chunks/",
         "(?<set_name>.+?)",
         "/",
         "(?<chunk_id>[0-9]+)")
(matched <- str_match_named(files, pattern))

PDPA.peaks <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  regions.str <- paste0(set.name, "/", chunk.id)
  f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), f))
  load(f)
  count.file <- sub("PDPA.model", "counts", f)
  load(count.file)
  counts.by.sample <- split(counts, counts$sample.id)
  for(sample.id in names(PDPA.model)){
    fit <- PDPA.model[[sample.id]]
    count.df <- counts.by.sample[[sample.id]]
    is.feasible <- function(loss.vec){
      !any(diff(loss.vec) == 0, na.rm=TRUE)
    }
    seg.vec <- seq(1, 19, by=2)
    loss.df <- data.frame(
      segments=seg.vec,
      peaks=(seg.vec-1)/2,
      PoissonLoss=fit$loss.vec[seg.vec],
      feasible=apply(fit$mean.mat[seg.vec,], 1, is.feasible))
    feasible.df <- subset(loss.df, feasible)
    peaks.list <- list()
    for(n.segments in feasible.df$segments){
      break.vec <- if(n.segments==1){
        c()
      }else{
        fit$ends.mat[n.segments, 2:n.segments]
      }
      first <- c(1, break.vec+1)
      last <- c(break.vec, nrow(count.df))
      status.str <- rep(c("background", "peak"), l=n.segments)
      seg.df <- data.frame(
        mean=fit$mean.mat[n.segments, 1:n.segments],
        first,
        last,
        chromStart=count.df$chromStart[first],
        chromEnd=count.df$chromEnd[last],
        status=factor(status.str, c("background", "peak")),
        peaks=(n.segments-1)/2,
        segments=n.segments)
      n.peaks <- (n.segments-1)/2
      peaks.list[[paste(n.peaks)]] <- subset(seg.df, status=="peak")
    }
    PDPA.peaks[[regions.str]][[sample.id]] <- peaks.list
  }
}

save(PDPA.peaks, file="PDPA.peaks.RData")
