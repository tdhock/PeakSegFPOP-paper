load("PDPA.cDPA.compare.RData")

simplest <-
  PDPA.cDPA.compare$loss[, .SD[which.min(segments),],
                         by=.(chunk.name, sample.id)]

PDPA.cDPA.compare$loss[order(segments, -diffBases),]
PDPA.cDPA.compare$loss[order(segments, dp.fwd-PDPA),]
PDPA.cDPA.compare$loss[order(diffBases),]
PDPA.cDPA.compare$loss[order(dp.fwd-PDPA),]

nrow(PDPA.cDPA.compare$loss)

simplest[order(diffBases),]
simplest[order(dp.fwd-PDPA),]

PDPA.cDPA.compare$loss[peaks %in% c(1,2),]

chunk.summary <- PDPA.cDPA.compare$loss[, list(
  diffBases=sum(diffBases),
  diffLoss=sum(dp.fwd-PDPA)
  ), by=.(chunk.name, peaks)]
chunk.summary[order(peaks, -diffBases),]

## Total number of models.
nrow(PDPA.cDPA.compare$all.loss)

## Number of problems for which the DP does not recover all 10 models.
no.dp <- PDPA.cDPA.compare$all.loss[is.na(dp.fwd), ]
length(no.dp[, table(paste(chunk.name, sample.id))])

## Number of problems for which the DP recovers a sub-optimal model.
dp.suboptimal <- PDPA.cDPA.compare$all.loss[PDPA+1e-6 < dp.fwd,]
nrow(dp.suboptimal)

## Number of problems for which the DP recovers a sub-optimal model,
## and the PeakSeg solution exists.
nrow(dp.suboptimal[PDPA.feasible==TRUE,])
nrow(PDPA.cDPA.compare$loss)

## Number of problems for which the DP does not recover the optimal
## solution, coseg recovers the optimal solution, and Segmentor does
## too.
nrow(PDPA.cDPA.compare$loss[Seg.feasible==TRUE,])

ggplot()+
  scale_fill_brewer(palette="Blues")+
  geom_point(aes(log10(dp.fwd-PDPA), log10(diffBases), fill=factor(peaks)),
             shape=21,
             color="black",
             data=PDPA.cDPA.compare$loss)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("peaks")+
  geom_point(aes(log10(dp.fwd-PDPA), log10(diffBases),
                 fill=Seg.feasible),
             shape=21,
             color="black",
             data=PDPA.cDPA.compare$loss)
quantile(PDPA.cDPA.compare$loss$diffBases)
PDPA.cDPA.compare$loss[order(-diffBases),][1:20,]

ggplot()+
  scale_fill_brewer(palette="Blues")+
  geom_point(aes(dp.fwd-PDPA, diffBases, fill=factor(peaks)),
             shape=21,
             color="black",
             data=PDPA.cDPA.compare$loss)

to.show.list <- list(
  "H3K4me3_TDH_immune/14"=29,
  "H3K4me3_PGP_immune/30"=29,
  "H3K36me3_AM_immune/10"=c(25, 24),
  "H3K4me3_PGP_immune/22"=c(104, 59)
  )

loss.dt <- PDPA.cDPA.compare$all.loss
setkey(loss.dt, chunk.name)

for(chunk.name in names(to.show.list)){
  chunk.dir <- file.path("data", chunk.name)
  load(file.path(chunk.dir, "counts.RData"))
  load(file.path(chunk.dir, "PDPA.model.RData"))
  load(file.path(chunk.dir, "dp.model.RData"))
  counts.by.sample <- split(counts, counts$sample.id)
  for(n.peaks in 1:2){
    n.segs <- n.peaks*2+1
    chunk.loss <- loss.dt[chunk.name, ][segments==n.segs, ]
    chunk.loss[, diff := dp.fwd - PDPA]
    setkey(chunk.loss, sample.id)
    chunk.ord <- chunk.loss[order(PDPA.feasible, diff), ]
    simplify <- function(x)gsub("McGill0", "", x)
    levs <- simplify(chunk.ord$sample.id)
    counts$simple.id <- factor(simplify(counts$sample.id), levs)

    type.df <- unique(counts[, c("cell.type", "sample.id")])
    rownames(type.df) <- type.df$sample.id

    peaks.list <- list()
    for(sample.id in names(dp.model)){
      cell.type <- type.df[sample.id, "cell.type"]
      sample.dp <- dp.model[[sample.id]]
      sample.counts <- counts.by.sample[[sample.id]]
      PDPA <- PDPA.model[[sample.id]]
      mean.vec <- PDPA$mean.mat[n.segs, 1:n.segs]
      is.feasible <- all(diff(mean.vec) != 0)
      simple.id <- factor(simplify(sample.id), levs)
      if(is.feasible){
        change.vec <- PDPA$ends.mat[n.segs, 2:n.segs]
        last <- c(change.vec, nrow(sample.counts))
        first <- c(1, change.vec+1)
        peaks.list[[paste(sample.id, "coseg")]] <- data.table(
          algorithm="coseg",
          cell.type,
          sample.id,
          simple.id,
          first,
          last,
          chromStart=sample.counts$chromStart[first],
          chromEnd=sample.counts$chromEnd[last],
          peaks=n.peaks,
          segments=n.segs)[c(FALSE, TRUE),]
      }
      peaks.list[[paste(sample.id, "PeakSegDP")]] <- data.table(
        algorithm="PeakSegDP",
        cell.type,
        sample.id,
        simple.id,
        sample.dp$peaks[[paste(n.peaks)]])
    }
    peaks <- do.call(rbind, peaks.list)

    ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(simple.id ~ ., scales="free")+
      geom_step(aes(chromStart/1e3, coverage),
                color="grey50",
                data=counts)+
      ## geom_point(aes(chromStart/1e3, 0,
      ##                color=algorithm),
      ##            data=dp.peaks,
      ##            size=5,
      ##            shape=1)+
      scale_size_manual(values=c(coseg=3, PeakSegDP=1))+
      geom_segment(aes(chromStart/1e3, 0,
                       xend=chromEnd/1e3, yend=0,
                       size=algorithm,
                       color=algorithm),
                   data=peaks)

  }
}
