source("packages.R")
chunk.name <- "H3K4me3_TDH_immune/11";sample.id <- "McGill0009";peaks <- 3
chunk.name <- "H3K4me3_PGP_immune/16";sample.id <- "McGill0005";peaks <- 3
chunk.name <- "H3K4me3_TDH_immune/11";sample.id <- "McGill0008";peaks <- 3
chunk.name <- "H3K4me3_TDH_immune/21";sample.id <- "McGill0104";peaks <- 3
chunk.name <- "H3K4me3_PGP_immune/28";sample.id <- "McGill0005";peaks <- 3
chunk.name <- "H3K4me3_PGP_immune/19";sample.id <- "McGill0005";peaks <- 6
chunk.name <- "H3K4me3_PGP_immune/21";sample.id <- "McGill0322";peaks <- 5
chunk.name <- "H3K4me3_XJ_immune/12";sample.id <- "McGill0101";peaks <- 4

chunk.name <- "H3K4me3_TDH_immune/1";sample.id <- "McGill0025";peaks <- 8
chunk.name <- "H3K4me3_TDH_immune/4";sample.id <- "McGill0106";peaks <- 4
chunk.name <- "H3K4me3_TDH_immune/4";sample.id <- "McGill0005";peaks <- 4
##                                      CDPA PDPA.ignore PDPA.join
## 12:  H3K4me3_TDH_immune/1 McGill0029  7-8         3-5         7
## 13:  H3K4me3_TDH_immune/1 McGill0008    6         3-4         8
## 14:  H3K4me3_TDH_immune/1 McGill0102    8         3-4         8
## 15:  H3K4me3_TDH_immune/1 McGill0107  6-7         6-7         8
## 16:  H3K4me3_TDH_immune/1 McGill0009  3-8         3-4       8-9
## 17:  H3K4me3_TDH_immune/1 McGill0106    7         3-4       8-9
## 18:   H3K4me3_TDH_other/2 McGill0036  8-9         2-4       8-9
## 19:  H3K4me3_TDH_immune/1 McGill0005  3-4         3-4         9

chunk.dir <- file.path("../chip-seq-paper/chunks", chunk.name)
counts.RData <- file.path(chunk.dir, "counts.RData")
load(counts.RData)
counts.dt <- data.table(counts)[sample.id, on=list(sample.id)]
regions.RData <- file.path(chunk.dir, "regions.RData")
load(regions.RData)
regions.dt <- data.table(regions)[sample.id, on=list(sample.id)]

ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
gg <- ggplot()+
  theme_bw()+
  geom_tallrect(aes(
    xmin=chromStart/1e3, xmax=chromEnd/1e3,
    fill=annotation),
    alpha=0.5,
    color="grey",
    size=0.5,
    data=regions.dt)+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart/1e3, coverage),
    color="grey50",
    data=counts.dt)

load("PDPA.infeasible.error.RData")
## load("PDPA.peaks.error.RData")
## load("dp.peaks.error.RData")
## CDPA.error.list <- list()
## chunk.error <- dp.peaks.error[[chunk.name]]
## for(cell.type in names(chunk.error)){
##   type.dt <- data.table(
##     chunk.name,
##     chunk.error[[cell.type]])
##   names(type.dt)[3] <- "peaks"
##   CDPA.error.list[[paste(chunk.name, cell.type)]] <- type.dt
## }
## CDPA.error <- do.call(rbind, CDPA.error.list) 
all.error <- rbind(
##   data.table(algo="CDPA", CDPA.error),
##   data.table(algo="PDPA.ignore", PDPA.peaks.error),
   data.table(algo="PDPA.join", PDPA.infeasible.error))
select.dt <- data.table(
  chunk.name, sample.id, peaks)
all.error[, peaks.int := as.integer(paste(peaks))]
sample.error <- all.error[select.dt, on=list(
    chunk.name, sample.id, peaks.int <= peaks)]

load(file.path(chunk.dir, "PDPA.model.RData"))
fit <- PDPA.model[[sample.id]]
seg.vec <- seq(1, peaks*2+1, by=2)
peaks.list <- list()
segs.list <- list()
for(n.segments in seg.vec){
  n.peaks <- (n.segments-1)/2
  mean.vec <- fit$mean.mat[n.segments, 1:n.segments]
  diff.vec <- diff(mean.vec)
  break.vec <- if(n.segments==1){
    c()
  }else{
    fit$ends.mat[n.segments, 2:n.segments]
  }
  first <- c(1, break.vec+1)
  last <- c(break.vec, nrow(counts.dt))
  peak.dt <- data.table(
    mean=mean.vec,
    first,
    last,
    chromStart=counts.dt$chromStart[first],
    chromEnd=counts.dt$chromEnd[last],
    is.peak=(seq_along(mean.vec)-1) %% 2,
    diff.before=c(Inf, diff.vec),
    diff.after=c(diff.vec, Inf),
    peaks=n.peaks,
    segments=n.segments,
    peak.i=0L)
  segs.list[[paste(n.peaks)]] <- peak.dt
  peaks.list[[paste(n.peaks)]] <- if(n.segments==1){
    data.table()
  }else{
    status.str <- rep(c("background", "peak"), l=n.segments)
    peak.dt[is.peak==0 & (diff.after==0|diff.before==0), is.peak := 1]
    peak.dt[, peak.i := cumsum(is.peak==0)]
    rbind(peak.dt[is.peak==1, list(
      rule="join",
      chromStart=counts.dt$chromStart[min(first)],
      chromEnd=counts.dt$chromEnd[max(last)]),
      by=list(peak.i, peaks, segments)],
      peak.dt[diff.after != 0 & diff.before != 0 & is.peak==1, list(
        rule="remove",
        chromStart,
        chromEnd,
        peak.i=NA, peaks, segments)])
  }
}
peaks.dt <- do.call(rbind, peaks.list)
segs.dt <- do.call(rbind, segs.list)
changes.dt <- segs.dt[1 < first]
changes.dt[, constraint := ifelse(diff.before==0, "equality", "inequality")]

peak.y <- c(
  join=-20,
  remove=-60)
h <- 3
some <- function(dt){
  dt[1<segments]
}
peak.color <- "deepskyblue"
text.dt <- peaks.dt[, .SD[chromEnd==max(chromEnd)], by=list(segments)]
text.dt[, Rule := paste0(toupper(substr(rule, 1, 1)), substr(rule, 2, 100))]
ggm <- ggplot()+
  theme_bw()+
  geom_text(aes(
    138312, peak.y[rule], label=paste0(Rule, " rule")),
    color=peak.color,
    hjust=1,
    size=3,
    data=text.dt)+
  ## geom_tallrect(aes(
  ##   xmin=chromStart/1e3, xmax=chromEnd/1e3,
  ##   fill=annotation),
  ##   alpha=0.5,
  ##   color="grey",
  ##   size=0.5,
  ##   data=regions.dt)+
  ## scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart/1e3, coverage),
    color="grey50",
    data=counts.dt)+
  facet_grid(segments ~ ., labeller=function(df){
    df$segments <- sprintf("$K=%d$
segments", df$segments)
    df
  })+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_vline(aes(
    xintercept=chromStart/1e3,
    linetype=constraint),
    color="blue",
    data=changes.dt)+
  scale_linetype_manual(
    values=c(
      equality="solid",
      inequality="dotted"))+
  geom_segment(aes(
    chromStart/1e3, mean,
    xend=chromEnd/1e3, yend=mean),
    color="blue",
    data=some(segs.dt))+
  geom_segment(aes(
    chromStart/1e3, peak.y[rule],
    xend=chromEnd/1e3, yend=peak.y[rule]),
    color=peak.color,
    size=2,
    data=some(peaks.dt))+
  scale_x_continuous(
    "Position on chromosome")+
  scale_y_continuous(
    "Aligned DNA sequence reads",
    breaks=seq(0, 120, by=40),
    limits=c(-70, 120))
  ## geom_rect(aes(
  ##   xmin=chromStart/1e3, xmax=chromEnd/1e3,
  ##   ymin=peak.y-h, ymax=peak.y+h),
  ##   data=sample.error)
tikz("figure-infeasible-error.tex", 6, 2.75)
print(ggm)
dev.off()
