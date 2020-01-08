library(data.table)
library(PeakSegDP)
library(PeakError)
library(Segmentor3IsBack)

data(chr11ChIPseq)

sample.id <- "McGill0322"
coverage.list <- split(chr11ChIPseq$coverage, chr11ChIPseq$coverage$sample.id)
compressed <- coverage.list[[sample.id]]
compressed$bases <- with(compressed, chromEnd-chromStart)
count <- with(compressed, rep(count, bases))
base <- with(compressed, (chromStart[1]+1):chromEnd[length(chromEnd)])
counts <- data.frame(base, count)

max.segments <- 7
maxPeaks <- as.integer((max.segments-1)/2)
fit <- Segmentor(count, Kmax=max.segments)
Segmentor.segs <- NULL
Segmentor.breaks <- NULL
for(model.i in seq(1, max.segments, by=2)){
  peaks <- (model.i-1)/2
  last.i <- fit@breaks[model.i, 1:model.i]
  break.after <- last.i[-model.i]
  first.i <- c(1, break.after+1)
  param <- fit@parameters[model.i, 1:model.i]
  for(segment.i in seq_along(first.i)){
    from <- first.i[[segment.i]]
    to <- last.i[[segment.i]]
    seg.mean <- param[[segment.i]]
    y <- count[from:to]
  }
  signs <- as.integer(sign(diff(param)))
  feasible <- all(cumsum(signs) %in% c(0,1))
  meta <-
    data.frame(peaks,
               feasible,
               segments=model.i)
  if(length(break.after)){
    Segmentor.breaks <- rbind(Segmentor.breaks, {
      data.frame(meta, chromEnd=base[break.after])
    })
  }
  Segmentor.segs <- rbind(Segmentor.segs, {
    data.frame(meta,
               chromStart=base[first.i],
               chromEnd=base[last.i],
               mean=param,
               row.names=NULL)
  })
}

dp.fit <- PeakSegDP(compressed, maxPeaks=maxPeaks)
pdpa <- PeakSegOptimal::PeakSegPDPAchrom(compressed, maxPeaks)

cfac <- function(x)factor(x, c("unconstrained", "constrained"))
seg.cols <- c("segments", "chromStart", "chromEnd", "mean")
both.segs <-
  rbind(data.frame(Segmentor.segs[,seg.cols], model=cfac("unconstrained")),
        data.frame(dp.fit$segments[,seg.cols], model=cfac("constrained")))
break.cols <- c("segments", "chromEnd")
both.breaks <-
  rbind(data.frame(Segmentor.breaks[,break.cols], model=cfac("unconstrained")),
        data.frame(dp.fit$breaks[,break.cols], model=cfac("constrained")))

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

## Confront the PeakSegDP model with the annotated regions.
region.list <- split(chr11ChIPseq$regions, chr11ChIPseq$regions$sample.id)
regions <- region.list[[sample.id]]
error.list <- list()
for(peaks in names(dp.fit$peaks)){
  peak.df <- dp.fit$peaks[[peaks]]
  segments <- as.integer(peaks)*2 + 1
  error.list[[peaks]] <-
    data.frame(PeakErrorChrom(peak.df, regions),
               peaks, segments, model=cfac("constrained"))
}
error.regions <- do.call(rbind, error.list)

library(ggplot2)
segs.regions <-
  ggplot()+
  scale_fill_manual("label", values=ann.colors, 
                    breaks=names(ann.colors))+
  scale_linetype_manual("error type", 
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  scale_y_continuous(breaks=seq(0, 50, by=25), minor_breaks=NULL)+
  scale_x_continuous("position on chr11 (kilo bases)",
                     breaks=seq(118080, 118120, by=20))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(segments ~ model, scales="free")+
  geom_step(aes(chromStart/1e3, count),
            data=compressed, color="grey40")+
  geom_segment(aes((chromStart-1/2)/1e3, mean,
                   xend=(chromEnd+1/2)/1e3, yend=mean),
               data=both.segs, color="green", alpha=3/4, size=1)+
  geom_vline(aes(xintercept=(chromEnd+1/2)/1e3),
             data=both.breaks,
             color="green")+
  ylab("count of aligned reads")

png("figure-Segmentor-PeakSeg-noregions.png",
    units="in", res=200, width=6, height=3)
print(segs.regions)
dev.off()

three.segs <- data.table(both.segs)[model=="constrained" & segments==3]
three.segs[, seg.i := 1:.N]
three.breaks <- data.table(both.breaks)[model=="constrained" & segments==3]
gg <- ggplot()+
  scale_fill_manual("label", values=ann.colors, 
                    breaks=names(ann.colors))+
  scale_linetype_manual("error type", 
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  scale_y_continuous(
    ##"read count",
    "$z_i, u_k$",
    breaks=seq(0, 50, by=25),
    minor_breaks=NULL)+
  scale_x_continuous(
    "",
    ##"position on chr11 (kilo bases)",
    breaks=three.segs[, c(chromStart, chromEnd[.N])],
    labels=sprintf("$t_%s$", paste0(0:3, c("=0", "", "", "=n")))
  )+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_step(aes(chromStart, count),
            data=compressed, color="grey40")+
  geom_segment(aes((chromStart-1/2), mean,
                   xend=(chromEnd+1/2), yend=mean),
               data=three.segs,
               color="green", alpha=3/4, size=2)+
  geom_text(aes(
    ifelse(seg.i==1, chromStart, chromEnd),
                mean,
                hjust=ifelse(seg.i==1, 1, 0),
    label=sprintf("$u_%d$", seg.i)),
    vjust=0,
            color="green",
    data=three.segs)+
  geom_text(aes(x, y, label="$S=3$ segments"),
            data=data.frame(x=118090000, y=50),
            hjust=0,
            size=3,
            color="green")+
  geom_vline(aes(
    xintercept=chromEnd),
    data=three.breaks,
    color="green",
    size=0.5)+
  geom_text(aes(
    ifelse(seg.i==1, chromStart, chromEnd),
    30,
    label=sprintf("$z_%s$", ifelse(seg.i==1, 1, "n"))),
    data=three.segs[seg.i != 2])
library(tikzDevice)
tikz("figure-PeakSeg.tex", width=4, height=1)
print(gg)
dev.off()
gg <- ggplot()+
  scale_fill_manual("label", values=ann.colors, 
                    breaks=names(ann.colors))+
  scale_linetype_manual("error type", 
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  scale_y_continuous(
    ##"read count",
    "Data $z_i$
Mean $u_k$
aligned DNA
sequences",
    breaks=seq(0, 50, by=25),
    minor_breaks=NULL)+
  scale_x_continuous(
    "Position on chromosome",
    ##"position on chr11 (kilo bases)",
    breaks=three.segs[, c(chromStart, chromEnd[.N])],
    labels=sprintf("$t_%s$", paste0(0:3, c("=0", "", "", "=N")))
  )+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_step(aes(chromStart, count),
            data=compressed, color="grey40")+
  geom_segment(aes((chromStart-1/2), mean,
                   xend=(chromEnd+1/2), yend=mean),
               data=three.segs,
               color="blue", alpha=3/4, size=2)+
  geom_text(aes(
    ifelse(seg.i==1, chromStart, chromEnd),
                mean,
                hjust=ifelse(seg.i==1, 1, 0),
    label=sprintf("$u_%d$", seg.i)),
    vjust=0,
            color="blue",
    data=three.segs)+
  geom_text(aes(x, y, label="$K=3$ segments"),
            data=data.frame(x=118090000, y=60),
            hjust=0,
            size=3,
            color="blue")+
  geom_vline(aes(
    xintercept=chromEnd),
    data=three.breaks,
    color="blue",
    linetype="dashed",
    size=0.5)+
  geom_text(aes(
    ifelse(seg.i==1, chromStart, chromEnd),
    20,
    label=sprintf("$z_%s$", ifelse(seg.i==1, 1, "N"))),
    color="grey40",
    data=three.segs[seg.i != 2])
tikz("figure-PeakSeg-big.tex", width=6, height=1.5)
print(gg)
dev.off()

five.segs <- data.table(both.segs)[model=="unconstrained" & segments==5]
three.breaks <- data.table(both.breaks)[model=="unconstrained" & segments==5]
lab.loc <- data.table(x=118090, y=50)
gg <- ggplot()+
  scale_fill_manual("label", values=ann.colors, 
                    breaks=names(ann.colors))+
  scale_linetype_manual("error type", 
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  scale_y_continuous(
    ##"read count",
    "$z_i, u_k$",
    breaks=seq(0, 50, by=25),
    minor_breaks=NULL)+
  scale_x_continuous(
    "", breaks=NULL
  )+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_step(aes(chromStart/1e3, count),
            data=compressed, color="grey40")+
  geom_segment(aes((chromStart-1/2)/1e3, mean,
                   xend=(chromEnd+1/2)/1e3, yend=mean),
               data=five.segs,
               color="green", alpha=3/4, size=2)+
  geom_text(aes(x, y, label="$S=5$ segments\nunconstrained"),
            data=lab.loc,
            hjust=0,
            size=3,
            color="green")+
  geom_vline(aes(
    xintercept=chromEnd/1e3),
    data=three.breaks,
    color="green",
    size=0.5)
tikz("figure-PeakSeg-unconstrained.tex", width=4, height=1)
print(gg)
dev.off()

five.segs <- data.table(both.segs)[model=="constrained" & segments==5]
three.breaks <- data.table(both.breaks)[model=="constrained" & segments==5]
gg <- ggplot()+
  scale_fill_manual("label", values=ann.colors, 
                    breaks=names(ann.colors))+
  scale_linetype_manual("error type", 
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  scale_y_continuous(
    ##"read count",
    "$z_i, u_k$",
    breaks=seq(0, 50, by=25),
    minor_breaks=NULL)+
  scale_x_continuous(
    "", breaks=NULL
  )+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_step(aes(chromStart/1e3, count),
            data=compressed, color="grey40")+
  geom_segment(aes((chromStart-1/2)/1e3, mean,
                   xend=(chromEnd+1/2)/1e3, yend=mean),
               data=five.segs,
               color="green", alpha=3/4, size=2)+
  geom_segment(aes(
    (chromStart-1/2)/1e3, 70,
    xend=(chromEnd+1/2)/1e3, yend=70),
    data=five.segs[seq_along(mean) %in% c(2, 4)],
    color="deepskyblue", alpha=3/4, size=4)+
  geom_text(aes(x, y, label="$S=5$ segments\nconstrained"),
            data=lab.loc,
            hjust=0,
            size=3,
            color="green")+
  geom_vline(aes(
    xintercept=chromEnd/1e3),
    data=three.breaks,
    color="green",
    size=0.5)
tikz("figure-PeakSeg-constrained.tex", width=4, height=1)
print(gg)
dev.off()


