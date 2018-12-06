source("packages.R")

data(chr11ChIPseq)

sample.id <- "McGill0322"
coverage.list <- split(chr11ChIPseq$coverage, chr11ChIPseq$coverage$sample.id)
compressed <- coverage.list[[sample.id]]
compressed$bases <- with(compressed, chromEnd-chromStart)
count <- with(compressed, rep(count, bases))
i.vec <- with(compressed, rep(seq_along(bases), bases))
base <- with(compressed, (chromStart[1]+1):chromEnd[length(chromEnd)])
counts <- data.frame(base, count)

max.segments <- 5
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
               first=i.vec[first.i],
               last=i.vec[last.i],
               chromStart=base[first.i],
               chromEnd=base[last.i],
               mean=param,
               row.names=NULL)
  })
}

dp.fit <- PeakSegDP(compressed, maxPeaks=maxPeaks)
pdpa.fit <- PeakSegPDPAchrom(compressed, 2L)
pdpa.fit$loss
dp.fit$error

cfac <- function(x)factor(x, c("unconstrained", "constrained"), c("Unconstrained", "Up-down constrained"))
seg.cols <- c("segments", "chromStart", "chromEnd", "mean", "first", "last")
both.segs <-
  rbind(data.frame(Segmentor.segs[,seg.cols], model=cfac("unconstrained")),
        data.frame(dp.fit$segments[,seg.cols], model=cfac("constrained")))
break.cols <- c("segments", "chromEnd")
both.breaks <-
  rbind(data.frame(Segmentor.breaks[,break.cols], model=cfac("unconstrained")),
        data.frame(dp.fit$breaks[,break.cols], model=cfac("constrained")))

n.segs <- 5
show.segs <- data.table(both.segs)[segments==n.segs]
setkey(show.segs, model)
show.breaks <- data.table(both.breaks)[segments==n.segs]
show.breaks[, feasible := {
  m <- model
  show.segs[J(m), {
    ifelse(sign(diff(mean)) == c(1, -1), "yes", "no")
  }]
}, by=list(model)]
lik.dt <- show.segs[, {
  compressed$mean <- rep(mean, last-first+1)
  lik.vec <- with(compressed, dpois(count, mean, log=TRUE))
  data.table(
    log.lik=sum(lik.vec * compressed$bases),
    loss=with(compressed, PoissonLoss(count, mean, bases)))
}, by=list(model)]
model.color <- "blue"
segs.regions <-
  ggplot()+
  scale_linetype_manual("feasible\nfor up-down\nconstraint", 
                        values=c(yes="dotted", no="solid"))+
  scale_y_continuous(breaks=seq(0, 50, by=25), minor_breaks=NULL)+
  scale_x_continuous("position on chr11 (kilo bases)",
                     breaks=seq(118080, 118120, by=20))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(model ~ .)+
  geom_text(aes(
    118100, 50,
    label=sprintf("LogLik=%.1f", log.lik)),
            data=lik.dt)+
  geom_step(aes(chromStart/1e3, count),
            data=compressed, color="grey40")+
  geom_segment(aes((chromStart-1/2)/1e3, mean,
                   xend=(chromEnd+1/2)/1e3, yend=mean),
               data=show.segs,
               color=model.color, alpha=3/4, size=1)+
  geom_vline(aes(xintercept=(chromEnd+1/2)/1e3, linetype=feasible),
             data=show.breaks,
             color=model.color)+
  ylab("count of aligned reads")
print(segs.regions)
png("figure-data-models.png",
    units="in", res=200, width=6, height=3.5)
print(segs.regions)
dev.off()

