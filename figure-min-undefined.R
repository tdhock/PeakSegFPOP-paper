source("packages.R")

dt <- data.table(
  count=as.integer(c(13, 14, 10, 1)),
  chromStart=0:3,
  chromEnd=1:4)
pdpa <- PeakSegPDPAchrom(dt, 1L)
cdpa <- PeakSegDP(dt, 1L)

loss.list <- list()
segs.list <- list()

guess <- function(...){
  mean.vec <- c(...)
  mean.name <- paste(mean.vec, collapse=",")
  segs.list[[mean.name]] <<- data.table(mean.name, mean=mean.vec, pos=1:4)
  segs <- do.call(rbind, segs.list)
  Poisson.Loss <- PoissonLoss(dt$count, mean.vec)
  plik <- sum(dpois(dt$count, mean.vec, log=TRUE))
  diff.vec <- diff(mean.vec)
  changes <- sum(diff.vec != 0)
  loss.list[[mean.name]] <<- data.table(
    mean.name, Poisson.Loss,
    segments=factor(changes+1),
    feasible=all(cumsum(sign(diff.vec)) %in% c(0, 1)))
  loss <- do.call(rbind, loss.list)
  loss.ord <- loss[order(Poisson.Loss),]
  loss.ord[, model.i := 1:.N]
  loss.tall <- melt(
    loss.ord,
    id.vars=c("mean.name", "model.i", "feasible", "segments"))
  viz <- list(
    data=ggplot()+
      theme_bw()+
      geom_text(aes(
        2.5, 0,
        label=sprintf("Poisson Loss = %f", Poisson.Loss),
        showSelected=mean.name),
        data=loss.ord)+
      geom_segment(aes(
        pos-0.5, mean, xend=pos+0.5, yend=mean, showSelected=mean.name),
        color="green",
        size=1,
        data=segs)+
      geom_point(aes(chromEnd, count),
                 size=4,
                 data=dt)+
      xlab("position on chromosome")+
      theme(panel.grid.minor=element_blank())+
      scale_y_continuous("align read counts", breaks=dt$count),
    loss=ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(variable ~ ., scales="free")+
      geom_tallrect(aes(
        xmin=model.i-0.5, xmax=model.i+0.5, clickSelects=mean.name),
        alpha=0.5,
        data=loss.ord)+
      scale_color_manual(values=c("TRUE"="white", "FALSE"="black"))+
      guides(colour=guide_legend(override.aes=list(fill="grey")))+
      geom_point(aes(model.i, value, color=feasible, fill=segments),
                 size=4,
                 shape=21,
                 data=loss.tall),
    first=list(mean.name=mean.name))
  print(viz$loss)
  animint2dir(viz, "figure-min-undefined")
}

dist.vec <- 10^seq(0, -3, by=-1)
mid <- 37/3
for(fig.i in seq_along(dist.vec)){
  dist.to.opt <- dist.vec[[fig.i]]
  lo <- mid-dist.to.opt
  hi <- mid+dist.to.opt
  mean.vec <- c(lo,lo,hi,1)
  Poisson.Loss <- PoissonLoss(dt$count, mean.vec)
  loss <- data.table(Poisson.Loss)
  segs <- data.table(mean=mean.vec, pos=1:4)
  mean.text <- data.table(
    mean=c(lo, hi),
    x=c(1.5, 3),
    vjust=c(1.5, -0.5))
  
  gg <- ggplot()+
    theme_bw()+
    geom_text(aes(
      x, mean, label=sprintf("%.3f", mean), vjust=vjust),
      data=mean.text)+
    geom_text(aes(
      2.5, 0,
      label=sprintf("Poisson Loss = %f", Poisson.Loss)),
      data=loss)+
    geom_segment(aes(
      pos-0.5, mean, xend=pos+0.5, yend=mean),
      color="green",
      size=1,
      data=segs)+
    geom_point(aes(chromEnd, count),
               size=4,
               data=dt)+
    xlab("position on chromosome")+
    theme(panel.grid.minor=element_blank())+
    scale_y_continuous("align read counts", breaks=dt$count)

  pdf(print(sprintf("figure-min-undefined-%d.pdf", fig.i)), w=5, h=3)
  print(gg)
  dev.off()
  
}

mean.vec <- c(mid, mid, mid, ,1)
Poisson.Loss <- PoissonLoss(dt$count, mean.vec)
loss <- data.table(Poisson.Loss)
segs <- data.table(mean=mean.vec, pos=1:4)
mean.text <- data.table(
  mean=c(mid, mid),
  x=c(1.5, 3),
  vjust=c(1.5, -0.5))
gg <- ggplot()+
  theme_bw()+
  geom_text(aes(
    x, mean, label=sprintf("%f", mean), vjust=vjust),
    data=mean.text)+
  geom_text(aes(
    2.5, 0,
    label=sprintf("Poisson Loss = %f", Poisson.Loss)),
    data=loss)+
  geom_segment(aes(
    pos-0.5, mean, xend=pos+0.5, yend=mean),
    color="green",
    size=1,
    data=segs)+
  geom_point(aes(chromEnd, count),
             size=4,
             data=dt)+
  xlab("position on chromosome")+
  theme(panel.grid.minor=element_blank())+
  scale_y_continuous("align read counts", breaks=dt$count)
pdf("figure-min-undefined.pdf", w=5, h=3)
print(gg)
dev.off()
