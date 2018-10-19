##          chunk.name  sample.id peaks      PDPA cDPA.forward cDPA.reverse
## H3K4me3_XJ_immune/1 McGill0024     9  -25744.2    -25903.03    -25903.03
source("packages.R")
load("../chip-seq-paper/chunks/H3K4me3_XJ_immune/1/counts.RData")
(objs <- load("../chip-seq-paper/chunks/H3K4me3_XJ_immune/1/PDPA.model/McGill0024.RData"))
input.dt <- data.table(counts)[sample.id=="McGill0024",]
input.dt[, weight := chromEnd-chromStart]
input.dt[, count := coverage]
(objs <- load("../chip-seq-paper/chunks/H3K4me3_XJ_immune/1/dp.model.RData"))
(loss.mat <- rbind(
  cDPA=dp.model$McGill0024$error$error,
  PDPA=result$model$models$min.cost,
  min=NA))

ggplot()+
  geom_step(aes(chromStart/1e3, count),
            color="grey50",
            data=input.dt)
(objs <- load("problematic.some.RData"))

getMinMean <- function(dt){
  dt[, -Log/Linear]
}
ploss <- function(dt, x){
  ## need to make a new data table, otherwise ifelse may only get one
  ## element, and return only one element.
  new.dt <- data.table(dt, x)
  new.dt[, ifelse(Log==0, 0, Log*log(x)) + Linear*x + Constant]
}
Minimize <- function(dt, from=min(dt$min.mean), to=max(dt$max.mean)){
  stopifnot(from < to)
  is.before <- dt$max.mean < from
  is.after <- to < dt$min.mean
  feasible <- dt[!(is.before | is.after),]
  feasible$min.mean[1] <- from
  feasible$max.mean[nrow(feasible)] <- to
  feasible$fun.min.mean <- getMinMean(feasible)
  feasible[, min.cost.mean := ifelse(
    fun.min.mean < min.mean, min.mean,
    ifelse(max.mean < fun.min.mean,
           max.mean,
           fun.min.mean))]
  feasible[, min.cost := ploss(feasible, min.cost.mean)]
  feasible[which.min(min.cost),]
}

timestep <- nrow(input.dt)
for(col.i in 1:10){
  peaks <- col.i-1
  total.segments <- peaks*2+1
  cost.model <- cost.models.list[[paste(total.segments, timestep)]]
  loss.mat["min", col.i] <- Minimize(cost.model)$min.cost
}
loss.mat

fit <- cDPA(input.dt$count, input.dt$weight, maxSegments=19)

PDPA.loss <- matrix(NA, 19, 376)
for(total.segments in 1:19){
  for(timestep in total.segments:nrow(input.dt)){
    cost.model <- cost.models.list[[paste(total.segments, timestep)]]
    PDPA.loss[total.segments, timestep] <- Minimize(cost.model)$min.cost
  }
}

scatter.dt.list <- list()
for(total.segments in 1:nrow(PDPA.loss)){
  scatter.dt.list[[paste(total.segments)]] <- data.table(
    total.segments,
    cDPA=fit$loss[total.segments,],
    PDPA=PDPA.loss[total.segments,],
    data.i=1:ncol(PDPA.loss))
}
scatter.dt <- do.call(rbind, scatter.dt.list)

scatter.dt[, should.be.positive := cDPA-PDPA]
scatter.dt[, problem := ifelse(should.be.positive > 0, 0, should.be.positive)]
abline.dt <- data.table(slope=1, intercept=0)
ggplot()+
  geom_point(aes(cDPA, PDPA, color=problem), shape=1, data=scatter.dt)+
  theme_bw()+
  facet_wrap("total.segments")+
  theme(panel.margin=grid::unit(0, "lines"))+
  coord_equal()+
  geom_abline(aes(slope=slope, intercept=intercept), data=abline.dt)

ggplot()+
  geom_point(aes(cDPA, PDPA, color=problem), shape=1,
             data=scatter.dt[total.segments==14,])+
  theme_bw()+
  facet_wrap("total.segments")+
  theme(panel.margin=grid::unit(0, "lines"))+
  coord_equal()+
  geom_abline(aes(slope=slope, intercept=intercept), data=abline.dt)

scatter.dt[total.segments==14 & PDPA > cDPA,]
scatter.dt[total.segments==12 & PDPA > cDPA,]
scatter.dt[PDPA > cDPA,]
scatter.dt[should.be.positive < -1e-10,]

## After examining MinEnv gg.pruning plots in
## figure-constrained-PDPA-poisson-real.R, it is clear that the bug in
## the MinEnv computation starts at total.segments=14, timestep=299.

data.lines.list <- list()
data.minima.list <- list()
data.infeasible.list <- list()
data.cost.list <- list()
data.intervals.list <- list()
## for animint:
for(total.segments in 1:max.segments){
  for(timestep in total.segments:length(input.dt$count)){
    cat(sprintf("decoding %4d / %4d segments %4d / %4d data points\n",
                total.segments, max.segments, timestep, length(input.dt$count)))
    data.i <- timestep
    seg.i <- total.segments
    no.constraint <- data.table(
      min.mean,
      max.mean,
      data.i=NA)
    constraint <- no.constraint
    segment.end <- timestep
    while(0 < seg.i && length(data.i)==1){
      unconstrained.fun <- cost.models.list[[paste(seg.i, data.i)]]
      if(seg.i==total.segments){
        data.intervals.list[[paste(total.segments, timestep)]] <- data.table(
          total.segments, timestep, intervals=nrow(unconstrained.fun))
      }
      if(nrow(unconstrained.fun)){
        show.lines <- getLines(unconstrained.fun)
        ##if(seg.i>1)show.lines$data.i <- show.lines$data.i+1L
        if(seg.i==1)show.lines$data.i <- 0
        data.lines.list[[paste(total.segments, timestep, seg.i)]] <-
          data.table(total.segments, timestep, seg.i,
                     show.lines)
      }
      min.dt <- Minimize(
        unconstrained.fun,
        constraint$min.mean,
        constraint$max.mean)
      if(seg.i==total.segments){
        data.cost.list[[paste(total.segments, timestep)]] <-
          data.table(total.segments, timestep,
                     optimal.cost=min.dt$min.cost,
                     constraint="inactive")
      }
      ggplot()+
        geom_point(aes(min.cost.mean, min.cost),
                   data=min.dt)+
        geom_line(aes(mean, cost, color=factor(data.i)),
                  data=getLines(unconstrained.fun))
      min.dt$segment.end <- segment.end
      min.dt[, segment.start := ifelse(seg.i==1, 1, 1+data.i)]
      segment.end <- min.dt$data.i
      if(min.mean < constraint$min.mean){
        data.infeasible.list[[paste(
          total.segments, timestep, seg.i)]] <-
          data.table(total.segments, timestep, seg.i,
                     min.mean=-Inf,
                     max.mean=constraint$min.mean)
      }
      if(constraint$max.mean < max.mean){
        data.infeasible.list[[paste(
          total.segments, timestep, seg.i)]] <-
          data.table(total.segments, timestep, seg.i,
                     min.mean=constraint$max.mean,
                     max.mean=Inf)
      }
      min.dt$constraint <- if(min.dt[, fun.min.mean == min.cost.mean]){
        "inactive"
      }else{
        data.cost.list[[paste(
          total.segments, timestep)]]$constraint <- "active"
        "active"
      }
      show.min <- data.table(
        total.segments, timestep, seg.i,
        min.dt)
      if(seg.i==1)show.min$data.i <- 0
      data.minima.list[[paste(total.segments, timestep, seg.i)]] <-
        show.min
      constraint <- no.constraint
      constraint.side <- if(seg.i %% 2){
        constraint$min.mean <- min.dt$min.cost.mean
      }else{
        constraint$max.mean <- min.dt$min.cost.mean
      }
      data.i <- min.dt$data.i
      seg.i <- seg.i-1
    }#while(...
  }#for(timestep
}#for(total.segments

data.intervals <- do.call(rbind, data.intervals.list)
data.dt <- data.table(
  count=input.dt$count,
  position=seq_along(input.dt$count),
  data.i.fac=factor(seq_along(input.dt$count)))
data.lines <- do.call(rbind, data.lines.list)
data.lines[, data.i.fac := factor(data.i)]
data.lines[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
left.of.intervals <-
  data.lines[, .SD[1,],
             by=.(total.segments, timestep, minimization, seg.i, piece.i)]
between.intervals <- left.of.intervals[min.mean != min(input.dt$count),]
data.minima <- do.call(rbind, data.minima.list)
data.minima[, data.i.fac := factor(data.i)]
data.minima[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
data.infeasible <- do.call(rbind, data.infeasible.list)
data.infeasible[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
data.cost <- do.call(rbind, data.cost.list)
data.cost[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
addY <- function(dt, y){
  data.table(dt, y=factor(y, c("data value", "intervals", "segments")))
}
largest.constant <- envelope[Linear==0, max(Constant)]
viz <- list(
  title=paste(
    "Constrained Pruned Dynamic Programming Algorithm"),
  funModels=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=800, height=300)+
    ##coord_cartesian(ylim=c(0, max(between.intervals$cost)))+
    geom_line(aes(mean, cost,
                  key=paste(timestep, total.segments),
                  showSelected=total.segments, showSelected2=timestep),
              color="grey",
              size=8,
              data=data.table(envelope, seg.i="pruning"))+
    geom_line(aes(mean, cost, color=data.i.fac,
                  group=paste(piece.i, data.i),
                  key=paste(cost.type, min.mean, max.mean),
                  showSelected=total.segments, showSelected2=timestep),
              data=data.table(cost.lines, seg.i="pruning"))+
    ## geom_point(aes(min.cost.mean, min.cost, color=data.i.fac,
    ##                showSelected=total.segments, showSelected2=timestep),
    ##            size=5,
    ##            data=data.table(minima, seg.i="pruning"))+    
    facet_grid(. ~ seg.i, scales="free", labeller=function(var, val){
      paste(ifelse(val!="pruning", "segment", ""), val)
    })+
    geom_tallrect(aes(xmin=min.mean, xmax=max.mean,
                      key=timestep,
                      showSelected=total.segments,
                      showSelected2=timestep),
                  fill="grey",
                  alpha=0.5,
                  color=NA,
                  data=data.infeasible)+
    geom_line(aes(mean, cost, color=data.i.fac,
                  group=piece.i,
                  key=paste(min.mean, max.mean),
                  showSelected=total.segments,
                  showSelected2=timestep),
              data=data.lines)+
    guides(color="none"),
  data=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=800, height=300)+
    facet_grid(y ~ ., scales="free")+
    geom_point(aes(position, count),
               data=addY(data.dt, "data value"))+
    geom_vline(aes(xintercept=segment.start-0.5,
                   key=seg.i,
                   showSelected=total.segments,
                   showSelected2=timestep),
               data=addY(data.minima[1<segment.start,], "data value"),
               color="green",
               linetype="dashed")+
    geom_segment(aes(segment.start-0.45, min.cost.mean,
                     showSelected=total.segments,
                     showSelected2=timestep,
                     key=seg.i,
                     xend=segment.end+0.45, yend=min.cost.mean),
                 data=addY(data.minima, "data value"),
                 color="green")+
    guides(color="none")+
    geom_tallrect(aes(xmin=timestep-0.5, xmax=timestep+0.5,
                      clickSelects=timestep),
                  data=data.table(timestep=seq_along(input.dt$count)),
                  alpha=0.5)+
    geom_line(aes(timestep, intervals, group=total.segments,
                  clickSelects=total.segments),
               data=addY(data.intervals, "intervals"))+
    geom_tile(aes(timestep, total.segments, fill=optimal.cost),
              data=addY(data.cost, "segments"))+
    geom_widerect(aes(ymin=total.segments-0.5, ymax=total.segments+0.5,
                      clickSelects=total.segments),
                  alpha=0.5,
                  data=addY(
                    data.table(total.segments=1:max.segments), "segments"))+
    ylab("")+
    scale_x_continuous(
      "data point",
      breaks=unique(c(seq(1, length(input.dt$count), by=10), length(input.dt$count)))),
  time=list(variable="timestep", ms=2000),
  duration=list(timestep=2000)
)
minima.active <- data.minima[constraint=="active",]
if(nrow(minima.active)){
  viz$funModels <- viz$funModels+
    geom_point(aes(min.cost.mean, min.cost,
                   key=data.i,
                   showSelected=total.segments, showSelected2=timestep),
               size=6,
               shape=21,
               fill="white",
               color="black",
               data=minima.active)
}
viz$funModels <- viz$funModels+
  geom_point(aes(min.cost.mean, min.cost, color=data.i.fac,
                 tooltip=paste(
                   "minimum cost =",
                   round(min.cost, 4),
                   "with",
                   constraint,
                   "constraint at mean =",
                   round(min.cost.mean, 4),
                   "for",
                   seg.i,
                   "segment model up to data point",
                   segment.end,
                   "previous segment end =",
                   data.i
                 ),
                 key=data.i,
                 showSelected=total.segments, showSelected2=timestep),
             size=5,
             data=data.minima)+
  geom_point(aes(mean, cost,
                 key=mean,
                 showSelected=total.segments, showSelected2=timestep),
             data=between.intervals)  
cost.active <- data.cost[constraint=="active",]
if(nrow(cost.active)){
  viz$data <- viz$data+
    geom_point(aes(timestep, total.segments),
               shape=21,
               color="black",
               fill="white",
               data=addY(cost.active, "segments"))
}
animint2dir(viz, "figure-constrained-PDPA-poisson-real")

intervalsPlot <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(total.segments ~ .)+
  geom_point(aes(timestep, intervals),
             data=data.intervals)

pdf("figure-constrained-PDPA-poisson-real.pdf")
print(intervalsPlot)
dev.off()

## FunctionalPruning <- list(
##   envelope=data.frame(envelope),
##   cost.lines=data.frame(cost.lines),
##   minima=data.frame(minima),
##   grid=data.frame(data.cost))
## save(FunctionalPruning, file="~/R/animint/data/FunctionalPruning.RData")
## prompt(FunctionalPruning, file="~/R/animint/man/FunctionalPruning.Rd")
