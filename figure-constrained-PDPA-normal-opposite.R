source("packages.R")

## from the other direction, 13 14 10 1

## The first comparison is C12< versus C11<+gamma_2
data.vec <- rev(c(1, 10, 14, 13))
min.mean <- min(data.vec)
max.mean <- max(data.vec)
gamma.dt <- data.table(
  quadratic=1,
  linear=-2*data.vec,
  constant=data.vec * data.vec)
C1.dt <- cumsum(gamma.dt)
gamma.dt$min.mean <- C1.dt$min.mean <- min.mean
gamma.dt$max.mean <- C1.dt$max.mean <- max.mean
gamma.dt$data.i <- C1.dt$data.i <- seq_along(data.vec)
quad <- function(dt, x){
  dt[, quadratic*x*x + linear*x + constant]
}
getLines <- function(dt.name){
  dt <- get(dt.name)
  line.list <- list()
  for(piece.i in 1:nrow(dt)){
    piece <- dt[piece.i,]
    mean.vec <- piece[, seq(min.mean, max.mean, l=100)]
    line.list[[piece.i]] <- data.table(
      piece.i,
      piece,
      mean=mean.vec,
      cost=quad(piece, mean.vec))
  }
  data.table(dt.name, do.call(rbind, line.list))
}
getMinMean <- function(dt){
  dt[, -linear/(2*quadratic)]
}
getLessEqualMin <- function(dt){
  if(1 < nrow(dt)){
    stop("TODO implement more general less min computation")
  }
  mu <- getMinMean(dt)
  cost <- quad(dt, mu)
  data.table(
    quadratic=c(dt$quadratic, 0),
    linear=c(dt$linear, 0),
    constant=c(dt$constant, cost),
    min.mean=c(min.mean, mu),
    max.mean=c(mu, max.mean))
}
getLessMin <- function(dt){
  if(1 < nrow(dt)){
    stop("TODO implement more general less min computation")
  }
  mu <- getMinMean(dt)
  cost <- quad(dt, mu)
  data.table(
    quadratic=0,
    linear=0,
    constant=cost,
    min.mean=mu,
    max.mean=dt$max.mean,
    data.i=NA)
}
getMoreMin <- function(dt){
  if(1 < nrow(dt)){
    stop("TODO implement more general more min computation")
  }
  mu <- getMinMean(dt)
  ## TODO check if mu is inside the interval. If it is not then return
  ## nothing (representing infinite cost on the entire interval).
  cost <- quad(dt, mu)
  data.table(
    quadratic=0,
    linear=0,
    constant=cost,
    min.mean=mu,
    max.mean=dt$max.mean,
    data.i=NA)
}
AddFuns <- function(dt1, dt2){
  mean.vec <- sort(unique(rbind(dt1, dt2)[, c(min.mean, max.mean)]))
  i1 <- 1
  i2 <- 1
  new.dt.list <- list()
  while(i1 <= nrow(dt1) && i2 <= nrow(dt2)){
    row1 <- dt1[i1,]
    row2 <- dt2[i2,]
    this.min <- max(row1$min.mean, row2$min.mean)
    this.max <- if(row1$max.mean < row2$max.mean){
      i1 <- i1+1
      row1$max.mean
    }else{
      i2 <- i2+1
      row2$max.mean
    }
    new.dt.list[[paste(i1, i2)]] <- data.table(
      quadratic=row1$quadratic+row2$quadratic,
      linear=row1$linear+row2$linear,
      constant=row1$constant+row2$constant,
      min.mean=this.min,
      max.mean=this.max,
      data.i=NA)
    this.min <- this.max
  }
  do.call(rbind, new.dt.list)
}

C11.less <- getLessMin(C1.dt[1,])
C22 <- AddFuns(gamma.dt[2,], C11.less)
C12.less <- getLessMin(C1.dt[2,])

##dput(RColorBrewer::brewer.pal(Inf, "Set1"))
break.colors <- c("1"="#E41A1C",
  "2"="#377EB8",
  "#4DAF4A",
  "3"="#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999")

first.choice.lines <- rbind(
  data.table(break.after=factor(2), getLines("C12.less")),
  data.table(break.after=factor(1), getLines("C22")))
first.choice.ends <- first.choice.lines[mean==min.mean | mean==max.mean,]
ggplot()+
  scale_x_continuous(breaks=first.choice.ends[, unique(c(min.mean, max.mean))])+
  scale_color_manual(values=break.colors)+
  geom_point(aes(mean, cost, group=dt.name, color=break.after),
            data=first.choice.ends)+
  geom_line(aes(mean, cost, group=dt.name, color=break.after),
            data=first.choice.lines)
## C22 always has a lower cost than C12.less, so prune. This means
## that we discard the possibility of a break after 2.

C23 <- AddFuns(C22, gamma.dt[3,])
C23.lines <- rbind(
  data.table(break.after=factor(1), getLines("C23")))
C23.ends <- C23.lines[mean==min.mean | mean==max.mean,]
breaks.vec <- C23.ends[, unique(c(min.mean, max.mean))]
min.cost.candidates <- data.table(break.after=factor(1), C23)
min.cost.candidates$min.cost.mean <- getMinMean(C23)
min.cost.mean.dt <-
  min.cost.candidates[min.mean < min.cost.mean & min.cost.mean < max.mean,]
min.cost.mean.dt$min.cost <-
  quad(min.cost.mean.dt, min.cost.mean.dt$min.cost.mean)
ggplot()+
  ##coord_cartesian(ylim=c(0,20))+
  scale_x_continuous(breaks=breaks.vec)+
  scale_color_manual(values=break.colors)+
  geom_line(aes(mean, cost, color=break.after),
            data=C23.lines)+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(C23.ends, what="limit"))+
  ## geom_point(aes(min.cost.mean, min.cost, color=break.after, fill=what),
  ##            shape=21,
  ##            data=data.table(min.cost.mean.dt, what="min"))+
  scale_fill_manual(values=c(min="black", limit=NA))

C13.less <- getLessMin(C1.dt[3,])
second.choice.lines <- rbind(
  data.table(break.after=factor(3), getLines("C13.less")),
  data.table(break.after=factor(1), getLines("C23")))
second.choice.ends <- second.choice.lines[mean==min.mean | mean==max.mean,]
breaks.vec <- second.choice.ends[, unique(c(min.mean, max.mean))]
discriminant <- C23$linear^2 - 4*C23$quadratic*(C23$constant-C13.less$constant)
numerator <- -C23$linear + c(-1,1)*sqrt(discriminant)
denominator <- 2*C23$quadratic
mean.at.equal.cost <- numerator/denominator
cost.at.equal.cost <- quad(C23, mean.at.equal.cost)
equal.cost <- data.table(
  mean=mean.at.equal.cost,
  cost=cost.at.equal.cost)
ggplot()+
  scale_color_manual(values=break.colors)+
  scale_x_continuous(
    labels=sprintf("%.1f", breaks.vec),
    breaks=breaks.vec)+
  scale_fill_manual(values=c(limit=NA, equality="grey50"))+
  geom_line(aes(mean, cost, group=dt.name, color=break.after),
            data=second.choice.lines)+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(second.choice.ends, what="limit"))+
  geom_point(aes(mean, cost, fill=what),
             shape=21,
             data=data.table(equal.cost, what="equality"))
## Clear winner is break.after 3, prune the break after 1.

C24 <- AddFuns(C13.less, gamma.dt[4,])
C24$break.after <- factor(3)
C24.lines <- getLines("C24")
C24.ends <- C24.lines[mean==min.mean | mean==max.mean,]
breaks.vec <- C24.ends[, unique(c(min.mean, max.mean))]
C24$min.cost.mean <- getMinMean(C24)
min.cost.mean.dt <-
  C24[min.mean < min.cost.mean & min.cost.mean < max.mean,]
min.cost.mean.dt$min.cost <-
  quad(min.cost.mean.dt, min.cost.mean.dt$min.cost.mean)
ggplot()+
  scale_x_continuous(
    labels=sprintf("%.1f", breaks.vec),
    breaks=breaks.vec)+
  scale_color_manual(values=break.colors)+
  scale_fill_manual(values=c(limits=NA, min="black"))+
  geom_line(aes(mean, cost, color=break.after),
            data=C24.lines)+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(C24.ends, what="limits"))
quad(C24, C24$min.mean)
## Since no minimum occurs on the inside of an interval, we conclude
## that the solution for 2 segments up to data point 4 is undefined.

C33 <- gamma.dt[3,]
C33$min.mean <- 13
## Since C23 does not attain its minimum in any interval, we only need
## to consider C33.
C33.lines <- rbind(
  ##data.table(break.after=factor(3), getLines("C23")),
  data.table(break.after=factor(2), getLines("C33")))
C33.ends <- C33.lines[mean==min.mean | mean==max.mean,]
breaks.vec <- C33.ends[, unique(c(min.mean, max.mean))]
min.cost.candidates <- data.table(break.after=factor(1), C33)
min.cost.candidates$min.cost.mean <- getMinMean(C33)
min.cost.mean.dt <-
  min.cost.candidates[min.mean < min.cost.mean & min.cost.mean < max.mean,]
min.cost.mean.dt$min.cost <-
  quad(min.cost.mean.dt, min.cost.mean.dt$min.cost.mean)
ggplot()+
  ##coord_cartesian(ylim=c(0,20))+
  scale_x_continuous(breaks=breaks.vec)+
  scale_color_manual(values=break.colors)+
  geom_line(aes(mean, cost, color=break.after),
            data=C33.lines)+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(C33.ends, what="limit"))+
  ## geom_point(aes(min.cost.mean, min.cost, color=break.after, fill=what),
  ##            shape=21,
  ##            data=data.table(min.cost.mean.dt, what="min"))+
  scale_fill_manual(values=c(min="black", limit=NA))

C34 <- AddFuns(C33, gamma.dt[4,])
C34.lines <- rbind(
  ##data.table(break.after=factor(3), getLines("C23")),
  data.table(break.after=factor(2), getLines("C34")))
C34.ends <- C34.lines[mean==min.mean | mean==max.mean,]
breaks.vec <- C34.ends[, unique(c(min.mean, max.mean))]
min.cost.candidates <- data.table(break.after=factor(1), C34)
min.cost.candidates$min.cost.mean <- getMinMean(C34)
min.cost.mean.dt <-
  min.cost.candidates[min.mean < min.cost.mean & min.cost.mean < max.mean,]
min.cost.mean.dt$min.cost <-
  quad(min.cost.mean.dt, min.cost.mean.dt$min.cost.mean)
ggplot()+
  ##coord_cartesian(ylim=c(0,20))+
  scale_x_continuous(breaks=breaks.vec)+
  scale_color_manual(values=break.colors)+
  geom_line(aes(mean, cost, color=break.after),
            data=C34.lines)+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(C34.ends, what="limit"))+
  ## geom_point(aes(min.cost.mean, min.cost, color=break.after, fill=what),
  ##            shape=21,
  ##            data=data.table(min.cost.mean.dt, what="min"))+
  scale_fill_manual(values=c(min="black", limit=NA))
