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
getLessMin <- function(dt){
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
AddFuns <- function(dt1, dt2){
  mean.vec <- sort(unique(rbind(dt1, dt2)[, c(min.mean, max.mean)]))
  i1 <- 1
  i2 <- 1
  this.min <- min.mean
  new.dt.list <- list()
  while(i1 <= nrow(dt1) && i2 <= nrow(dt2)){
    row1 <- dt1[i1,]
    row2 <- dt2[i2,]
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
      max.mean=this.max)
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
  coord_cartesian(ylim=c(0,10))+
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
  coord_cartesian(ylim=c(0,20))+
  scale_x_continuous(breaks=breaks.vec)+
  scale_color_manual(values=break.colors)+
  geom_point(aes(mean, cost, group=dt.name, color=break.after),
            data=C23.ends)+
  geom_point(aes(min.cost.mean, min.cost, group=dt.name, color=break.after),
            data=min.cost.mean.dt)+
  geom_line(aes(mean, cost, group=dt.name, color=break.after),
            data=C23.lines)

C13.less <- getLessMin(C1.dt[3,])
second.choice.lines <- rbind(
  data.table(break.after=factor(3), getLines("C13.less")),
  data.table(break.after=factor(1), getLines("C23")))
second.choice.ends <- second.choice.lines[mean==min.mean | mean==max.mean,]
breaks.vec <- second.choice.ends[, unique(c(min.mean, max.mean))]
ggplot()+
  scale_color_manual(values=break.colors)+
  scale_x_continuous(
    labels=sprintf("%.1f", breaks.vec),
    breaks=breaks.vec)+
  coord_cartesian(ylim=c(0,20))+
  geom_point(aes(mean, cost, group=dt.name, color=break.after),
            data=second.choice.ends)+
  geom_line(aes(mean, cost, group=dt.name, color=break.after),
            data=second.choice.lines)
## No clear winner, so do not prune.

break.after.1 <- AddFuns(C23, gamma.dt[4,])
break.after.3 <- AddFuns(C13.less, gamma.dt[4,])
C24.lines <- rbind(
  data.table(break.after=factor(1), getLines("break.after.1")),
  data.table(break.after=factor(3), getLines("break.after.3")))
C24.ends <- C24.lines[mean==min.mean | mean==max.mean,]
breaks.vec <- C24.ends[, unique(c(min.mean, max.mean))]
ggplot()+
  scale_x_continuous(
    labels=sprintf("%.1f", breaks.vec),
    breaks=breaks.vec)+
  scale_color_manual(values=break.colors)+
  coord_cartesian(ylim=c(0,200))+
  geom_point(aes(mean, cost, group=dt.name, color=break.after),
            data=C24.ends)+
  geom_line(aes(mean, cost, group=dt.name, color=break.after),
            data=C24.lines)

getMinMean(break.after.1)
getMinMean(break.after.3)
