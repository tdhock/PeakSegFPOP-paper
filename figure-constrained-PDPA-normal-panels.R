source("packages.R")

## The first comparison is C12< versus C11<+gamma_2
data.list <- list(
  "1,10,14,13"=c(1, 10, 14, 13),
  "13,14,10,1"=c(13, 14, 10, 1))

quad <- function(dt, x){
  dt[, quadratic*x*x + linear*x + constant]
}
getLines <- function(dt){
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
  do.call(rbind, line.list)
}
getMinMean <- function(dt){
  dt[, -linear/(2*quadratic)]
}
AddFuns <- function(dt1, dt2){
  data.i <- min(dt1$data.i, dt2$data.i)
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
      data.i)
    this.min <- this.max
  }
  do.call(rbind, new.dt.list)
}

less.more.min.list <- list(not.strict=list(
  less=function(dt){
    if(1 < nrow(dt)){
      stop("TODO implement more general less equal min computation")
    }
    mu <- getMinMean(dt)
    cost <- quad(dt, mu)
    data.table(
      quadratic=c(dt$quadratic, 0),
      linear=c(dt$linear, 0),
      constant=c(dt$constant, cost),
      min.mean=c(min.mean, mu),
      max.mean=c(mu, max.mean),
      data.i=dt$data.i)[min.mean!=max.mean,]
  },
  more=function(dt){
    if(1 < nrow(dt)){
      stop("TODO implement more general more equal min computation")
    }
    mu <- getMinMean(dt)
    cost <- quad(dt, mu)
    data.table(
      quadratic=c(0, dt$quadratic),
      linear=c(0, dt$linear),
      constant=c(cost, dt$constant),
      min.mean=c(min.mean, mu),
      max.mean=c(mu, max.mean),
      data.i=dt$data.i)[min.mean!=max.mean,]
  }), strict=list(
    less=function(dt){
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
        data.i=dt$data.i)[min.mean!=max.mean,]
    },
    more=function(dt){
      if(1 < nrow(dt)){
        stop("TODO implement more general less min computation")
      }
      mu <- getMinMean(dt)
      cost <- quad(dt, mu)
      data.table(
        quadratic=0,
        linear=0,
        constant=cost,
        min.mean=dt$min.mean,
        max.mean=mu,
        data.i=dt$data.i)[min.mean!=max.mean,]
    }))

Minimize <- function(dt){
  dt$min.cost.mean <- getMinMean(dt)
  dt$min.cost <- quad(dt, dt$min.cost.mean)
  dt[is.finite(min.cost.mean) &
       min.mean < min.cost.mean & min.cost.mean < max.mean,]
}

MinEnvelope <- function(dt1, dt2){
  stopifnot(dt1[, min.mean != max.mean])
  stopifnot(dt2[, min.mean != max.mean])
  ## First we have to figure out which function starts out with a lower cost. 
  row1 <- dt1[1,]
  row2 <- dt2[1,]
  same.name.vec <- c("min.mean", "quadratic", "linear", "constant")
  same.mat <-
    row1[, same.name.vec, with=FALSE] ==
    row2[, same.name.vec, with=FALSE]
  is.row2 <- if(all(same.mat)){
    ## if they start out the same, then pick whichever one has the
    ## lower cost after the smallest of the two max.mean values.
    if(row1$max.mean < row2$max.mean){
      row.diff <- dt1[2,]-row2
      mean.at.last.equal.cost <- row1$max.mean
    }else{
      row.diff <- row1-dt2[2,]
      mean.at.last.equal.cost <- row2$max.mean
    }
    slope <- row.diff[, 2*quadratic*mean.at.last.equal.cost + linear]
    if(slope==0){
      0 < row.diff[, 2*quadratic] #hessian/2nd derivative.
    }else{
      0 < slope
    }
  }else if(same.mat[, "min.mean"]){
    ## They could also be different functions that start at the same
    ## min.mean.
    m <- row1$min.mean
    row1.cost <- quad(row1, m)
    row2.cost <- quad(row2, m)
    if(row1.cost == row2.cost){
      ## If the different functions have the same value at min.mean,
      ## then look at the derivative to determine which one to start
      ## with.
      row.diff <- row1-row2
      slope <- row.diff[, 2*min.mean*quadratic + linear]
      slope < 0
    }else{
      row2.cost < row1.cost
    }
  }else{
    row2$min.mean < row1$min.mean 
  }
  i1 <- 1
  i2 <- 1
  new.dt.list <- list()
  while(i1 <= nrow(dt1) && i2 <= nrow(dt2)){
    row1 <- dt1[i1,]
    row2 <- dt2[i2,]
    row.diff <- row1-row2
    discriminant <- row.diff[, linear^2 - 4*quadratic*constant]
    numerator <- -row.diff$linear + c(-1,1)*sqrt(discriminant)
    denominator <- 2*row.diff$quadratic
    mean.at.equal.cost <- numerator/denominator
    slope.at.equal.cost <- row.diff[, 2*mean.at.equal.cost*quadratic + linear]
    ## If the slope is negative, then row2 is minimal on the left and
    ## row1 is optimal on the right (and vice versa).
    in.row1 <-
      row1$min.mean < mean.at.equal.cost &
      mean.at.equal.cost < row1$max.mean
    in.row2 <-
      row2$min.mean < mean.at.equal.cost &
      mean.at.equal.cost < row2$max.mean
    in.both <- in.row1 & in.row2
    cross.dt <- data.table(
      mean=mean.at.equal.cost,
      slope=slope.at.equal.cost)[in.both,]
    new.dt.list[[paste(i1, i2)]] <- if(nrow(cross.dt)==0){
      ## No crossing points, so use current.
      if(is.row2)row2 else row1
    }else if(nrow(cross.dt)==1){
      r <- if(is.row2)rbind(row2,row1) else rbind(row1,row2)
      r$max.mean[1] <- r$min.mean[2] <- cross.dt$mean
      r
    }else if(nrow(cross.dt)==2){
      r <- if(is.row2)rbind(row2,row1,row2) else rbind(row1,row2,row1)
      r$max.mean[1:2] <- cross.dt$mean
      r$min.mean[2:3] <- cross.dt$mean
    }else stop("more than two crossing points")
    this.max <- if(row1$max.mean < row2$max.mean){
      row1$max.mean
    }else{
      row2$max.mean
    }
    if(row1$max.mean == this.max)i1 <- i1+1
    if(row2$max.mean == this.max)i2 <- i2+1
  }
  do.call(rbind, new.dt.list)
}

cost.lines.list <- list()
minima.list <- list()
envelope.list <- list()
for(data.name in names(data.list)){
  data.vec <- data.list[[data.name]]
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
  cost.models.list <- list()
  for(data.i in 1:nrow(C1.dt)){
    cost.models.list[[paste(1, data.i)]] <- C1.dt[data.i,]
  }
  for(min.type in names(less.more.min.list)){
    min.fun.list <- less.more.min.list[[min.type]]
    for(n.segments in 2:3){
      prev.cost.model <- cost.models.list[[paste(n.segments-1, n.segments-1)]]
      min.fun.name <- ifelse(n.segments %% 2, "more", "less")
      min.fun <- min.fun.list[[min.fun.name]]
      first.min <- min.fun(prev.cost.model)
      first.data <- gamma.dt[n.segments,]
      cost.model <- AddFuns(first.min, first.data)
      cost.models.list[[paste(n.segments, n.segments)]] <- cost.model
      for(timestep in (n.segments+1):length(data.vec)){
        cost.models.list[[paste(n.segments, timestep)]] <- cost.model
        prev.cost.model <- cost.models.list[[paste(n.segments-1, timestep-1)]]
        compare.cost <- min.fun(prev.cost.model)
        cost.minima <- Minimize(cost.model)
        compare.minima <- Minimize(compare.cost)
        cost.lines.list[[paste(data.name, min.type, n.segments, timestep)]] <-
          data.table(data.name, min.type, n.segments, timestep, rbind(
            getLines(cost.model), getLines(compare.cost)))
        one.env <- MinEnvelope(compare.cost, cost.model)
        envelope.list[[paste(data.name, min.type, n.segments, timestep)]] <-
          data.table(data.name, min.type, n.segments, timestep,
                     getLines(one.env))
        if(nrow(cost.minima)){
          minima.list[[paste(data.name, min.type, n.segments, timestep)]] <-
            data.table(data.name, min.type, n.segments, timestep, rbind(
              cost.minima, compare.minima))
        }
        ## Now that we are done with this step, we can perform the
        ## recursion by setting the new model of the cost to the min
        ## envelope that we just computed.
        cost.model <- AddFuns(one.env, gamma.dt[timestep,])
      }
    }
  }
}
cost.lines <- do.call(rbind, cost.lines.list)
cost.lines[, data.i.fac := factor(data.i)]
minima <- do.call(rbind, minima.list)
minima[, data.i.fac := factor(data.i)]
envelope <- do.call(rbind, envelope.list)
envelope[, data.i.fac := factor(data.i)]

break.colors <- c("1"="#E41A1C",
  "2"="#377EB8",
  "#4DAF4A",
  "3"="#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999")
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(data.name + min.type ~ timestep + n.segments, scales="free")+
  ##coord_cartesian(xlim=c(10,14), ylim=c(0,10))+
  scale_color_manual(values=break.colors)+
  geom_line(aes(mean, cost, group=data.i.fac),
            color="grey",
            size=2,
            data=envelope)+
  geom_line(aes(mean, cost, color=data.i.fac),
            data=cost.lines)+
  geom_point(aes(min.cost.mean, min.cost, color=data.i.fac),
             data=minima)
    
C11.less <- getLessMin(C1.dt[1,])
C12.less <- getLessMin(C1.dt[2,])

##dput(RColorBrewer::brewer.pal(Inf, "Set1"))

first.choice.lines <- rbind(
  data.table(break.after=factor(2), getLines("C12.less")),
  data.table(break.after=factor(1), getLines("C22")))
first.choice.ends <- first.choice.lines[mean==min.mean | mean==max.mean,]
min.cost.candidates <- data.table(break.after=factor(1), C22)
min.cost.candidates$min.cost.mean <- getMinMean(C22)
min.cost.mean.dt <-
  min.cost.candidates[min.mean < min.cost.mean & min.cost.mean < max.mean,]
min.cost.mean.dt$min.cost <-
  quad(min.cost.mean.dt, min.cost.mean.dt$min.cost.mean)
ggplot()+
  ggtitle("first choice")+
  xlab("mean of second segment")+
  scale_x_continuous(breaks=first.choice.ends[, unique(c(min.mean, max.mean))])+
  scale_color_manual(values=break.colors)+
  scale_fill_manual(values=c(min="black", limit=NA))+
  geom_line(aes(mean, cost, group=dt.name, color=break.after),
            data=first.choice.lines)+
  geom_point(aes(mean, cost, group=dt.name, color=break.after, fill=what),
             shape=21,
            data=data.table(first.choice.ends, what="limit"))+
  geom_point(aes(min.cost.mean, min.cost, color=break.after, fill=what),
             shape=21,
             data=data.table(min.cost.mean.dt, what="min"))
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
  ggtitle("Cost in two segments up to data point 3")+
  xlab("mean of second segment")+
  ##coord_cartesian(ylim=c(0,20))+
  scale_x_continuous(breaks=breaks.vec)+
  scale_color_manual(values=break.colors)+
  geom_line(aes(mean, cost, color=break.after),
            data=C23.lines)+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(C23.ends, what="limit"))+
  geom_point(aes(min.cost.mean, min.cost, color=break.after, fill=what),
             shape=21,
             data=data.table(min.cost.mean.dt, what="min"))+
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
  ggtitle("second choice")+
  scale_color_manual(values=break.colors)+
  scale_x_continuous(
    "mean of second segment",
    labels=sprintf("%.1f", breaks.vec),
    breaks=breaks.vec)+
  scale_fill_manual(values=c(limit=NA, equality="grey50"))+
  geom_line(aes(mean, cost, group=dt.name, color=break.after),
            data=second.choice.lines)+
  ## geom_point(aes(mean, cost, fill=what),
  ##            shape=21,
  ##            data=data.table(equal.cost, what="equality"))+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(second.choice.ends, what="limit"))
## Clear winner is break.after 1, prune the break after 3.

C24 <- AddFuns(C23, gamma.dt[4,])
C24$break.after <- factor(1)
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
  geom_point(aes(min.cost.mean, min.cost, color=break.after, fill=what),
             shape=21,
             data=data.table(min.cost.mean.dt, what="min"))+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(C24.ends, what="limits"))
## This curve achieves its minimum within the interval, so we conclude
## that we have found the minimum error model with two segments.

C33 <- gamma.dt[3,]
C33$max.mean <- 10
C23more <- getMoreMin(C23)
C33.lines <- rbind(
  data.table(break.after=factor(3), getLines("C23more")),
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
  ggtitle("last choice")+
  xlab("mean of third segment")+
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
## We can prune the break after 2 now. C23more is always less costly
## than C33.

C34 <- AddFuns(C23more, gamma.dt[4,])
C34.lines <- rbind(
  ##data.table(break.after=factor(3), getLines("C23more")),
  data.table(break.after=factor(3), getLines("C34")))
C34.ends <- C34.lines[mean==min.mean | mean==max.mean,]
breaks.vec <- C34.ends[, unique(c(min.mean, max.mean))]
min.cost.candidates <- data.table(break.after=factor(1), C34)
min.cost.candidates$min.cost.mean <- getMinMean(C34)
min.cost.candidates$min.cost <-
  quad(min.cost.candidates, min.cost.candidates$min.cost.mean)
ggplot()+
  ##coord_cartesian(ylim=c(0,20))+
  scale_x_continuous(breaks=breaks.vec)+
  scale_color_manual(values=break.colors)+
  geom_line(aes(mean, cost, color=break.after),
            data=C34.lines)+
  geom_point(aes(mean, cost, color=break.after, fill=what),
             shape=21,
             data=data.table(C34.ends, what="limit"))+
  geom_point(aes(min.cost.mean, min.cost, color=break.after, fill=what),
             shape=21,
             data=data.table(min.cost.candidates, what="min"))+
  scale_fill_manual(values=c(min="black", limit=NA))
## Since the minimum occurs outside of the interval, we can conclude
## that the optimal model is undefined.
