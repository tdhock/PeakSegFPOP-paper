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
  ## NOTE asymmetrey of dt1 and dt2 -- dt1 is used for data.i
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
      data.i=row1$data.i)
    this.min <- this.max
  }
  do.call(rbind, new.dt.list)
}

less.more.min.list <- list(
  not.strict=list(
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
        browser()
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
less.more.test.list <- list(
  ## list(fun=data.table(
  ##   quadratic=c(0, 1, 0),
  ##   linear=c(0, -2, 0),
  ##   constant=c(1, 1, 1),
  ##   min.mean=c(-2, -1, 1),
  ##   max.mean=c(-1, 1, 2)))
  list(input=data.table(
    quadratic=c(2,1),
    linear=c(-54, -28),
    constant=c(365, 196),
    min.mean=c(1,13),
    max.mean=c(13,14),
    data.i=1
  ), output=list(
    ## strict=list(
    ##   less=data.table(),
    ##   more=data.table(
    ##     quadratic=0,
    ##     linear=0,
    ##     constant=0,
    ##     min.mean=1,
    ##     max.mean=14
    ##   ),
    not.strict=list(
      ## less=data.table(
      ##   quadratic=c(2,1),
      ##   linear=c(-54, -28),
      ##   constant=c(365, 196),
      ##   min.mean=c(1,13),
      ##   max.mean=c(13,14)),
      more=data.table(
        quadratic=0,
        linear=0,
        constant=0,
        min.mean=1,
        max.mean=14,
        data.i=1
      )))
  )
)
for(test.case in less.more.test.list){
  for(min.type in names(test.case$output)){
    type.list <- test.case$output[[min.type]]
    for(min.fun.name in names(type.list)){
      expected <- type.list[[min.fun.name]]
      min.fun <- less.more.min.list[[min.type]][[min.fun.name]]
      computed <- min.fun(test.case$input)
      stopifnot(all(expected==computed))
    }
  }
}

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

all.cost.models <- list()
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
      first.data$data.i <- n.segments-1
      cost.model <- AddFuns(first.data, first.min)
      cost.models.list[[paste(n.segments, n.segments)]] <- cost.model
      for(timestep in (n.segments+1):length(data.vec)){
        prev.cost.model <- cost.models.list[[paste(n.segments-1, timestep-1)]]
        compare.cost <- min.fun(prev.cost.model)
        compare.cost$data.i <- timestep-1
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
        ## envelope, plus a new data point.
        cost.model <- AddFuns(one.env, gamma.dt[timestep,])
        cost.models.list[[paste(n.segments, timestep)]] <- cost.model
      }#for(timestep
    }#for(n.segments
    all.cost.models[[data.name]][[min.type]] <- cost.models.list
  }#for(min.type
}#for(data.name
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
    
