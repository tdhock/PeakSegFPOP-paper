source("packages.R")
data(intreg, package="animint")

options(warn=2)

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
  if(nrow(dt1)==0)return(dt1)
  if(nrow(dt2)==0)return(dt2)
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
      min.mean=c(dt$min.mean[1], mu),
      max.mean=c(mu, dt$max.mean[nrow(dt)]),
      data.i=dt$data.i)[min.mean!=max.mean,]
  },
  more=function(dt){
    mean.at.min.vec <- getMinMean(dt)
    min.achieved <-
      dt[, min.mean < mean.at.min.vec & mean.at.min.vec < max.mean]
    stopifnot(nrow(min.achieved)==1)
    new.dt.list <- list()
    overall.min.cost <- NULL
    row.i <- nrow(dt)
    while(is.null(overall.min.cost)){
      r <- dt[row.i,]
      mean.at.min <- getMinMean(r)
      cost.at.min <- quad(r, mean.at.min)
      if(r$max.mean <= mean.at.min){
        min.position <- "right"
        decreasing.to.right <- TRUE
        overall.min.cost <- quad(r, r$max.mean)
        new.dt.list[[paste(row.i)]] <- data.table(
          quadratic=0,
          linear=0,
          constant=overall.min.cost,
          min.mean=dt$min.mean[1],
          max.mean=r$max.mean,
          data.i=r$data.i)
      }else if(mean.at.min <= r$min.mean){
        min.position <- "left"
        decreasing.to.right <- FALSE
        new.dt.list[[paste(row.i)]] <- r
      }else{
        min.position <- "inside"
        overall.min.cost <- cost.at.min
        decreasing.to.right <- FALSE
        new.dt.list[[paste(row.i)]] <- data.table(
          quadratic=c(0, r$quadratic),
          linear=c(0, r$linear),
          constant=c(overall.min.cost, r$constant),
          min.mean=c(dt$min.mean[1], mean.at.min),
          max.mean=c(mean.at.min, r$max.mean),
          data.i=r$data.i)
      }
      min.text <- data.table(
        mean.at.min, cost.at.min, min.position,
        decreasing.to.right)
      gg.funs <- ggplot()+
        geom_line(aes(mean, cost),
                  data=getLines(dt),
                  size=2,
                  color="grey")+
        geom_line(aes(mean, cost),
                  data=getLines(r),
                  size=1,
                  color="black")+
        geom_point(aes(mean.at.min, cost.at.min),
                   data=min.text,
                   shape=21,
                   fill="white")+
        geom_text(aes(mean.at.min, cost.at.min, label=paste(
          min.position,
          ifelse(decreasing.to.right, "decreasing", "increasing"))),
          vjust=1.5,
          data=min.text)
      ##print(gg.funs)
      row.i <- row.i-1
    }
    do.call(rbind, rev(new.dt.list))
  })
less.more.test.list <- list(
  real=list(input=data.table(
    quadratic=c(2,1),
    linear=c(-54, -28),
    constant=c(365, 196),
    min.mean=c(1,13),
    max.mean=c(13,14),
    data.i=1
  ), output=list(
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
      ))
  )
)
for(test.case.name in names(less.more.test.list)){
  test.case <- less.more.test.list[[test.case.name]]
  input <- test.case$input
  for(min.fun.name in names(test.case$output)){
    expected <- test.case$output[[min.fun.name]]
    min.fun <- less.more.min.list[[min.fun.name]]
    computed <- min.fun(input)
    test.lines.list <- list()
    for(object.name in c("input","expected","computed")){
      dt <- get(object.name)
      if(nrow(dt)){
        test.lines.list[[object.name]] <- data.table(
          object.name=factor(object.name, c("input", "expected", "computed")),
          getLines(dt))
      }
    }
    test.lines <- do.call(rbind, test.lines.list)
    gg.test <- ggplot()+
      ggtitle(paste(test.case.name, min.fun.name))+
      scale_size_manual(values=c(input=3, expected=2, computed=1))+
      scale_color_manual(values=c(
        input="grey", expected="black", computed="red"))+
      geom_line(aes(mean, cost, color=object.name, size=object.name),
                data=test.lines)
    if(!identical(expected, computed)){
      print(gg.test)
      print(list(input=input, expected=expected, computed=computed))
      stop("expected not identical to computed")
    }
  }
}

Minimize <- function(dt){
  dt$min.cost.mean <- getMinMean(dt)
  dt$min.cost <- quad(dt, dt$min.cost.mean)
  dt[is.finite(min.cost.mean) &
       min.mean <= min.cost.mean & min.cost.mean <= max.mean,]
}

MinEnvelope <- function(dt1, dt2){
  stopifnot(dt1[, min.mean < max.mean])
  stopifnot(dt2[, min.mean < max.mean])
  if(nrow(dt1)==0)return(dt2)
  if(nrow(dt2)==0)return(dt1)
  ## First we have to figure out which function starts out with a lower cost. 
  row1 <- dt1[1,]
  row2 <- dt2[1,]
  same.name.vec <- c("min.mean", "quadratic", "linear", "constant")
  same.mat <-
    row1[, same.name.vec, with=FALSE] ==
    row2[, same.name.vec, with=FALSE]
  ggplot()+
    geom_line(aes(mean, cost, color=fun.i),
              size=2,
              data=data.table(getLines(dt1), fun.i=factor(1)))+
    geom_line(aes(mean, cost, color=fun.i),
              size=2,
              data=data.table(getLines(dt2), fun.i=factor(2)))+
    geom_line(aes(mean, cost),
              data=getLines(row1))+
    geom_line(aes(mean, cost),
              data=getLines(row2))
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
  this.row <- if(is.row2)row2 else row1
  other.row <- if(is.row2)row1 else row2
  last.min.mean <- this.row$min.mean
  i1 <- 1
  i2 <- 1
  new.dt.list <- list()
  while(i1 <= nrow(dt1) && i2 <= nrow(dt2)){
    row1 <- dt1[i1,]
    row2 <- dt2[i2,]
    this.row <- if(is.row2)row2 else row1
    other.row <- if(is.row2)row1 else row2
    row.diff <- row1-row2
    row.diff$min.mean <- min(row1$min.mean, row2$min.mean)
    row.diff$max.mean <- max(row1$max.mean, row2$max.mean)
    ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(y ~ ., scales="free")+
      geom_line(aes(mean, cost, color=fun.i),
                data.table(getLines(row1),fun.i=factor(1),y="cost"))+
      geom_line(aes(mean, cost, color=fun.i),
                data.table(getLines(row2),fun.i=factor(2),y="cost"))+
      geom_line(aes(mean, cost),
                data.table(getLines(row.diff),y="diff"))
    discriminant <- row.diff[, linear^2 - 4*quadratic*constant]
    cross.dt <- if(0 < discriminant){
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
      data.table(
        mean=mean.at.equal.cost,
        slope=slope.at.equal.cost)[in.both,]
    }else{
      ## one or no intersection points, so one function is always above
      ## the other.
      data.table()
    }
    if(nrow(cross.dt)==1){
      before.cross <- this.row
      before.cross$min.mean <- last.min.mean
      before.cross$max.mean <- cross.dt$mean
      last.min.mean <- cross.dt$mean
      new.dt.list[[paste(i1, i2, "cross")]] <- before.cross
      is.row2 <- !is.row2
      tmp.row <- this.row
      this.row <- other.row
      other.row <- tmp.row
    }else if(nrow(cross.dt)==2){
      r <- rbind(this.row, other.row)
      r$max.mean <- cross.dt$mean
      r$min.mean <- c(last.min.mean, cross.dt$mean[1])
      new.dt.list[[paste(i1, i2, "cross")]] <- r
      last.min.mean <- cross.dt$mean[2]
    }
    if(this.row$max.mean < other.row$max.mean){
      ## this function piece stops but the other function continues to
      ## a larger mean, so store it and move on to the next piece of
      ## this function.
      r <- this.row
      r$min.mean <- last.min.mean
      new.dt.list[[paste(i1,i2)]] <- r
      last.min.mean <- this.row$max.mean
      if(is.row2){
        i2 <- i2+1
      }else{
        i1 <- i1+1
      }
    }else{
      ## this function piece continues to a larger mean.
      if(is.row2){
        i1 <- i1+1
      }else{
        i2 <- i2+1
      }
    }
  }
  this.row$min.mean <- last.min.mean
  new.dt.list[["last"]] <- this.row
  new.dt <- do.call(rbind, new.dt.list)
  browser(expr=any(table(new.dt$min.mean)>1))
  browser(expr=any(table(new.dt$max.mean)>1))
  new.dt
}

all.cost.models <- list()
cost.lines.list <- list()
minima.list <- list()
envelope.list <- list()
data.vec <- -subset(intreg$signals, signal=="4.2")$logratio[80:200]
min.mean <- min(data.vec)
max.mean <- max(data.vec)
gamma.dt <- data.table(
  quadratic=1,
  linear=-2*data.vec,
  constant=data.vec * data.vec)
C1.dt <- cumsum(gamma.dt)
gamma.dt$min.mean <- C1.dt$min.mean <- min.mean
gamma.dt$max.mean <- C1.dt$max.mean <- max.mean
gamma.dt$data.i <- C1.dt$data.i <- 0
cost.models.list <- list()
for(data.i in 1:nrow(C1.dt)){
  cost.models.list[[paste(1, data.i)]] <- C1.dt[data.i,]
}
max.segments <- 3
for(n.segments in 2:max.segments){
  prev.cost.model <- cost.models.list[[paste(n.segments-1, n.segments-1)]]
  min.fun.name <- ifelse(n.segments %% 2, "more", "less")
  min.fun <- less.more.min.list[[min.fun.name]]
  first.min <- min.fun(prev.cost.model)
  first.data <- gamma.dt[n.segments,]
  first.data$data.i <- n.segments-1
  cost.model <- AddFuns(first.data, first.min)
  cost.models.list[[paste(n.segments, n.segments)]] <- cost.model
  for(timestep in (n.segments+1):length(data.vec)){
    cat(sprintf("%4d / %4d segments %4d / %4d data points %d intervals\n",
                n.segments, max.segments, timestep, length(data.vec),
                nrow(cost.model)))
    prev.cost.model <- cost.models.list[[paste(n.segments-1, timestep-1)]]
    compare.cost <- min.fun(prev.cost.model)
    compare.cost$data.i <- timestep-1
    cost.minima <- Minimize(cost.model)
    compare.minima <- Minimize(compare.cost)
    if(nrow(compare.cost)){
      cost.lines.list[[
        paste(n.segments, timestep, "compare")]] <-
        data.table(n.segments, timestep,
                   compare.cost.lines <- getLines(compare.cost))
    }
    if(nrow(cost.model)){ # may be Inf over entire interval.
      cost.lines.list[[paste(n.segments, timestep)]] <-
        data.table(n.segments, timestep,
                   cost.model.lines <- getLines(cost.model))
    }
    one.env <- MinEnvelope(compare.cost, cost.model)
    if(nrow(one.env)){
      envelope.list[[paste(n.segments, timestep)]] <-
        data.table(n.segments, timestep,
                   env.lines <- getLines(one.env))
    }
    if(nrow(cost.minima)){
      minima.list[[paste(n.segments, timestep)]] <-
        data.table(n.segments, timestep, rbind(
          cost.minima, compare.minima))
    }
    ## Now that we are done with this step, we can perform the
    ## recursion by setting the new model of the cost to the min
    ## envelope, plus a new data point.
    ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      geom_line(aes(mean, cost),
                color="grey",
                size=2,
                data=env.lines)+
      geom_line(aes(mean, cost, color=factor(data.i)),
                data=cost.model.lines)+
      geom_line(aes(mean, cost, color=factor(data.i)),
                data=compare.cost.lines)
    cost.model <- AddFuns(one.env, gamma.dt[timestep,])
    cost.models.list[[paste(n.segments, timestep)]] <- cost.model
  }#for(timestep
}#for(n.segments
cost.lines <- do.call(rbind, cost.lines.list)
cost.lines[, data.i.fac := factor(data.i)]
minima <- do.call(rbind, minima.list)
minima[, data.i.fac := factor(data.i)]
envelope <- do.call(rbind, envelope.list)
envelope[, data.i.fac := factor(data.i)]

gg.pruning <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(min.type ~ timestep + n.segments, scales="free",
             labeller=function(var, val){
               if(var %in% c("n.segments", "timestep")){
                 paste(var, "=", val)
               }else{
                 paste(val)
               }
             })+
  geom_line(aes(mean, cost, group=data.i.fac),
            color="grey",
            size=2,
            data=envelope)+
  geom_line(aes(mean, cost, color=data.i.fac),
            data=cost.lines)+
  geom_point(aes(min.cost.mean, min.cost, color=data.i.fac),
             data=minima)
pdf("figure-constrained-PDPA-normal-panels-pruning.pdf")
print(gg.pruning)
dev.off()

data.lines.list <- list()
data.minima.list <- list()
data.infeasible.list <- list()
data.cost.list <- list()
data.intervals.list <- list()
## for animint:
for(total.segments in 1:max.segments){
  for(timestep in total.segments:length(data.vec)){
    cat(sprintf("decoding %4d / %4d segments %4d / %4d data points\n",
                total.segments, max.segments, timestep, length(data.vec)))
    data.i <- timestep
    seg.i <- total.segments
    no.constraint <- data.table(
      quadratic=0,
      linear=0,
      constant=0,
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
      min.dt <- Minimize(unconstrained.fun)
      min.dt$segment.end <- segment.end
      min.dt[, segment.start := ifelse(seg.i==1, 1, 1+data.i)]
      segment.end <- min.dt$data.i
      if(min.mean < constraint$min.mean){
        data.infeasible.list[[paste(total.segments, timestep, seg.i)]] <-
          data.table(total.segments, timestep, seg.i,
                     min.mean, max.mean=constraint$min.mean)
      }
      if(constraint$max.mean < max.mean){
        data.infeasible.list[[paste(total.segments, timestep, seg.i)]] <-
          data.table(total.segments, timestep, seg.i,
                     min.mean=constraint$max.mean,
                     max.mean)
      }
      optimal.cost <- if(nrow(min.dt)==0){
        Inf
      }else{
        constraint.status <- "inactive"
        if(min.dt$min.cost.mean < constraint$min.mean){
          min.dt$min.cost.mean <- constraint$min.mean
          constraint.status <- "active"
        }
        if(constraint$max.mean < min.dt$min.cost.mean){
          min.dt$min.cost.mean <- constraint$max.mean
          constraint.status <- "active"
        }
        constrained.fun <- unconstrained.fun[
          min.mean <= min.dt$min.cost.mean &
            min.dt$min.cost.mean <= max.mean,]
        min.dt$min.cost <- quad(constrained.fun, min.dt$min.cost.mean)
        data.cost.list[[
          paste(total.segments, timestep)]]$constraint <- constraint.status
        min.dt$constraint <- constraint.status
        show.min <- data.table(total.segments, timestep, seg.i,
          min.dt)
        ##if(seg.i>1)show.min$data.i <- show.min$data.i+1L
        if(seg.i==1)show.min$data.i <- 0
        data.minima.list[[paste(total.segments, timestep, seg.i)]] <-
          show.min
        constraint <- no.constraint
        constraint.side <- if(seg.i %% 2){
          constraint$min.mean <- min.dt$min.cost.mean
        }else{
          constraint$max.mean <- min.dt$min.cost.mean
        }
        min.dt$min.cost
      }
      if(seg.i==total.segments){
        data.cost.list[[paste(total.segments, timestep)]] <-
          data.table(total.segments, timestep, optimal.cost,
                     constraint="inactive")
      }
      data.i <- min.dt$data.i
      seg.i <- seg.i-1
    }#while(...
  }#for(timestep
}#for(total.segments

data.intervals <- do.call(rbind, data.intervals.list)
data.dt <- data.table(
  count=data.vec,
  position=seq_along(data.vec),
  data.i.fac=factor(seq_along(data.vec)))
data.lines <- do.call(rbind, data.lines.list)
data.lines[, data.i.fac := factor(data.i)]
data.lines[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
limits <- data.lines[, .SD[1,],
  by=.(total.segments, timestep, minimization, seg.i, piece.i)]
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
  data.table(dt, y=factor(y, c("count", "intervals", "segments")))
}
viz <- list(
  title=paste(
    "Constrained Pruned Dynamic Programming Algorithm"),
  funModels=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=500, height=300)+
    facet_grid(. ~ seg.i, scales="free", labeller=function(var, val){
      paste("segment", val)
    })+
    geom_tallrect(aes(xmin=min.mean, xmax=max.mean,
                      showSelected=minimization),
                  fill="grey",
                  alpha=0.5,
                  color=NA,
                  data=data.infeasible)+
    geom_line(aes(mean, cost, color=data.i.fac,
                  group=piece.i,
                  showSelected=minimization),
              data=data.lines)+
    geom_point(aes(min.cost.mean, min.cost,
                   showSelected=minimization),
               size=5,
               shape=21,
               fill="white",
               color="black",
               data=data.minima[constraint=="active",])+
    guides(color="none")+
    geom_point(aes(min.cost.mean, min.cost, color=data.i.fac,
                   tooltip=paste(
                     "minimum cost =",
                     round(min.cost, 2),
                     "at mean =",
                     round(min.cost.mean, 2),
                     "for",
                     seg.i,
                     "segment model up to data point",
                     timestep
                   ),
                   showSelected=minimization),
               size=4,
               data=data.minima)+
    geom_point(aes(mean, cost, 
                  showSelected=minimization),
              data=limits),
  data=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=800, height=300)+
    facet_grid(y ~ ., scales="free")+
    geom_point(aes(position, count),
               data=addY(data.dt, "count"))+
    geom_segment(aes(segment.start-0.5, min.cost.mean,
                     showSelected=minimization,
                     xend=segment.end+0.5, yend=min.cost.mean),
                 data=addY(data.minima, "count"),
                 color="green")+
    guides(color="none")+
    geom_point(aes(timestep, intervals),
               data=addY(data.intervals[total.segments==2,],"intervals"))+
    geom_tile(aes(timestep, total.segments, fill=optimal.cost,
                  clickSelects=minimization),
              data=addY(data.cost, "segments"))+
    geom_point(aes(timestep, total.segments, 
                   clickSelects=minimization),
               shape=21,
               color="black",
               fill="white",
               data=addY(data.cost[constraint=="active",], "segments"))
)
animint2dir(viz, "figure-constrained-PDPA-normal-real")

intervalsPlot <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(total.segments ~ .)+
  geom_point(aes(timestep, intervals),
             data=data.intervals)

pdf("figure-constrained-PDPA-normal-real.pdf")
print(intervalsPlot)
dev.off()
