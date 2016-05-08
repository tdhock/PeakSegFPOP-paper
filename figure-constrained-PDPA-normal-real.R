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

Minimize <- function(dt, from=min(dt$min.mean), to=max(dt$max.mean)){
  stopifnot(from < to)
  is.before <- dt$max.mean < from
  is.after <- to < dt$min.mean
  feasible <- dt[!(is.before | is.after),]
  feasible$min.mean[1] <- from
  feasible$max.mean[nrow(feasible)] <- to
  feasible$quad.min.mean <- getMinMean(feasible)
  feasible[, min.cost.mean := ifelse(
    quad.min.mean < min.mean, min.mean,
    ifelse(max.mean < quad.min.mean,
           max.mean,
           quad.min.mean))]
  feasible[, min.cost := quad(feasible, min.cost.mean)]
  feasible[which.min(min.cost),]
}

CompareRows <- function(row1, row2, insignificant.cost.diff){
  last.min.mean <- if(row1$min.mean < row2$min.mean){
    row2$min.mean
  }else{
    row1$min.mean
  }
  first.max.mean <- if(row1$max.mean < row2$max.mean){
    row1$max.mean
  }else{
    row2$max.mean
  }
  stopifnot(last.min.mean < first.max.mean)
  if(row1$quadratic==row2$quadratic){
    return(data.table(
      mean=c(first.min.mean, last.max.mean),
      status="equal"))
  }
  first.max.points <- 
    data.table(mean=first.max.mean, cost=c(cost1, cost2))
  ggplot()+
    ##coord_cartesian(xlim=c(-0.05, -0.1), ylim=c(0,1e-3))+
    geom_line(aes(mean, cost, color=fun.i),
              size=2,
              data=data.table(getLines(dt1), fun.i=factor(1)))+
    geom_line(aes(mean, cost, color=fun.i),
              size=1,
              data=data.table(getLines(dt2), fun.i=factor(2)))+
    geom_line(aes(mean, cost),
              linetype="dashed",
              data=getLines(row1))+
    geom_line(aes(mean, cost),
              linetype="dotted",
              size=1,
              data=getLines(row2))+
    geom_point(aes(mean, cost),
               shape=1,
               data=first.max.points)
}
MinEnvelope <- function(dt1, dt2){
  stopifnot(dt1[, min.mean < max.mean])
  stopifnot(dt2[, min.mean < max.mean])
  stopifnot(dt1[1, min.mean]==dt2[1, min.mean])
  stopifnot(dt1[.N, max.mean]==dt2[.N, max.mean])
  insignificant.list <- list(
    ##machine=.Machine$double.eps
  )
  for(dt.i in 1:2){
    if(dt.i==1){
      dt <- dt1
      other.dt <- dt2
    }else{
      dt <- dt2
      other.dt <- dt1
    }
    row.i <- 1
    while(row.i < nrow(dt)){
      this.row <- dt[row.i,]
      next.row <- dt[row.i+1,]
      mean.value <- this.row$max.mean
      this.cost <- quad(this.row, mean.value)
      next.cost <- quad(next.row, mean.value)
      insignificant.cost.diff <- abs(this.cost-next.cost)
      row.i <- row.i+1
    }
  }
  ## First we have to figure out which function starts out with a
  ## lower cost. To do that we scan from the left until we find a
  ## difference.
  i1 <- 1
  i2 <- 1
  is.row2 <- NULL
  while(is.null(is.row2)){
    if(nrow(dt1) < i1 || nrow(dt2) < i2){
      ## Both are the same over the entire region, so just pick the
      ## one with the fewer number of intervals.
      is.row2 <- nrow(dt2) < nrow(dt1)
    }else{
      row1 <- dt1[i1,]
      row2 <- dt2[i2,]
      first.max.mean <- if(row1$max.mean < row2$max.mean){
        row1$max.mean
      }else{
        row2$max.mean
      }
      cost1 <- quad(row1, first.max.mean)
      cost2 <- quad(row2, first.max.mean)
      if(row1$quadratic==row2$quadratic){
        if(row1$max.mean==first.max.mean){
          i1 <- i1+1
        }
        if(row2$max.mean==first.max.mean){
          i2 <- i2+1
        }
      }else{
        is.row2 <- cost2 < cost1
      }   
      first.max.points <- 
        data.table(mean=first.max.mean, cost=c(cost1, cost2))
      ggplot()+
        ##coord_cartesian(xlim=c(-0.05, -0.1), ylim=c(0,1e-3))+
        geom_line(aes(mean, cost, color=fun.i),
                  size=2,
                  data=data.table(getLines(dt1), fun.i=factor(1)))+
        geom_line(aes(mean, cost, color=fun.i),
                  size=1,
                  data=data.table(getLines(dt2), fun.i=factor(2)))+
        geom_line(aes(mean, cost),
                  linetype="dashed",
                  data=getLines(row1))+
        geom_line(aes(mean, cost),
                  linetype="dotted",
                  size=1,
                  data=getLines(row2))+
        geom_point(aes(mean, cost),
                   shape=1,
                   data=first.max.points)
    }
  }
  ## At this point we know that if(is.row2) then dt2 has a lower cost
  ## at least until the first max.mean in i1,i2.
  new.dt.list <- list()
  last.min.mean <- dt1[1, min.mean]
  add.i <- 1
  if(is.row2){
    while(add.i < i2){
      new.dt.list[[paste(add.i)]] <- dt2[add.i,]
      last.min.mean <- dt2[add.i, max.mean]
      add.i <- add.i+1
    }
    if(row2$max.mean < row1$max.mean){
      last.min.mean <- row2$max.mean
      new.dt.list[[paste(i2)]] <- row2
      i2 <- i2+1
    }
  }else{
    while(add.i < i1){
      new.dt.list[[paste(add.i)]] <- dt1[add.i,]
      last.min.mean <- dt1[add.i, max.mean]
      add.i <- add.i+1
    }
    if(row1$max.mean < row2$max.mean){
      last.min.mean <- row1$max.mean
      new.dt.list[[paste(i1)]] <- row1
      i1 <- i1+1
    }
  }
  ## Now go through the rest of the pieces looking for crossing
  ## points.
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
    cross.dt <- if(.Machine$double.eps < discriminant){
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
    if(this.row$max.mean == other.row$max.mean){
      r <- this.row
      r$min.mean <- last.min.mean
      new.dt.list[[paste(i1,i2)]] <- r
      last.min.mean <- this.row$max.mean
      i1 <- i1+1
      i2 <- i2+1
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
      ## this function piece continues to a larger mean, and the other
      ## function stops before.
      if(is.row2){
        i1 <- i1+1
      }else{
        i2 <- i2+1
      }
    }
  }
  this.row$min.mean <- last.min.mean
  new.dt <- do.call(rbind, new.dt.list)
  browser(expr=any(table(new.dt$min.mean)>1))
  browser(expr=any(table(new.dt$max.mean)>1))
  new.dt
}
## A test case for MinEnvelope:
result <- data.table(
  quadratic = c(3, 2, 1),
  linear = c(0.0757353883540654, -0.0842153650871264, -0.132854723477952),
  constant = c(0.0114001020653698, 0.00507885284652773, 0.00703123698800913),
  min.mean = c(-0.0922074380970889, -0.0713260031733877, 0.026116725688849),
  max.mean = c(-0.0713260031733877, 0.026116725688849, 0.0954195650786824),
  data.i = c(1, 2, 3))
## list format is: input1, input2, output.
min.env.test.list <- list(not.quasiconvex=list(
  data.table(
    quadratic = c(4, 0),
    linear = c(0.260150264548243, 0),
    constant = c(0.0199023137057983, 0.0156724286967657),
    min.mean = c(-0.0922074380970889, -0.0325187830685304),
    max.mean = c(-0.0325187830685304, 0.0954195650786824),
    data.i = c(4, 4)),
  result,
  result),
  equal.at.max=list(
    data.table(
      quadratic = c(1, 7, 0),
      linear = c(-0.02585234888854, 0.580695481248607, 0),
      constant = c(0.0664360533749589, 0.0548316584743037, 0.0427885426906412),
      min.mean = c(-0.15521264992094, -0.117545121200233, -0.0414782486606148),
      max.mean = c(-0.117545121200233, -0.0414782486606148, 0.0922074380970889),
      data.i = c(10, 10, 10)),
    data.table(
      quadratic = c(7, 4, 2),
      linear = c(0.580695481248607, -0.0534236722285902, -0.155618052058249),
      constant = c(0.0548316584743037, 0.0213227334070839, 0.0277237775077624),
      min.mean = c(-0.15521264992094, -0.105686525579533, 0.0365259791823832),
      max.mean = c(-0.105686525579533, 0.0365259791823832, 0.0922074380970889),
      data.i = c(6, 6, 8)),
    data.table(
      quadratic = c(1,
        7, 4, 2),
      linear = c(-0.02585234888854,
        0.580695481248607, -0.0534236722285902, -0.155618052058249),
      constant = c(0.0664360533749589,
        0.0548316584743037, 0.0213227334070839, 0.0277237775077624),
      min.mean = c(-0.15521264992094,
        -0.117545121200233, -0.105686525579533, 0.0365259791823832),
      max.mean = c(-0.117545121200233,
        -0.105686525579533, 0.0365259791823832, 0.0922074380970889),
      data.i = c(10,
        6, 6, 8))))
for(min.env.test.name in names(min.env.test.list)){
  test <- min.env.test.list[[min.env.test.name]]
  dt1 <- test[[1]]
  dt2 <- test[[2]]
  computed <- MinEnvelope(dt1, dt2)
  expected <- test[[3]]
  gg.test <- ggplot()+
    ggtitle(min.env.test.name)+
    geom_line(aes(mean, cost),
              color="grey",
              size=4,
              data=getLines(expected))+
    geom_line(aes(mean, cost, color=fun.i),
              size=2.5,
              data=data.table(getLines(dt1), fun.i=factor(1)))+
    geom_line(aes(mean, cost, color=fun.i),
              size=1.5,
              data=data.table(getLines(dt2), fun.i=factor(2)))+
    geom_line(aes(mean, cost),
              size=0.5,
              data=getLines(computed))
  if(!identical(expected, computed)){
    print(gg.test)
    print(list(expected=expected, computed=computed))
    stop("expected not identical to computed")
  }
}

all.cost.models <- list()
cost.lines.list <- list()
minima.list <- list()
envelope.list <- list()
data.vec <- -subset(intreg$signals, signal=="4.2")$logratio[80:200]
## TODO: increase the number of data points and see where the bug is
## coming from.
data.vec <- data.vec[1:11]
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
  if(n.segments %% 2){
    min.fun.name <- "more"
    env.transform <- function(dt){
      new.max <- -dt$min.mean
      new.min <- -dt$max.mean
      dt$min.mean <- new.min
      dt$max.mean <- new.max
      dt[, linear := -linear]
      dt[.N:1,]
    }
  }else{
    min.fun.name <- "less"
    env.transform <- identity
  }
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
    cost.model <- cost.models.list[[paste(n.segments, timestep-1)]]
    cost.minima <- Minimize(cost.model)
    compare.minima <- Minimize(compare.cost)
    gg <- ggplot()+
      ggtitle(paste(n.segments, "segments,", timestep, "data points"))
    if(nrow(compare.cost)){
      cost.lines.list[[
        paste(n.segments, timestep, "compare")]] <-
        data.table(n.segments, timestep,
                   compare.cost.lines <- getLines(compare.cost))
      gg <- gg+
        geom_line(aes(mean, cost, color=factor(data.i)),
                  compare.cost.lines)
    }
    if(nrow(cost.model)){ # may be Inf over entire interval.
      cost.lines.list[[paste(n.segments, timestep)]] <-
        data.table(n.segments, timestep,
                   cost.model.lines <- getLines(cost.model))
      gg <- gg+
        geom_line(aes(mean, cost, color=factor(data.i)),
                  cost.model.lines)
    }
    transformed.compare <- env.transform(compare.cost)
    transformed.model <- env.transform(cost.model)
    transformed.env <-
      MinEnvelope(transformed.compare, transformed.model)
    one.env <- env.transform(transformed.env)
    if(nrow(one.env)){
      envelope.list[[paste(n.segments, timestep)]] <-
        data.table(n.segments, timestep,
                   env.lines <- getLines(one.env))
      gg <- gg+
        geom_line(aes(mean, cost),
                  env.lines,
                  color="grey",
                  size=2,
                  alpha=0.5)
    }
    if(nrow(cost.minima)){
      minima.list[[paste(n.segments, timestep)]] <-
        data.table(n.segments, timestep, rbind(
          cost.minima, compare.minima))
    }
    ## Now that we are done with this step, we can perform the
    ## recursion by setting the new model of the cost to the min
    ## envelope, plus a new data point.
    cost.model <- AddFuns(one.env, gamma.dt[timestep,])
    cost.models.list[[paste(n.segments, timestep)]] <- cost.model
  }#for(timestep
}#for(n.segments
cost.lines <- do.call(rbind, cost.lines.list)
cost.lines[, minimization := paste(
  n.segments, "segments up to data point", timestep)]
cost.lines[, data.i.fac := factor(data.i)]
minima <- do.call(rbind, minima.list)
minima[, minimization := paste(
  n.segments, "segments up to data point", timestep)]
minima[, data.i.fac := factor(data.i)]
envelope <- do.call(rbind, envelope.list)
envelope[, data.i.fac := factor(data.i)]
envelope[, minimization := paste(
  n.segments, "segments up to data point", timestep)]

gg.pruning <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(timestep ~ n.segments, scales="free",
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

ti <- 11
gg.pruning <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(timestep ~ n.segments, scales="free",
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
            data=envelope[timestep==ti,])+
  geom_line(aes(mean, cost, color=data.i.fac),
            data=cost.lines[timestep==ti,])+
  geom_point(aes(min.cost.mean, min.cost, color=data.i.fac),
             data=minima[timestep==ti,])
print(gg.pruning)

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
      constraint.status <- if(min.dt[, quad.min.mean == min.cost.mean]){
        "inactive"
      }else{
        if(constraint$min.mean == min.dt$min.cost.mean){
          data.infeasible.list[[paste(
            total.segments, timestep, seg.i)]] <-
            data.table(total.segments, timestep, seg.i,
                       min.mean,
                       max.mean=constraint$min.mean)
        }else{
          data.infeasible.list[[paste(
            total.segments, timestep, seg.i)]] <-
            data.table(total.segments, timestep, seg.i,
                       min.mean=constraint$max.mean,
                       max.mean)
        }
        "active"
      }
      if(constraint.status=="active"){
        data.cost.list[[paste(
          total.segments, timestep)]]$constraint <- constraint.status
      }
      min.dt$constraint <- constraint.status
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
  count=data.vec,
  position=seq_along(data.vec),
  data.i.fac=factor(seq_along(data.vec)))
data.lines <- do.call(rbind, data.lines.list)
data.lines[, data.i.fac := factor(data.i)]
data.lines[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
left.of.intervals <-
  data.lines[, .SD[1,],
             by=.(total.segments, timestep, minimization, seg.i, piece.i)]
between.intervals <- left.of.intervals[min.mean != min(data.vec),]
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
largest.constant <- envelope[quadratic==0, max(constant)]
viz <- list(
  title=paste(
    "Constrained Pruned Dynamic Programming Algorithm"),
  funModels=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=500, height=300)+
    coord_cartesian(ylim=c(0, max(between.intervals$cost)))+
    geom_line(aes(mean, cost,
                  showSelected=minimization),
              color="grey",
              size=8,
              data=data.table(envelope, seg.i="pruning"))+
    geom_line(aes(mean, cost, color=data.i.fac,
                  group=paste(piece.i, data.i),
                  showSelected=minimization),
              data=data.table(cost.lines, seg.i="pruning"))+
    geom_point(aes(min.cost.mean, min.cost, color=data.i.fac,
                   showSelected=minimization),
               size=5,
               data=data.table(minima, seg.i="pruning"))+    
    facet_grid(. ~ seg.i, scales="free", labeller=function(var, val){
      paste(ifelse(val!="pruning", "segment", ""), val)
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
    guides(color="none"),
  data=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=800, height=300)+
    facet_grid(y ~ ., scales="free")+
    geom_point(aes(position, count),
               data=addY(data.dt, "count"))+
    geom_segment(aes(segment.start-0.45, min.cost.mean,
                     showSelected=minimization,
                     xend=segment.end+0.45, yend=min.cost.mean),
                 data=addY(data.minima, "count"),
                 color="green")+
    guides(color="none")+
    geom_point(aes(timestep, intervals),
               data=addY(data.intervals[total.segments==2,],"intervals"))+
    geom_tile(aes(timestep, total.segments, fill=optimal.cost,
                  clickSelects=minimization),
              data=addY(data.cost, "segments"))
)
minima.active <- data.minima[constraint=="active",]
if(nrow(minima.active)){
  viz$funModels <- viz$funModels+
    geom_point(aes(min.cost.mean, min.cost,
                   showSelected=minimization),
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
                   "at mean =",
                   round(min.cost.mean, 4),
                   "for",
                   seg.i,
                   "segment model up to data point",
                   timestep,
                   "change after",
                   data.i
                 ),
                 showSelected=minimization),
             size=5,
             data=data.minima)+
  geom_point(aes(mean, cost, 
                 showSelected=minimization),
             data=between.intervals)  
cost.active <- data.cost[constraint=="active",]
if(nrow(cost.active)){
  viz$data <- viz$data+
    geom_point(aes(timestep, total.segments, 
                   clickSelects=minimization),
               shape=21,
               color="black",
               fill="white",
               data=addY(cost.active, "segments"))
}
animint2dir(viz, "figure-constrained-PDPA-normal-real")

mlcost <- function(d.vec){
  m <- mean(d.vec)
  sum((d.vec-m)^2)
}
m34 <- mean(data.vec[3:4])
mlcost(data.vec[1:3])+(data.vec[4]-m34)
mlcost(data.vec[1:2])+mlcost(data.vec[3:4])
mlcost(data.vec[1:2])+mlcost(data.vec[3])
mlcost(data.vec[2:3])+mlcost(data.vec[1])

intervalsPlot <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(total.segments ~ .)+
  geom_point(aes(timestep, intervals),
             data=data.intervals)

pdf("figure-constrained-PDPA-normal-real.pdf")
print(intervalsPlot)
dev.off()

## FunctionalPruning <- list(
##   envelope=data.frame(envelope),
##   cost.lines=data.frame(cost.lines),
##   minima=data.frame(minima),
##   grid=data.frame(data.cost))
## save(FunctionalPruning, file="~/R/animint/data/FunctionalPruning.RData")
## prompt(FunctionalPruning, file="~/R/animint/man/FunctionalPruning.Rd")
