source("packages.R")
data(intreg, package="animint")

##options(warn=2)

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
    new.dt.list <- list()
    prev.min.cost <- NULL
    row.i <- 1
    prev.min.mean <- dt$min.mean[1]
    while(row.i <= nrow(dt)){
      this.row <- dt[row.i,]
      gg <- ggplot()+
        geom_vline(xintercept=prev.min.mean)+
        geom_line(aes(mean, cost),
                  color="grey",
                  size=4,
                  data=getLines(dt))+
        geom_line(aes(mean, cost),
                  data=getLines(this.row))
      if(is.numeric(prev.min.cost)){
        gg <- gg+
          geom_hline(yintercept=prev.min.cost)
      }
      if(is.null(prev.min.cost)){
        ## Look for min achieved in this interval.
        mu <- getMinMean(this.row)
        if(mu <= this.row$min.mean){
          ## The minimum is achieved before this interval, so this
          ## function is always increasing in this interval. We don't
          ## need to store it.
          prev.min.cost <- quad(this.row, this.row$min.mean)
          prev.data.i <- this.row$data.i
        }else if(mu < this.row$max.mean){
          ## Minimum in this interval.
          new.row <- this.row
          new.row$min.mean <- prev.min.mean
          new.row$max.mean <- mu
          new.dt.list[[paste(row.i)]] <- new.row
          prev.min.mean <- mu
          prev.min.cost <- quad(this.row, mu)
          prev.data.i <- this.row$data.i
        }else{
          ## Minimum after this interval, so this function is
          ## decreasing on this entire interval, and so we can just
          ## store it as is.
          new.row <- this.row
          new.row$min.mean <- prev.min.mean
          new.dt.list[[paste(row.i)]] <- new.row
          prev.min.mean <- this.row$max.mean
        }
      }else{
        ## Look for a function with prev.min.cost.
        discriminant <- this.row[, linear^2-4*quadratic*(constant-prev.min.cost)]
        if(0 < discriminant){
          mu <- this.row[, (-linear-sqrt(discriminant))/(2*quadratic)]
          if(this.row[, min.mean < mu & mu < max.mean]){
            new.dt.list[[paste(row.i, "constant")]] <- data.table(
              quadratic=0,
              linear=0,
              constant=prev.min.cost,
              min.mean=prev.min.mean,
              max.mean=mu,
              data.i=prev.data.i)
            prev.min.cost <- NULL
            prev.min.mean <- mu
            row.i <- row.i-1
          }
        }
      }
      row.i <- row.i+1
    }
    new.dt.list[["last"]] <- data.table(
      quadratic=0,
      linear=0,
      constant=prev.min.cost,
      min.mean=prev.min.mean,
      max.mean=this.row$max.mean,
      data.i=prev.data.i)
    do.call(rbind, new.dt.list)
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
mirror <- function(dt){
  new.max <- -dt$min.mean
  new.min <- -dt$max.mean
  dt$min.mean <- new.min
  dt$max.mean <- new.max
  dt[, linear := -linear]
  dt[.N:1,]
}
less.more.min.list$more <- function(dt){
  mirror.input <- mirror(dt)
  mirror.output <- less.more.min.list$less(mirror.input)
  mirror(mirror.output)
}
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
### Return just the minimum interval on the intersection of row1 and
### row2.
CompareRows <- function(dt1, dt2, i1, i2){
  row1 <- dt1[i1,]
  row2 <- dt2[i2,]
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
              data=getLines(row2))
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
    ## The functions are exactly equal over the entire interval so we
    ## can return either one of them.
    row1$min.mean <- last.min.mean
    row1$max.mean <- first.max.mean
    return(row1)
  }
  ## They are not equal over the entire interval, but they may be
  ## equal on the left or right
  insignificant.list <- list(
    ##machine=.Machine$double.eps
  )
  cost1.right <- quad(row1, first.max.mean)
  if(row1$max.mean < row2$max.mean){
    cost1.next <- quad(dt1[i1+1,], first.max.mean)
    insignificant.list[["row1 next"]] <- abs(cost1.right-cost1.next)
  }
  cost2.right <- quad(row2, first.max.mean)
  if(row2$max.mean < row1$max.mean){
    cost2.next <- quad(dt2[i2+1,], first.max.mean)
    insignificant.list[["row2 next"]] <- abs(cost2.right-cost2.next)
  }
  cost1.left <- quad(row1, last.min.mean)
  if(row2$min.mean < row1$min.mean){
    cost1.prev <- quad(dt1[i1-1,], last.min.mean)
    insignificant.list[["row1 prev"]] <- abs(cost1.left-cost1.prev)
  }
  cost2.left <- quad(row2, last.min.mean)
  if(row1$min.mean < row2$min.mean){
    cost2.prev <- quad(dt2[i2-1,], last.min.mean)
    insignificant.list[["row2 prev"]] <- abs(cost2.prev-cost2.left)
  }
  insignificant.vec <- unlist(insignificant.list)
  thresh <- max(insignificant.vec)
  row.diff <- row1-row2
  row.diff$min.mean <- last.min.mean
  row.diff$max.mean <- first.max.mean
  discriminant <- row.diff[, linear^2 - 4*quadratic*constant]
  right.diff <- cost1.right-cost2.right
  if(thresh < abs(right.diff)){
    row1.min.on.right <- right.diff < 0
  }else{
    ## They are equal on the right limit, so use the first and second
    ## derivatives to see which is minimal just before the right
    ## limit. Do we need to check if they intersect before the right
    ## limit? Only if one is going up and the other is going down. And
    ## in that case we just need to check the -sqrt of the
    ## difference (since we know the +sqrt is on the right limit).
    deriv1.right <- row1[, 2*quadratic*first.max.mean + linear]
    deriv2.right <- row2[, 2*quadratic*first.max.mean + linear]
    sign1 <- sign(deriv1.right)
    sign2 <- sign(deriv2.right)
    if(sign1 != 0 && sign2 != 0 && sign1 != sign2){
      ## There could be a crossing point to the left.
      numerator <- - row.diff$linear - sqrt(discriminant)
      denominator <- 2*row.diff$quadratic
      mean.at.equal.cost <- numerator/denominator
      in.interval <-
        last.min.mean < mean.at.equal.cost &
        mean.at.equal.cost < first.max.mean
      if(in.interval){
        stop("intersection to the left of the equal right limit")
      }
    }
    row1.min.before.right <- if(sign1==sign2){
      row2$quadratic < row1$quadratic
    }else{
      deriv2.right < deriv1.right 
    }
    this.row <- if(row1.min.before.right)row1 else row2
    this.row$min.mean <- last.min.mean
    this.row$max.mean <- first.max.mean
    return(this.row)
  }
  left.diff <- cost1.left-cost2.left
  if(thresh < abs(left.diff)){
    row1.min.on.left <- left.diff < 0
  }else{
    ## Equal on the left.
    deriv1.left <- row1[, 2*quadratic*last.min.mean + linear]
    deriv2.left <- row2[, 2*quadratic*last.min.mean + linear]
    sign1 <- sign(deriv1.left)
    sign2 <- sign(deriv2.left)
    if(sign1 != 0 && sign2 != 0 && sign1 != sign2){
      ## There could be a crossing point to the right.
      numerator <- -row.diff$linear + sqrt(discriminant)
      denominator <- 2*row.diff$quadratic
      mean.at.equal.cost <- numerator/denominator
      in.interval <-
        last.min.mean < mean.at.equal.cost &
        mean.at.equal.cost < first.max.mean
      if(in.interval){
        stop("intersection to the right of the equal left limit")
      }
    }
    row1.min.after.left <- if(sign1==sign2){
      row1$quadratic < row2$quadratic
    }else{
      deriv1.left < deriv2.left
    }
    this.row <- if(row1.min.after.left)row1 else row2
    this.row$min.mean <- last.min.mean
    this.row$max.mean <- first.max.mean
    return(this.row)
  }
  ## The only remaining case is that the curves are equal neither on
  ## the left nor on the right of the interval. However they may be
  ## equal inside the interval, so let's check for that.
  mean.in.interval <- if(0 < discriminant){
    numerator <- -row.diff$linear + c(-1,1)*sqrt(discriminant)
    denominator <- 2*row.diff$quadratic
    mean.at.equal.cost <- numerator/denominator
    in.interval <-
      last.min.mean < mean.at.equal.cost &
      mean.at.equal.cost < first.max.mean
    mean.at.equal.cost[in.interval]
  }
  new.intervals <- if(length(mean.in.interval)==1){
    stopifnot(row1.min.on.left != row2.min.on.right)
    if(row1.min.on.left)rbind(row1,row2) else rbind(row2, row1)
  }else if(length(mean.in.interval)==2){
    stopifnot(row1.min.on.left == row2.min.on.right)
    if(row1.min.on.left){
      rbind(row1, row2, row1)
    }else{
      rbind(row2, row1, row2)
    }
  }else{
    ## functions do not cross in this interval.
    stopifnot(row1.min.on.right==row1.min.on.left)
    if(row1.min.on.right)row1 else row2
  }
  new.intervals$min.mean <- c(last.min.mean, mean.in.interval)
  new.intervals$max.mean <- c(mean.in.interval, first.max.mean)
  new.intervals
}
## Another implementation that checks equality of quadratic
## coefficients rather than cost function values.
sameFuns <- function(row1, row2){
  row1$quadratic==row2$quadratic &&
    row1$constant==row2$constant
}
CompareRows <- function(dt1, dt2, i1, i2){
  row1 <- dt1[i1,]
  row2 <- dt2[i2,]
  ggplot()+
    ##coord_cartesian(xlim=c(-0.1, 0), ylim=c(0,0.2))+
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
              data=getLines(row2))
  if(row1$min.mean < row2$min.mean){
    prev2 <- dt2[i2-1, ]
    same.at.left <- sameFuns(prev2, row1)
    last.min.mean <- row2$min.mean
  }else{
    last.min.mean <- row1$min.mean
    prev1 <- dt1[i1-1, ]
    same.at.left <- if(row2$min.mean < row1$min.mean){
      sameFuns(prev1, row2)
    }else{
      if(i1==1 && i2==1){
        FALSE
      }else{
        prev2 <- dt2[i2-1, ]
        sameFuns(prev1, prev2)
      }
    }
  }
  if(row1$max.mean < row2$max.mean){
    next1 <- dt1[i1+1,]
    same.at.right <- sameFuns(next1, row2)
    first.max.mean <- row1$max.mean
  }else{
    first.max.mean <- row2$max.mean
    same.at.right <- if(row2$max.mean < row1$max.mean){
      next2 <- dt2[i2+1,]
      sameFuns(row1, next2)
    }else{
      if(i1==nrow(dt1) && i2==nrow(dt2)){
        FALSE
      }else{
        next1 <- dt1[i1+1,]
        next2 <- dt2[i2+1,]
        sameFuns(next1, next2)
      }
    }
  }
  stopifnot(last.min.mean < first.max.mean)
  if(sameFuns(row1, row2)){
    ## The functions are exactly equal over the entire interval so we
    ## can return either one of them.
    row1$min.mean <- last.min.mean
    row1$max.mean <- first.max.mean
    return(row1)
  }
  row.diff <- row1-row2
  row.diff$min.mean <- last.min.mean
  row.diff$max.mean <- first.max.mean
  discriminant <- row.diff[, linear^2 - 4*quadratic*constant]
  cost2.left <- quad(row2, last.min.mean)
  cost1.left <- quad(row1, last.min.mean)
  cost1.right <- quad(row1, first.max.mean)
  cost2.right <- quad(row2, first.max.mean)
  if(!same.at.right){
    row1.min.on.right <- cost1.right < cost2.right
  }else{
    ## They are equal on the right limit, so use the first and second
    ## derivatives to see which is minimal just before the right
    ## limit. Do we need to check if they intersect before the right
    ## limit? Only if the sign of the slope of the more curved
    ## function is positive. And in that case we just need to check
    ## the -sqrt of the difference (since we know the +sqrt is on the
    ## right limit).
    deriv1.right <- row1[, 2*quadratic*first.max.mean + linear]
    deriv2.right <- row2[, 2*quadratic*first.max.mean + linear]
    sign1 <- sign(deriv1.right)
    sign2 <- sign(deriv2.right)
    maybe.cross <-
      (row1$quadratic < row2$quadratic && 0 < sign2) ||
      (row2$quadratic < row1$quadratic && 0 < sign1)
    if(1e-10 < discriminant && maybe.cross){
      ## There could be a crossing point to the left.
      numerator <- row.diff[, -linear - sign(quadratic)*sqrt(discriminant)]
      denominator <- 2*row.diff$quadratic
      mean.at.equal.cost <- numerator/denominator
      in.interval <-
        last.min.mean < mean.at.equal.cost &
        mean.at.equal.cost < first.max.mean
      if(in.interval){
        stop("intersection to the left of the equal right limit")
      }
    }
    row1.min.before.right <- if(sign1==sign2){
      cost1.left < cost2.left
    }else{
      deriv2.right < deriv1.right 
    }
    this.row <- if(row1.min.before.right)row1 else row2
    this.row$min.mean <- last.min.mean
    this.row$max.mean <- first.max.mean
    return(this.row)
  }
  if(!same.at.left){
    row1.min.on.left <- cost1.left < cost2.left
  }else{
    ## Equal on the left.
    deriv1.left <- row1[, 2*quadratic*last.min.mean + linear]
    deriv2.left <- row2[, 2*quadratic*last.min.mean + linear]
    sign1 <- sign(deriv1.left)
    sign2 <- sign(deriv2.left)
    maybe.cross <-
      (row1$quadratic < row2$quadratic && sign2 < 0) ||
      (row2$quadratic < row1$quadratic && sign1 < 0)
    if(1e-10 < discriminant && maybe.cross){
      ## There could be a crossing point to the right.
      numerator <- row.diff[, -linear + sign(quadratic)*sqrt(discriminant)]
      denominator <- 2*row.diff$quadratic
      mean.at.equal.cost <- numerator/denominator
      in.interval <-
        last.min.mean < mean.at.equal.cost &
        mean.at.equal.cost < first.max.mean
      if(in.interval){
        stop("intersection to the right of the equal left limit")
      }
    }
    row1.min.after.left <- if(sign1==sign2){
      cost1.right < cost2.right
    }else{
      deriv1.left < deriv2.left
    }
    this.row <- if(row1.min.after.left)row1 else row2
    this.row$min.mean <- last.min.mean
    this.row$max.mean <- first.max.mean
    return(this.row)
  }
  ## The only remaining case is that the curves are equal neither on
  ## the left nor on the right of the interval. However they may be
  ## equal inside the interval, so let's check for that.
  mean.in.interval <- if(0 < discriminant){
    numerator <- row.diff[, -linear + sign(quadratic)*c(-1,1)*sqrt(discriminant)]
    denominator <- 2*row.diff$quadratic
    mean.at.equal.cost <- numerator/denominator
    in.interval <-
      last.min.mean < mean.at.equal.cost &
      mean.at.equal.cost < first.max.mean
    mean.at.equal.cost[in.interval]
  }
  new.intervals <- if(length(mean.in.interval)==1){
    if(row1.min.on.left)rbind(row1,row2) else rbind(row2, row1)
  }else if(length(mean.in.interval)==2){
    if(row1.min.on.left){
      rbind(row1, row2, row1)
    }else{
      rbind(row2, row1, row2)
    }
  }else{
    ## functions do not cross in this interval.
    if(row1.min.on.right)row1 else row2
  }
  new.intervals$min.mean <- c(last.min.mean, mean.in.interval)
  new.intervals$max.mean <- c(mean.in.interval, first.max.mean)
  new.intervals
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
    cross.dt <- if(0 < discriminant){
      numerator <- -row.diff$linear + c(-1,1)*sqrt(discriminant)
      denominator <- 2*row.diff$quadratic
      mean.at.equal.cost <- numerator/denominator
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
MinEnvelope <- function(dt1, dt2){
  stopifnot(dt1[, min.mean < max.mean])
  stopifnot(dt2[, min.mean < max.mean])
  stopifnot(dt1[1, min.mean]==dt2[1, min.mean])
  stopifnot(dt1[.N, max.mean]==dt2[.N, max.mean])
  i1 <- 1
  i2 <- 1
  new.dt.list <- list()
  row.to.add <- NULL
  while(i1 <= nrow(dt1) && i2 <= nrow(dt2)){
    row1 <- dt1[i1,]
    row2 <- dt2[i2,]
    new.rows <- CompareRows(dt1, dt2, i1, i2)
    for(row.i in 1:nrow(new.rows)){
      new.row <- new.rows[row.i,]
      if(is.null(row.to.add)){
        row.to.add <- new.row
      }else{
        if(new.row$quadratic==row.to.add$quadratic){
          ## this is the same function as the previous one, so just make
          ## it extend further to the right.
          row.to.add$max.mean <- new.row$max.mean
        }else{
          ## new.row is a different function than the previous
          ## row.to.add, so now is the time to store it and move on.
          new.dt.list[[paste(i1, i2, row.i)]] <- row.to.add
          row.to.add <- new.row
        }
      }
    }
    if(row1$max.mean == new.row$max.mean){
      i1 <- i1+1
    }
    if(row2$max.mean == new.row$max.mean){
      i2 <- i2+1
    }
  }
  new.dt.list[["last"]] <- row.to.add
  do.call(rbind, new.dt.list)
}
## Test cases for MinEnvelope: list format is: input1, input2,
## [output]. If output is omitted then it is the same as input2.
min.env.test.list <- list(not.quasiconvex=list(
  data.table(
    quadratic = c(4, 0),
    linear = c(0.260150264548243, 0),
    constant = c(0.0199023137057983, 0.0156724286967657),
    min.mean = c(-0.0922074380970889, -0.0325187830685304),
    max.mean = c(-0.0325187830685304, 0.0954195650786824),
    data.i = c(4, 4)),
  data.table(
    quadratic = c(3, 2, 1),
    linear = c(0.0757353883540654, -0.0842153650871264, -0.132854723477952),
    constant = c(0.0114001020653698, 0.00507885284652773, 0.00703123698800913),
    min.mean = c(-0.0922074380970889, -0.0713260031733877, 0.026116725688849),
    max.mean = c(-0.0713260031733877, 0.026116725688849, 0.0954195650786824),
    data.i = c(1, 2, 3))),
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
                 10, 6, 8))),
  same.on.right=list(
    data.table(
      quadratic = c(0, 6),
      linear = c(0, -0.606547830137147),
      constant = c(0.0393353945617033, 0.05466457248854),
      min.mean = c(-0.0922074380970889, 0.0505456525114289),
      max.mean = c(0.0505456525114289, 0.15521264992094),
      data.i = c(9, 9)),
    data.table(
      quadratic = c(1, 3, 6),
      linear = c(0.129765703169709, 0.0275713233400502, -0.606547830137147),
      constant = c(0.0275566915219987, 0.0211556474213202, 0.05466457248854),
      min.mean = c(-0.0922074380970889, -0.0365259791823832, 0.105686525579533),
      max.mean = c(-0.0365259791823832, 0.105686525579533, 0.15521264992094),
      data.i = c(8, 6, 9))),
  same.on.left=list(
    data.table(
      quadratic = c(6, 0),
      linear = c(0.606547830137147, 0),
      constant = c(0.05466457248854, 0.0393353945617033),
      min.mean = c(-0.15521264992094, -0.0505456525114289),
      max.mean = c(-0.0505456525114289, 0.0922074380970889),
      data.i = c(9, 9)),
    data.table(
      quadratic = c(6, 3, 1),
      linear = c(0.606547830137147, -0.0275713233400502, -0.129765703169709),
      constant = c(0.05466457248854, 0.0211556474213202, 0.0275566915219987),
      min.mean = c(-0.15521264992094, -0.105686525579533, 0.0365259791823832),
      max.mean = c(-0.105686525579533, 0.0365259791823832, 0.0922074380970889),
      data.i = c(9, 6, 8))))
for(min.env.test.name in names(min.env.test.list)){
  test <- min.env.test.list[[min.env.test.name]]
  dt1 <- test[[1]]
  dt2 <- test[[2]]
  computed <- MinEnvelope(dt1, dt2)
  expected <- if(length(test)==3)test[[3]] else test[[2]]
  gg.test <- ggplot()+
    ggtitle(min.env.test.name)+
    scale_color_manual(values=c(
      input1="#4DAF4A",
      input2="#FF7F00", #orange
      computed="black",
      expected="grey"))+
    geom_line(aes(mean, cost, color=fun.i),
              size=4,
              data=data.table(getLines(expected), fun.i="expected"))+
    geom_line(aes(mean, cost, color=fun.i),
              size=2.5,
              data=data.table(getLines(dt1), fun.i="input1"))+
    geom_line(aes(mean, cost, color=fun.i),
              size=1.5,
              data=data.table(getLines(dt2), fun.i="input2"))+
    geom_line(aes(mean, cost, color=fun.i),
              size=0.5,
              data=data.table(getLines(computed), fun.i="computed"))
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
##data.vec <- data.vec[1:60]
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
max.segments <- 5
for(total.segments in 2:max.segments){
  prev.cost.model <- cost.models.list[[paste(total.segments-1, total.segments-1)]]
  if(total.segments %% 2){
    min.fun.name <- "more"
    env.transform <- mirror
  }else{
    min.fun.name <- "less"
    env.transform <- identity
  }
  min.fun <- less.more.min.list[[min.fun.name]]
  first.min <- min.fun(prev.cost.model)
  first.data <- gamma.dt[total.segments,]
  first.data$data.i <- total.segments-1
  cost.model <- AddFuns(first.data, first.min)
  cost.models.list[[paste(total.segments, total.segments)]] <- cost.model
  for(timestep in (total.segments+1):length(data.vec)){
    cat(sprintf("%4d / %4d segments %4d / %4d data points %d intervals\n",
                total.segments, max.segments, timestep, length(data.vec),
                nrow(cost.model)))
    prev.cost.model <- cost.models.list[[paste(total.segments-1, timestep-1)]]
    compare.cost <- min.fun(prev.cost.model)
    compare.cost$data.i <- timestep-1
    cost.model <- cost.models.list[[paste(total.segments, timestep-1)]]
    cost.minima <- Minimize(cost.model)
    compare.minima <- Minimize(compare.cost)
    gg <- ggplot()+
      ggtitle(paste(total.segments, "segments,", timestep, "data points"))
    if(nrow(compare.cost)){
      cost.lines.list[[
        paste(total.segments, timestep, "compare")]] <-
        data.table(total.segments, timestep,
                   cost.type="compare",
                   compare.cost.lines <- getLines(compare.cost))
      gg <- gg+
        geom_line(aes(mean, cost,
                      color=factor(data.i),
                      group=piece.i),
                  compare.cost.lines)
    }
    if(nrow(cost.model)){ # may be Inf over entire interval.
      cost.lines.list[[paste(total.segments, timestep)]] <-
        data.table(total.segments, timestep,
                   cost.type="model",
                   cost.model.lines <- getLines(cost.model))
      gg <- gg+
        geom_line(aes(mean, cost,
                      group=piece.i,
                      color=factor(data.i)),
                  cost.model.lines)
    }
    ## transformed.compare <- env.transform(compare.cost)
    ## transformed.model <- env.transform(cost.model)
    ## transformed.env <-
    ##   MinEnvelope(transformed.compare, transformed.model)
    ## one.env <- env.transform(transformed.env)
    one.env <- MinEnvelope(compare.cost, cost.model)
    stopifnot(one.env[, min.mean < max.mean])
    if(nrow(one.env)){
      envelope.list[[paste(total.segments, timestep)]] <-
        data.table(total.segments, timestep,
                   env.lines <- getLines(one.env))
      gg <- gg+
        geom_line(aes(mean, cost),
                  env.lines,
                  color="grey",
                  size=2,
                  alpha=0.5)
    }
    if(nrow(cost.minima)){
      minima.list[[paste(total.segments, timestep)]] <-
        data.table(total.segments, timestep, rbind(
          cost.minima, compare.minima))
    }
    ## Now that we are done with this step, we can perform the
    ## recursion by setting the new model of the cost to the min
    ## envelope, plus a new data point.
    cost.model <- AddFuns(one.env, gamma.dt[timestep,])
    cost.models.list[[paste(total.segments, timestep)]] <- cost.model
  }#for(timestep
}#for(total.segments
cost.lines <- do.call(rbind, cost.lines.list)
cost.lines[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
cost.lines[, data.i.fac := factor(data.i)]
minima <- do.call(rbind, minima.list)
minima[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
minima[, data.i.fac := factor(data.i)]
envelope <- do.call(rbind, envelope.list)
envelope[, data.i.fac := factor(data.i)]
envelope[, minimization := paste(
  total.segments, "segments up to data point", timestep)]

gg.pruning <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(timestep ~ total.segments, scales="free",
             labeller=function(var, val){
               if(var %in% c("total.segments", "timestep")){
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

save.image("figure-constrained-PDPA-normal-real.RData")

ti <- 35
gg.pruning <- ggplot()+
  ##coord_cartesian(xlim=c(-0.2, 0), ylim=c(0, 0.3))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(timestep ~ total.segments, scales="free",
             labeller=function(var, val){
               if(var %in% c("total.segments", "timestep")){
                 paste(var, "=", val)
               }else{
                 paste(val)
               }
             })+
  geom_line(aes(mean, cost),
            color="grey",
            size=2,
            data=envelope[timestep==ti,])+
  geom_line(aes(mean, cost, group=paste(data.i.fac, piece.i),
                color=data.i.fac),
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
      min.dt$constraint <- if(min.dt[, quad.min.mean == min.cost.mean]){
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
  data.table(dt, y=factor(y, c("data value", "intervals", "segments")))
}
largest.constant <- envelope[quadratic==0, max(constant)]
viz <- list(
  title=paste(
    "Constrained Pruned Dynamic Programming Algorithm"),
  funModels=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=800, height=300)+
    coord_cartesian(ylim=c(0, max(between.intervals$cost)))+
    geom_line(aes(mean, cost,
                  key=1,
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
                  ##key=quadratic-timestep,
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
                  data=data.table(timestep=seq_along(data.vec)),
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
      breaks=unique(c(seq(1, length(data.vec), by=10), length(data.vec)))),
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
