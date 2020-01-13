source("packages.R")

ploss <- function(dt, x){
  ## need to make a new data table, otherwise ifelse may only get one
  ## element, and return only one element.
  new.dt <- data.table(dt, x)
  new.dt[, ifelse(Log==0, 0, Log*log(x)) + Linear*x + Constant]
}
pderiv <- function(dt, x){
  dt[, Linear+Log/x]
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
      cost=ploss(piece, mean.vec))
  }
  do.call(rbind, line.list)
}
getMinMean <- function(dt){
  dt[, -Log/Linear]
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
      Linear=row1$Linear+row2$Linear,
      Log=row1$Log+row2$Log,
      Constant=row1$Constant+row2$Constant,
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
    prev.data.i <- NULL
    row.i <- 1
    prev.min.mean <- dt$min.mean[1]
    while(row.i <= nrow(dt)){
      this.row <- dt[row.i,]
      if(is.null(prev.min.cost)){
        ## Look for min achieved in this interval.
        mu <- getMinMean(this.row)
        if(mu <= this.row$min.mean){
          ## The minimum is achieved before this interval, so this
          ## function is always increasing in this interval. We don't
          ## need to store it.
          prev.min.cost <- ploss(this.row, this.row$min.mean)
          prev.data.i <- this.row$data.i
        }else if(mu < this.row$max.mean){
          ## Minimum in this interval.
          new.row <- this.row
          new.row$min.mean <- prev.min.mean
          new.row$max.mean <- mu
          new.dt.list[[paste(row.i)]] <- new.row
          prev.min.mean <- mu
          prev.min.cost <- ploss(this.row, mu)
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
        if(this.row$Log==0){
          ## degenerate linear case.
          if(this.row$Linear < 0){
            ## decreasing linear function
            stop("this should never happen")
          }else{
            ##increasing linear function, so will not intersect the
            ##constant below.
          }
        }else{
          discriminant <-
            this.row[, Linear/Log*exp((prev.min.cost-Constant)/Log)]
          if(-1/exp(1) < discriminant){
            ## Since the Log constant is negative, the principal branch
            ## W(,0) results in the smaller of the two mean values.
            mu <- this.row[, Log/Linear*LambertW::W(discriminant, 0)]
            if(this.row[, min.mean < mu & mu < max.mean]){
              new.dt.list[[paste(row.i, "constant")]] <- data.table(
                Linear=0,
                Log=0,
                Constant=prev.min.cost,
                min.mean=prev.min.mean,
                max.mean=mu,
                data.i=prev.data.i)
              prev.min.cost <- NULL
              prev.min.mean <- mu
              row.i <- row.i-1
            }
          }#if(there are two roots
        }#if(degenerate linear) else
      }#if(looking for a min)else
      row.i <- row.i+1
    }
    if(!is.null(prev.data.i)){
      new.dt.list[["last"]] <- data.table(
        Linear=0,
        Log=0,
        Constant=prev.min.cost,
        min.mean=prev.min.mean,
        max.mean=this.row$max.mean,
        data.i=prev.data.i)
    }
    do.call(rbind, new.dt.list)
  }, more=function(dt){
    new.dt.list <- list()
    prev.min.cost <- NULL
    prev.data.i <- NULL
    row.i <- nrow(dt)
    prev.max.mean <- dt$max.mean[row.i]
    while(1 <= row.i){
      this.row <- dt[row.i,]
      if(is.null(prev.min.cost)){
        ## Look for min achieved in this interval.
        mu <- getMinMean(this.row)
        if(this.row$max.mean <= mu){
          ## The minimum is achieved after this interval, so this
          ## function is always decreasing in this interval. We don't
          ## need to store it.
          prev.min.cost <- ploss(this.row, this.row$max.mean)
          prev.data.i <- this.row$data.i
        }else if(this.row$min.mean < mu){
          ## Minimum in this interval.
          new.row <- this.row
          new.row$max.mean <- prev.max.mean
          new.row$min.mean <- mu
          new.dt.list[[paste(row.i)]] <- new.row
          prev.max.mean <- mu
          prev.min.cost <- ploss(this.row, mu)
          prev.data.i <- this.row$data.i
        }else{
          ## Minimum before this interval, so this function is
          ## increasing on this entire interval, and so we can just
          ## store it as is.
          new.row <- this.row
          new.row$max.mean <- prev.max.mean
          new.dt.list[[paste(row.i)]] <- new.row
          prev.max.mean <- this.row$min.mean
        }
      }else{
        ## Look for a function with prev.min.cost.
        mu <- if(this.row$Log==0){
          ## degenerate case where the function is linear.
          this.row[, (prev.min.cost - Constant)/Linear]
        }else{
          discriminant <-
            this.row[, Linear/Log*exp((prev.min.cost-Constant)/Log)]
          if(-1/exp(1) < discriminant){
            ## Since the Log constant is negative, the non-principal
            ## branch W(,-1) results in the larger of the two mean
            ## values.
            browser(expr=!is.finite(discriminant))
            this.row[, Log/Linear*LambertW::W(discriminant, -1)]
          }#if(-1/e < discriminant
        }
        if(is.numeric(mu) && this.row[, min.mean < mu & mu < max.mean]){
          new.dt.list[[paste(row.i, "constant")]] <- data.table(
            Linear=0,
            Log=0,
            Constant=prev.min.cost,
            min.mean=mu,
            max.mean=prev.max.mean,
            data.i=prev.data.i)
          prev.min.cost <- NULL
          prev.max.mean <- mu
          row.i <- row.i+1
        }#if(mu in interval
      }#if(is.null(prev.min.cost)else
      row.i <- row.i-1
    }
    if(!is.null(prev.data.i)){
      new.dt.list[["last"]] <- data.table(
        Linear=0,
        Log=0,
        Constant=prev.min.cost,
        min.mean=this.row$min.mean,
        max.mean=prev.max.mean,
        data.i=prev.data.i)
    }
    do.call(rbind, rev(new.dt.list))
  })

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
sameFuns <- function(row1, row2){
  row1$Linear==row2$Linear &&
    row1$Log==row2$Log &&
    row1$Constant==row2$Constant
}
CompareRows <- function(dt1, dt2, i1, i2){
  row1 <- dt1[i1,]
  row2 <- dt2[i2,]
  ## ggplot()+
  ##   ##coord_cartesian(xlim=c(-0.1, 0), ylim=c(0,0.2))+
  ##   geom_line(aes(mean, cost, color=fun.i),
  ##             size=2,
  ##             data=data.table(getLines(dt1), fun.i=factor(1)))+
  ##   geom_line(aes(mean, cost, color=fun.i),
  ##             size=1,
  ##             data=data.table(getLines(dt2), fun.i=factor(2)))+
  ##   geom_line(aes(mean, cost),
  ##             linetype="dashed",
  ##             data=getLines(row1))+
  ##   geom_line(aes(mean, cost),
  ##             linetype="dotted",
  ##             size=1,
  ##             data=getLines(row2))
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
  is.same <- sameFuns(row1, row2)
  if(is.same){
    ## The functions are exactly equal over the entire interval so we
    ## can return either one of them.
    row1$min.mean <- last.min.mean
    row1$max.mean <- first.max.mean
    return(row1)
  }
  row.diff <- row1-row2
  row.diff$min.mean <- last.min.mean
  row.diff$max.mean <- first.max.mean
  if(row.diff$Log==0){
    if(row.diff$Linear==0){
      ## They are offset by a constant.
      new.row <- if(row.diff$Constant < 0)row1 else row2
      new.row$min.mean <- last.min.mean
      new.row$max.mean <- first.max.mean
      return(new.row)
    }
    if(row.diff$Constant==0){
      ## The only difference is the Linear term.
      new.row <- if(row.diff$Linear < 0)row1 else row2
      new.row$min.mean <- last.min.mean
      new.row$max.mean <- first.max.mean
      return(new.row)
    }
    mean.at.equal.cost <- row.diff[, -Constant/Linear]
    root.in.interval <-
      last.min.mean < mean.at.equal.cost &&
      mean.at.equal.cost < first.max.mean
    if(root.in.interval){
      new.rows <- if(0 < row.diff$Linear){
        rbind(row1, row2)
      }else{
        rbind(row2, row1)
      }
      new.rows$min.mean <- c(last.min.mean, mean.at.equal.cost)
      new.rows$max.mean <- c(mean.at.equal.cost, first.max.mean)
      return(new.rows)
    }else{
      new.row <- if(mean.at.equal.cost < last.min.mean)row1 else row2
      new.row$min.mean <- last.min.mean
      new.row$max.mean <- first.max.mean
      return(new.row)
    }
  }
  ## cost2.left <- ploss(row2, last.min.mean)
  ## cost1.left <- ploss(row1, last.min.mean)
  ## cost1.right <- ploss(row1, first.max.mean)
  ## cost2.right <- ploss(row2, first.max.mean)
  cost.diff.right <- ploss(row.diff, first.max.mean)
  cost.diff.left <- ploss(row.diff, last.min.mean)
  discriminant <- row.diff[, Linear/Log*exp(-Constant/Log)]
  ## The discriminant could be -Inf, if the exp argument is larger
  ## than about 700.
  two.roots <- -1/exp(1) + 1e-8 < discriminant
  if(!same.at.right){
    row1.min.on.right <- cost.diff.right < 0
  }else{
    ## They are equal on the right limit, so use the first and second
    ## derivatives to see which is minimal just before the right
    ## limit. Do we need to check if they intersect before the right
    ## limit? Only if the sign of the slope of the more curved
    ## function is positive. And in that case we just need to check
    ## the smaller root.
    deriv1.right <- pderiv(row1, first.max.mean)
    deriv2.right <- pderiv(row2, first.max.mean)
    sign1 <- sign(deriv1.right)
    sign2 <- sign(deriv2.right)
    maybe.cross <-
      (row2$Log < row1$Log && 0 < sign2) ||
      (row1$Log < row2$Log && 0 < sign1)
    if(two.roots && maybe.cross){
      ## There could be a crossing point to the left.
      root.left <- row.diff[, Log/Linear*LambertW::W(discriminant, 0)]
      mean.at.equal.cost <- root.left
      in.interval <-
        last.min.mean < mean.at.equal.cost &
        mean.at.equal.cost < first.max.mean
      if(in.interval){
        new.rows <- if(cost.diff.left < 0){
          rbind(row1, row2)
        }else{
          rbind(row2, row1)
        }
        new.rows$min.mean <- c(last.min.mean, mean.at.equal.cost)
        new.rows$max.mean <- c(mean.at.equal.cost, first.max.mean)
        return(new.rows)
      }
    }
    row1.min.before.right <- if(sign1==sign2){
      cost.diff.left < 0
    }else{
      deriv2.right < deriv1.right 
    }
    this.row <- if(row1.min.before.right)row1 else row2
    this.row$min.mean <- last.min.mean
    this.row$max.mean <- first.max.mean
    return(this.row)
  }
  if(!same.at.left){
    row1.min.on.left <- cost.diff.left < 0
  }else{
    ## Equal on the left.
    deriv1.left <- pderiv(row1, last.min.mean)
    deriv2.left <- pderiv(row2, last.min.mean)
    sign1 <- sign(deriv1.left)
    sign2 <- sign(deriv2.left)
    maybe.cross <-
      (row2$Log < row1$Log && sign2 < 0) ||
      (row1$Log < row2$Log && sign1 < 0)
    if(two.roots && maybe.cross){
      ## There could be a crossing point to the right.
      root.right <- row.diff[, Log/Linear*LambertW::W(discriminant, -1)]
      mean.at.equal.cost <- root.right
      in.interval <-
        last.min.mean < mean.at.equal.cost &
        mean.at.equal.cost < first.max.mean
      if(in.interval){
        new.rows <- if(cost.diff.right < 0){
          rbind(row2, row1)
        }else{
          rbind(row1, row2)
        }
        new.rows$min.mean <- c(last.min.mean, mean.at.equal.cost)
        new.rows$max.mean <- c(mean.at.equal.cost, first.max.mean)
        return(new.rows)
      }
    }
    row1.min.after.left <- if(sign1==sign2){
      cost.diff.right < 0
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
  mean.in.interval <- if(two.roots){
    root.right <- row.diff[, Log/Linear*LambertW::W(discriminant, -1)]
    root.left <- row.diff[, Log/Linear*LambertW::W(discriminant, 0)]
    mean.at.equal.cost <- sort(c(root.left, root.right))
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
        if(sameFuns(new.row, row.to.add)){
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


input.dt <- data.table(count=c(3, 9, 18, 15, 20, 2), weight=1)
library(animint)

cost.lines.list <- list()
minima.list <- list()
envelope.list <- list()
cost.models.list <- list()

min.mean <- min(input.dt$count)
max.mean <- max(input.dt$count)
gamma.dt <- input.dt[, data.table(
  Linear=weight,
  Log=-count*weight,
  ##Constant=ifelse(count==0, 0, weight*count*(log(count)-1))
  Constant=0
  )]
cum.weight <- input.dt[, cumsum(weight)]
C1.dt <- cumsum(gamma.dt)
gamma.dt$min.mean <- C1.dt$min.mean <- min.mean
gamma.dt$max.mean <- C1.dt$max.mean <- max.mean
gamma.dt$data.i <- C1.dt$data.i <- 0

data.colors <- c(
  "1"="#E41A1C",
  "2"="#377EB8",
  "3"="#4DAF4A",
  "4"="#984EA3",
  "5"="#FF7F00",
  "6"="#FFFF33",
  "0"="black",
  "#A65628", "#F781BF", "#999999")

for(data.i in 1:nrow(C1.dt)){
  cost.models.list[[paste(1, data.i)]] <- C1.dt[data.i,]
}
max.segments <- 5
pdftikz <- function(pre, g, w=3){
  pdf(paste0(pre, ".pdf"), 5, 3)
  print(g)
  dev.off()
  tikz(paste0(pre, ".tex"), width=w, height=1.5)
  print(g)
  dev.off()
}
compare.cost.list <- list()
for(total.segments in 2:max.segments){
  prev.cost.model <-
    cost.models.list[[paste(total.segments-1, total.segments-1)]]
  if(total.segments %% 2){
    min.fun.name <- "more"
  }else{
    min.fun.name <- "less"
  }
  min.fun <- less.more.min.list[[min.fun.name]]
  first.min.total <- min.fun(prev.cost.model)
  gg.prev <- ggplot()+
    ggtitle(paste(total.segments, "segments,", total.segments, "data points"))+
    scale_x_continuous()+
    scale_color_manual("prev seg end", values=data.colors)+
    geom_line(aes(mean, cost),
              getLines(first.min.total),
              color="grey",
              size=3)+
    geom_line(aes(mean, cost,
                  color=factor(data.i),
                  group=piece.i),
              data=getLines(prev.cost.model))
  i <- unique(prev.cost.model$data.i)
  if(length(i)==1 && i==0){
    gg.prev <- gg.prev+guides(color="none")
  }
  ##pdftikz(sprintf("figure-PeakSegPDPA-demo-minlessmore-%dsegments-%ddata", total.segments, total.segments), gg.prev)
  gg.prev <- ggplot()+
    ggtitle(paste(total.segments, "segments,", total.segments, "data points"))+
    scale_x_continuous()+
    scale_color_manual(values=c(
      "prev cost"="black",
      "constrained min"="grey",
      "unconstrained min"="red"))+
    geom_line(aes(mean, cost, color=fun.type),
              data=data.table(
                fun.type="constrained min",
                getLines(first.min.total)),
              size=3)+
    geom_line(aes(mean, cost,
                  color=fun.type,
                  group=piece.i),
              data=data.table(
                fun.type="prev cost",
                getLines(prev.cost.model)))
  first.min <- Minimize(prev.cost.model)
  if(nrow(first.min)){
    gg.prev <- gg.prev+    geom_hline(aes(yintercept=min.cost, color=fun.type),
               data=data.table(
                 fun.type="unconstrained min",
                 first.min))
  }
  ##pdftikz(sprintf("figure-PeakSegPDPA-demo-mincompare-%dsegments-%ddata", total.segments, total.segments), gg.prev)
  first.data <- gamma.dt[total.segments,]
  first.data$data.i <- total.segments-1
  cost.model.total <- AddFuns(first.data, first.min.total)
  cost.models.list[[paste(total.segments, total.segments)]] <- cost.model.total
  changepoint <- paste0(
    "change\npoint $t_", total.segments-1, "$")
  for(timestep in (total.segments+1):length(input.dt$count)){
    prev.cost.model <- cost.models.list[[paste(total.segments-1, timestep-1)]]
    gg.prev <- ggplot()+
      ggtitle(paste(total.segments, "segments,", timestep, "data points"))+
      scale_color_manual(changepoint, values=data.colors)+
      geom_line(aes(mean, cost,
                    color=factor(data.i),
                    group=piece.i),
                data=getLines(prev.cost.model))
    ## pdf(sprintf("figure-PeakSegPDPA-demo-cost-%dsegments-%ddata.pdf", total.segments, timestep), 5, 3)
    ## print(gg.prev)
    ## dev.off()
    compare.cost <- min.fun(prev.cost.model)
    compare.cost$data.i <- timestep-1
    cost.model <- cost.models.list[[paste(total.segments, timestep-1)]]
    cat(sprintf("%4d / %4d segments %4d / %4d data points %d intervals\n",
                total.segments, max.segments, timestep, length(input.dt$count),
                nrow(cost.model)))
    cost.minima <- Minimize(cost.model)
    compare.minima <- Minimize(compare.cost)
    gg <- ggplot()+
      scale_x_continuous()+
      scale_color_manual(changepoint, values=data.colors)+
      ggtitle(paste(total.segments, "segments,", timestep, "data points"))
    if(nrow(compare.cost)){
      compare.cost.list[[paste(total.segments, timestep)]] <- compare.cost
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
      gg.prev <- ggplot()+
        ggtitle(paste(total.segments, "segments,", timestep, "data points"))+
        scale_x_continuous()+
        scale_color_manual(changepoint, values=data.colors)+
        geom_line(aes(mean, cost),
                  compare.cost.lines,
                  color="grey",
                  size=3)+
        geom_line(aes(mean, cost,
                      color=factor(data.i),
                      group=piece.i),
                  data=getLines(prev.cost.model))
      i <- unique(prev.cost.model$data.i)
      if(length(i)==1 && i==0){
        gg.prev <- gg.prev+guides(color="none")
      }
      ##pdftikz(sprintf("figure-PeakSegPDPA-demo-minlessmore-%dsegments-%ddata", total.segments, timestep), gg.prev)
      gg.prev <- ggplot()+
        ggtitle(paste(total.segments, "segments,", timestep, "data points"))+
        scale_x_continuous()+
        scale_color_manual(values=c(
                             "prev cost"="black",
                             "constrained min"="grey",
                             "unconstrained min"="red"))+
        geom_line(aes(mean, cost, color=fun.type),
                  data=data.table(
                    fun.type="constrained min",
                    compare.cost.lines),
                  size=3)+
        geom_line(aes(mean, cost,
                      color=fun.type,
                      group=piece.i),
                  data=data.table(
                    fun.type="prev cost",
                    getLines(prev.cost.model)))
      if(nrow(compare.minima)){
        gg.prev <- gg.prev+          geom_hline(aes(yintercept=min.cost, color=fun.type),
                                                data=data.table(
                                                  fun.type="unconstrained min",
                                                  compare.minima))
      }
      ##pdftikz(sprintf("figure-PeakSegPDPA-demo-mincompare-%dsegments-%ddata", total.segments, timestep), gg.prev)
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
    one.env <- MinEnvelope(cost.model, compare.cost)
    stopifnot(one.env[, min.mean < max.mean])
    if(nrow(one.env)){
      envelope.list[[paste(total.segments, timestep)]] <-
        data.table(total.segments, timestep,
                   env.lines <- getLines(one.env))
      gg <- gg+
        geom_line(aes(mean, cost),
                  env.lines,
                  color="grey",
                  size=3)
    }
    gg <- ggplot()+
      geom_line(aes(mean, cost),
                env.lines,
                color="grey",
                size=3)+
      geom_line(aes(mean, cost,
                    color=factor(data.i),
                    group=piece.i),
                compare.cost.lines)+
      geom_line(aes(mean, cost,
                    group=piece.i,
                    color=factor(data.i)),
                cost.model.lines)+
      scale_x_continuous()+
      scale_color_manual(changepoint, values=data.colors)+
      ggtitle(paste(total.segments, "segments,", timestep, "data points"))
    ##pdftikz(sprintf("figure-PeakSegPDPA-demo-minenv-%dsegments-%ddata", total.segments, timestep), gg)
    if(nrow(cost.minima)){
      minima.list[[paste(total.segments, timestep)]] <-
        data.table(total.segments, timestep, rbind(
          cost.minima, compare.minima))
    }
    ## Now that we are done with this step, we can perform the
    ## recursion by setting the new model of the cost to the min
    ## envelope, plus a new data point.
    cost.model.total <- AddFuns(one.env, gamma.dt[timestep,])
    cost.lines.list[[paste(
      total.segments, timestep, "result")]] <- data.table(
        total.segments=total.segments+1, timestep=timestep+1,
        cost.type="result",
        getLines(cost.model.total))
    cost.models.list[[paste(total.segments, timestep)]] <- cost.model.total
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

PeakSegDP::cDPA(input.dt$count, maxSegments=5L)

tsegs <- 3
ti <- 4
print(compare.cost.list[[paste(tsegs, ti)]])#this is the grey C^\geq_{2,3} in the figure.
cost.lines[, total.cost := cost]
minima[, total.min.cost := min.cost]
gg.pruning <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_line(aes(mean, total.cost, group=paste(cost.type, piece.i), linetype=cost.type),
            size=3,
            color="grey80",
            data=cost.lines[timestep==ti & tsegs==total.segments & cost.type == "compare"])+
  geom_text(aes(
    mean, cost, label=label),
            data=data.table(
              mean=15, cost=-14.42*3),
              label=c(
                "best cost in\n2 segments\nand 3 data\n$C_{2,3}(u_2)$"))+
  geom_text(aes(
    mean, cost, label=label),
            color="red",
            data=data.table(
              mean=7, cost=-43.4),
              label=c(
                "CDPA discards possibility
of non-decreasing change after $t_1=2$"))+
  geom_segment(aes(x, y, xend=xend, yend=yend),
            color="red",
               linetype="dotted",
               data=data.table(
                 x=10, y=-14.48*3, xend=17, yend=-43.51),
               arrow=grid::arrow(length=grid::unit(0.1, "in"), type="closed"))+
  geom_text(aes(
    mean, cost, label=label),
            color="red",
            data=data.table(
              mean=7, cost=-14.5*3),
              label=c(
                "CDPA computes this scalar min cost,
only considering a possible
non-decreasing change after $t_1=1$"))+
  geom_segment(aes(x, y, xend=xend, yend=yend),
            color="red",
               linetype="dotted",
               data=data.table(
                 x=10, y=-14.5*3, xend=13, yend=-43.56),
               arrow=grid::arrow(length=grid::unit(0.1, "in"), type="closed"))+
  geom_text(aes(
    mean, cost, label=label),
            color="grey50",
            ##vjust=1,
            data=data.table(
              mean=11, cost=-43.6,
              label=c(
                "GPDPA considers both possible changes $t_1\\in\\{1,2\\}$ by computing $C^\\geq_{2,3}(u_3)$,
the functional min cost of a non-increasing change after $t_2=3$")))+
  geom_line(aes(mean, total.cost),
            size=1.25,
            data=cost.lines[timestep==ti & tsegs==total.segments & cost.type == "result"])+
  geom_line(aes(mean, total.cost, group=paste(cost.type, piece.i), linetype=cost.type,
                color=data.i.fac),
            size=1,
            data=cost.lines[timestep==ti & tsegs==total.segments & cost.type == "result"])+
  geom_point(aes(min.cost.mean, total.min.cost),
             size=3,
             shape=1,
             color="red",
             data=minima[timestep==ti & tsegs==total.segments,])+
  coord_cartesian(ylim=c(-43.625, -43.2))+
  scale_color_manual(values=c("1"="deepskyblue", "2"="violet"), guide="none")+
  scale_linetype_discrete(labels=c(result="$C_{2,3}$", compare="C^\\geq{2,3}"), guide="none")+
  geom_text(aes(
    mean, cost, color=data.i.fac, label=label),
            data=data.table(
              mean=c(9, 18), cost=c(-43.25,-43.35), data.i.fac=factor(c(1, 2)),
              label=c(
                "cost of non-decreasing\nchange after $t_1=1$\n$\\ell(y_3, u_2)+C_{2,2}(u_2)$",
                "cost of non-decreasing\nchange after $t_1=2$\n$\\ell(y_3, u_2)+C^\\leq_{1,2}(u_2)$")))+
  ylab("cost")+
  xlab("segment mean")
print(gg.pruning)
tikz("figure-CDPA-fails.tex", width=6, height=4)
print(gg.pruning)
dev.off()
##system("pdflatex aoas-supplementary")
