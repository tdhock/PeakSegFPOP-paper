source("packages.R")

cost.dt.list <- list(more=data.table(
  ## 3 segs, 4 intervals.
  Linear = c(5500, 5000, 4500, 3000, 1000, 500), 
  Log = c(-10192000, -10072000, -9641500, -7862000, -3606500, -2157500),
  Constant = c(68287393.5425206, 67749716.8717195, 65363783.303153,
               54686408.0843032, 26874159.7358435, 17010715.3104453),
  min.mean = c(197, 240.022373420469, 1556.83732507916, 
               1720.08889691817, 3388.58410993867, 4893.60479447855),
  max.mean = c(240.022373420469, 1556.83732507916, 1720.08889691817,
               3388.58410993867, 4893.60479447855, 21424),
  data.i = c(10, 1, 2, 5, 9, 10)),
  ## 4 segs, data point 30
  less=data.table(
    Linear = c(500, 2000, 8000, 8500, 1000),
    Log = c(-98500, -1560500, -8587000, -9588000, -665500),
    Constant = c(3544830.38422225, 11902005.4111316, 53071092.3660107,
                 59630489.0895849, 7377394.61123679),
    min.mean = c(197, 515.525665401257, 570.232907695563, 1438.59795952913,
                 2175.55624844772),
    max.mean = c(515.525665401257, 570.232907695563, 
                 1438.59795952913, 2175.55624844772, 21424),
    data.i = c(28, 25, 13, 12, 28)))

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

ploss <- function(dt, x){
  ## need to make a new data table, otherwise ifelse may only get one
  ## element, and return only one element.
  new.dt <- data.table(dt, x)
  new.dt[, ifelse(Log==0, 0, Log*log(x)) + Linear*x + Constant]
}
getMinMean <- function(dt){
  dt[, -Log/Linear]
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

cost.lines.list <- list()
for(min.fun.name in c("less", "more")){
  cost.dt <- cost.dt.list[[min.fun.name]]
  min.fun <- less.more.min.list[[min.fun.name]]
  min.dt <- min.fun(cost.dt)
  cost.lines.list[[paste("min", min.fun.name)]] <- data.table(
    min.fun.name, fun.type="min", getLines(min.dt))
  cost.lines.list[[paste("cost", min.fun.name)]] <- data.table(
    min.fun.name, fun.type="cost", getLines(cost.dt))
}
cost.lines <- do.call(rbind, cost.lines.list)
cost.lines[, fun.type.fac := factor(fun.type, c("min", "cost"))]


lab.vec <-
  c("min-less\n$C^{\\leq}(\\mu)$",
    "cost\n$C(\\mu)$")
gg.less <- ggplot()+
  ggtitle("min-less operator")+
  theme_bw()+
  xlab("segment mean $\\mu$")+
  ylab("cost")+
  scale_size_manual("function", values=c(cost=1, min=3), labels=lab.vec)+
  scale_color_manual("function", values=c(cost="black", min="grey"), labels=lab.vec)+
  guides(color=guide_legend(keyheight=3))+
  geom_line(aes(mean/100, cost/1e6, color=fun.type.fac, size=fun.type.fac),
            data=cost.lines[min.fun.name=="less",])+
  coord_cartesian(xlim=c(0, 25), ylim=c(1, 5))
tikz("figure-1-min-less-operator.tex", h=2, w=3)
print(gg.less)
dev.off()

lab.vec <-
  c("min-more\n$C^{\\geq}(\\mu)$",
    "cost\n$C(\\mu)$")
gg.more <- ggplot()+
  ggtitle("min-more operator")+
  theme_bw()+
  xlab("segment mean $\\mu$")+
  ylab("cost")+
  scale_size_manual("function", values=c(cost=1, min=3), labels=lab.vec)+
  scale_color_manual("function", values=c(cost="black", min="grey"), labels=lab.vec)+
  guides(color=guide_legend(keyheight=3))+
  geom_line(aes(mean/1000, cost/1e6, color=fun.type.fac, size=fun.type.fac),
            data=cost.lines[min.fun.name=="more",])+
  coord_cartesian(xlim=c(0, 15), ylim=c(0, 4))
tikz("figure-1-min-more-operator.tex", h=2, w=3)
print(gg.more)
dev.off()

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ min.fun.name)+
  scale_size_manual(values=c(cost=1, min=3))+
  scale_color_manual(values=c(cost="black", min="grey"))+
  geom_line(aes(mean, cost, color=fun.type.fac, size=fun.type.fac),
            data=cost.lines)

pdf("figure-1-min-operators.pdf")
print(gg)
dev.off()
