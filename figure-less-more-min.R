source("packages.R")

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
Inf.dt <- data.table(
  quadratic=numeric(),
  linear=numeric(),
  constant=numeric(),
  min.mean=numeric(),
  max.mean=numeric(),
  data.i=numeric())
less.more.min.list <- list(
  "cost $f$"=identity,
  ## less.min=function(dt){
  ##   if(1 < nrow(dt)){
  ##     stop("TODO implement more general less equal min computation")
  ##   }
  ##   mu <- getMinMean(dt)
  ##   cost <- quad(dt, mu)
  ##   data.table(
  ##     quadratic=c(dt$quadratic, 0),
  ##     linear=c(dt$linear, 0),
  ##     constant=c(dt$constant, cost),
  ##     min.mean=c(min.mean, mu),
  ##     max.mean=c(mu, max.mean),
  ##     data.i=dt$data.i)[min.mean!=max.mean,]
  ## },
  "$f^{\\geq}$"=function(dt){
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
          min.mean,
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
          min.mean=c(min.mean, mean.at.min),
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
min.mean <- 1
max.mean <- 20
cost.model <- data.table(
    quadratic=c(2,1),
    linear=c(-54, -28),
    constant=c(365, 196),
    min.mean=c(1,13),
    max.mean=c(13,20),
    data.i=1
)
min.lines.list <- list()
for(fun.type in names(less.more.min.list)){
  min.fun <- less.more.min.list[[fun.type]]
  fun.model <- min.fun(cost.model)
  min.lines.list[[fun.type]] <- data.table(
    fun.type, getLines(fun.model))
}
min.lines <- do.call(rbind, min.lines.list)

gg.cost <- ggplot()+
  xlab("segment mean $\\mu$")+
  ylab("function value $f(\\mu)$ or $f^{\\geq}(\\mu)$")+
  scale_size_manual(values=c("cost $f$"=2, "$f^{\\geq}$"=1))+
  geom_line(aes(mean, cost, color=fun.type, size=fun.type),
            data=min.lines)

tikz("figure-less-more-min.tex", width=4, height=3)
print(gg.cost)
dev.off()

