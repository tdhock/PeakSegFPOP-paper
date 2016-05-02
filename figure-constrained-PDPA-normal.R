source("packages.R")

data.vec <- c(1, 10, 14, 13)
gamma.dt <- data.table(
  quadratic=1,
  linear=-2*data.vec,
  constant=data.vec * data.vec)
C12 <- gamma.dt[1,]+gamma.dt[2,]
C13 <- C12+gamma.dt[3,]
C14 <- C13+gamma.dt[4,]
quad <- function(dt, x){
  dt[, quadratic*x*x + linear*x + constant]
}
## C22 = gamma_1^< + gamma_2 = gamma_2
C22 <- data.table(seg1.end=1, min.mean=1, max.mean=14, gamma.dt[2,])
## The first candidate to compare is (gamma_1+gamma_2)^<
C23.candidate <- rbind(
  data.table(seg1.end=2, min.mean=1, max.mean=5.5, C12),
  data.table(seg1.end=2,
    min.mean=5.5, max.mean=14, quadratic=0, linear=0, constant=quad(C12, 5.5)))
getLines <- function(dt.name){
  dt <- get(dt.name)
  line.list <- list()
  for(piece.i in 1:nrow(dt)){
    piece <- dt[piece.i,]
    mean.vec <- piece[, seq(min.mean, max.mean, l=50)]
    line.list[[piece.i]] <- data.table(
      piece.i,
      piece,
      mean=mean.vec,
      cost=quad(piece, mean.vec))
  }
  data.table(dt.name, do.call(rbind, line.list))
}

pruneNotMin <- function(dt1, dt2){
  row1 <- 1
  row2 <- 1
}

C23.lines <- rbind(
  getLines("C22"),
  getLines("C23.candidate"))
with.legend <- ggplot()+
  geom_line(aes(mean, cost, color=dt.name),
            data=C23.lines)
(with.labels <- direct.label(with.legend)+xlim(1, 20))
pdf("figure-constrained-PDPA-normal.pdf")
print(with.labels)
dev.off()

gamma3 <- data.table(0,0,0,gamma.dt[3,])
C23 <- C22+gamma3
C23.can <- C23.candidate+rbind(gamma3,gamma3)
C23.lines <- rbind(
  getLines("C23"),
  getLines("C23.can"))
with.legend <- ggplot()+
  geom_line(aes(mean, cost, color=dt.name),
            data=C23.lines)
(with.labels <- direct.label(with.legend)+xlim(1, 20))

C22.upper <- C22
C22.upper[, min.mean := 10]
C22more <- rbind(
  data.table(
    seg1.end=1,
    min.mean=1,
    max.mean=10,
    quadratic=0,
    linear=0,
    constant=0),
  C22.upper)
with.legend <- ggplot()+
  geom_line(aes(mean, cost, color=dt.name),
            data=getLines("C22more"))
direct.label(with.legend)

C33 <- C22more+data.table(0,0,0,rbind(gamma.dt[3,],gamma.dt[3,]))
ggplot()+
  geom_line(aes(mean, cost, color=piece.i),
            data=getLines("C33"))

C23.min.mean <- C23[,-linear/(2*quadratic)]
C23.min.cost <- quad(C23, C23.min.mean)
C23more <- data.table(
  seg1.end=NA,
  min.mean=c(1, C23.min.mean),
  max.mean=c(C23.min.mean, 14),
  quadratic=c(0, C23$quadratic),
  linear=c(0, C23$linear),
  constant=c(C23.min.cost, C23$constant))

C33.lines <- rbind(getLines("C33"),getLines("C23more"))
with.legend <- ggplot()+
  geom_line(aes(mean, cost, color=dt.name),
            data=C33.lines)
direct.label(with.legend)

## from the other direction, 13 14 10 1

## The first comparison is C12< versus C11<+gamma_2

