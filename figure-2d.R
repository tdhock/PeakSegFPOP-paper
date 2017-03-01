## two data sequences.
x1 <- c(1, 2, 3, 2)
x2 <- c(-2, -3, 3, 1)

## Cost of a mean level m1 for the first data sequence, and a mean
## level m2 for the second data sequence.
C1 <- function(m1,m2){
  (x1[1]-m1)^2+(x2[1]-m2)^2
}
C22 <- function(m1,m2){
  (x1[2]-m1)^2+(x2[2]-m2)^2
}
C12 <- function(m1,m2){
  sum(
    c(x1[1:2]-m1)^2,
    c(x2[1:2]-m2)^2)
}
C12(1.5, -2.5)#==1
M22 <- function(m1,m2){
  ##M22(u)=min{C22(u),min_x C12(x)}
  c22 <- C22(m1,m2)
  ifelse(c22 < 1, c22, 1)
}
M22_1 <- function(m1){
  c22 <- (x1[2]-m1)^2
  ifelse(c22<0.5, c22, 0.5)
}
M22_2 <- function(m2){
  c22 <- (x2[2]-m2)^2
  ifelse(c22<0.5, c22, 0.5)
}
library(data.table)
m.grid <- seq(-3,3,by=0.1)
grid.dt <- data.table(expand.grid(m1=m.grid, m2=m.grid))
grid.dt[, cost := M22(m1, m2)]
grid.dt[, cost1 := M22_1(m1)]
grid.dt[, cost2 := M22_2(m2)]

library(ggplot2)
## True min cost.
ggplot()+
  geom_tile(aes(m1, m2, fill=cost), data=grid.dt)

## Cost computed by adding min on each dimension.
ggplot()+
  geom_tile(aes(m1, m2, fill=cost1+cost2), data=grid.dt)

library(rgl)
