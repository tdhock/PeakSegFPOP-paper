source("packages.R")
h <- 2 #inches
y.min <- 0
y.max <- 4
sloss <- function(dt, x){
  dt[, Quadratic*x*x + Linear*x + Constant]
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
      cost=sloss(piece, mean.vec))
  }
  do.call(rbind, line.list)
}
funPiece <- function(Quadratic, Linear, Constant, min.mean, max.mean){
  data.table(Quadratic, Linear, Constant, min.mean, max.mean)
}
C11 <- rbind(
  funPiece(1, -4, 4, 0, 4)
  )
C11min <- rbind(
  funPiece(0, 0, 0, 0, 4))
C11minless <- rbind(
  funPiece(1, -4, 4, 0, 2),
  funPiece(0, 0, 0, 2, 4))
C22 <- rbind(
  funPiece(2, -6, 5, 0, 2),
  funPiece(1, -2, 1, 2, 4))
C22un <- rbind(
  funPiece(1, -2, 1, 0, 4))
C12minless <- rbind(
  funPiece(2, -6, 5, 0, 1.5),
  funPiece(0, 0, 0.5, 1.5, 4))
small.title <-   theme(plot.title=element_text(size=10))
gg <- ggplot()+
  ggtitle("Cost of 2 segments up to data point 2")+
  theme_bw()+small.title+
  xlab("segment mean $\\mu$")+
  geom_line(aes(mean, cost), data=getLines(C22), color="red", size=3)+
  geom_line(aes(mean, cost), data=getLines(C22un), color="grey50", size=1.5)+
  geom_text(aes(x=1.5,y=8.5,label="constrained"),color="red")+
  geom_text(aes(x=1.5,y=7,label="$C_{2,2}(\\mu)=$"),color="red")+
  geom_text(aes(x=3,y=1.1,label="$(\\mu-1)^2$"),color="grey50", hjust=0.5, vjust=0)+
  geom_text(aes(x=3,y=0,label="unconstrained"),color="grey50", hjust=0.5, vjust=0)+
  geom_text(aes(x=1.5,y=5.5,label="$C^\\leq_{1,1}(\\mu)+(\\mu-1)^2$"),color="red")
tikz("figure-compare-cost.tex", 2.3, h)
print(gg)
dev.off()
gg <- ggplot()+
  ggtitle("Min-less computation for data point 1")+
  scale_x_continuous("segment mean $\\mu$", limits=c(-0.5, 4))+
  geom_line(aes(mean, cost), data=getLines(C11minless), color="red", size=3)+
  geom_line(aes(mean, cost), data=getLines(C11min), color="grey50", size=1.5)+
  geom_line(aes(mean, cost), data=getLines(C11), color="black", size=1)+
  geom_text(aes(x=3.5,y=3.5,label="$C_{1,1}(\\mu) = (\\mu-2)^2$"),color="black", hjust=1)+
  geom_text(aes(x=-Inf,y=0.3,label="$\\min_\\mu C_{1,1}(\\mu)$"),color="grey50", hjust=0, vjust=0)+
  geom_text(aes(x=2,y=2.6,label="$C^\\leq_{1,1}(\\mu)=$"),color="red")+
  geom_text(aes(x=2,y=2,label="$\\min_{x\\leq \\mu} C_{1,1}(x)$"),color="red")+
  theme_bw()+small.title
tikz("figure-compare-unconstrained.tex", 2.3, h)
print(gg)
dev.off()



