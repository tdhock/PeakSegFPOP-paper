source("packages.R")

data.vec <- c(1, 10, 14, 13)
mean.vec <- seq(min(data.vec), max(data.vec), by=0.5)
mean.grid <- data.table(expand.grid(
  seg1.mean=mean.vec,
  seg2.mean=mean.vec))
mean.grid[, cost.diff := ifelse(
  seg1.mean<seg2.mean,
  (seg2.mean-10)^2 - (seg1.mean-10)^2,
  Inf)]
minima <- data.table(
  seg1.mean=c(1, 5.5),
  seg2.mean=c(12, 14))
ggplot()+
  scale_fill_gradient2()+
  geom_tile(aes(seg1.mean, seg2.mean, fill=cost.diff),
            data=mean.grid)+
  geom_point(aes(seg1.mean, seg2.mean),
             shape=21,
             data=mean.grid[which.min(cost.diff),])+
  geom_point(aes(seg1.mean, seg2.mean),
             data=minima)+
  coord_equal()

fun.list <- list(
  end1=function(seg1.mean, seg2.mean){
    rbind(seg1.mean, seg2.mean, seg2.mean)
  },
  end2=function(seg1.mean, seg2.mean){
    rbind(seg1.mean, seg1.mean, seg2.mean)
  })
cost.rects.list <- list()
cost.minima.list <- list()
setkey(mean.grid, seg1.mean, seg2.mean)
setkey(minima, seg1.mean, seg2.mean)
for(fun.name in names(fun.list)){
  fun <- fun.list[[fun.name]]
  residual.mat <- mean.grid[, fun(seg1.mean, seg2.mean)-c(1, 10, 14)]
  cost.vec <- colSums(residual.mat * residual.mat)
  is.feasible <- mean.grid[, seg1.mean<seg2.mean]
  cost.vec[!is.feasible] <- Inf
  mean.grid$cost <- cost.vec
  cost.minima.list[[fun.name]] <- data.table(fun.name, mean.grid[minima])
  cost.rects.list[[fun.name]] <- data.table(fun.name, mean.grid)
}
cost.rects <- do.call(rbind, cost.rects.list)
cost.minima <- do.call(rbind, cost.minima.list)
dput(RColorBrewer::brewer.pal(Inf, "Blues"))
c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", 
  "#2171B5", "#08519C", "#08306B")
with.legend <- ggplot()+
  ggtitle(paste(
    "segmentation cost up to 3 data",
    "red = break after data point 1 less costly",
    "blue = break after data point 2 less costly",
    sep="\n"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ fun.name, labeller=function(var, val){
    sub("end", "contours: cost of break after ", val)
  })+
  coord_equal()+
  scale_x_continuous(breaks=c(1, 5.5, 14))+
  scale_y_continuous(breaks=c(1, 10, 12, 14))+
  scale_fill_gradient2()+
  geom_tile(aes(seg1.mean, seg2.mean, fill=cost.diff),
            data=mean.grid)+
  geom_contour(aes(seg1.mean, seg2.mean, z=cost, color=..level..),
               breaks=seq(20, 160, by=20),
               data=cost.rects)+
  geom_text(aes(seg1.mean, seg2.mean, color=cost, label=cost),
            data=cost.minima,
            vjust=-0.5,
            hjust=1)+
  geom_point(aes(seg1.mean, seg2.mean, color=cost),
             data=cost.minima)
(with.labels <- direct.label(with.legend, "top.pieces"))
cost.rects[seg1.mean==5.5 & seg2.mean==14,]
## The end1 model has a lower cost at the optimum for the end2 model.

pdf("figure-constrained-PDPA-normal-grid.pdf", w=13)
print(with.labels)
dev.off()

mean3.grid <- data.table(expand.grid(
  seg1.mean=mean.vec,
  seg2.mean=mean.vec,
  seg3.mean=mean.vec))[seg1.mean < seg2.mean & seg2.mean > seg3.mean,]
break.grid <- data.table(expand.grid(
  end1=1:3,
  end2=1:3))[end1<end2,]
all.cost.list <- list()
for(break.i in 1:nrow(break.grid)){
  break.row <- break.grid[break.i,]
  segs <- break.row[, data.table(
    start=c(1, end1+1, end2+1),
    end=c(end1, end2, 4))]
  segs[, size := end-start+1]
  mean.mat <- matrix(NA, 4, nrow(mean3.grid))
  for(seg.i in 1:nrow(segs)){
    seg <- segs[seg.i,]
    row.i.vec <- seg[, start:end]
    col.name <- paste0("seg", seg.i, ".mean")
    for(row.i in row.i.vec){
      mean.mat[row.i,] <- mean3.grid[[col.name]]
    }
  }
  residual.mat <- data.vec-mean.mat
  cost.vec <- colSums(residual.mat * residual.mat)
  all.cost.list[[break.i]] <- data.table(break.row, mean3.grid, cost=cost.vec)
}
all.cost <- do.call(rbind, all.cost.list)
all.cost[order(cost),]
## with the constraints seg1.mean < seg2.mean > seg3.mean the minimum
## error model is undefined
all.cost[seg3.mean==13,][order(cost),]
all.cost[seg1.mean==5.5 & seg2.mean==14 & seg3.mean==13,][order(cost),]

## Evaluate C13(1,12+x,12), show that the minimum is undefined.
curve(colSums((c(1,10,14,13)-rbind(1,12+x,12+x,12))^2))
