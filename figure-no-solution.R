source("packages.R")
candidate.vec <- c(5.5, 14, 14, 13)
data.vec <- c(1, 10, 14, 13)
mu <- mean(data.vec[-1])
solution.vec <- c(1, mu, mu, mu)
a.vec <- seq(0, 1, l=20)
better.mat.list <- list()
loss.vec.list <- list()
for(a in a.vec){
  better <- candidate.vec*(1-a) + solution.vec*a
  better.mat.list[[paste(a)]] <- better
  loss.vec.list[[paste(a)]] <- PoissonLoss(data.vec, better)
}
better.mat <- do.call(rbind, better.mat.list)
plot(do.call(rbind, loss.vec.list))

