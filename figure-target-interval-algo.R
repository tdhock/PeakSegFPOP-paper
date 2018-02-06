library(penaltyLearning)
data(neuroblastoma, package="neuroblastoma", envir=environment())
one <- subset(neuroblastoma$profiles, profile.id==599 & chromosome=="14")
max.segments <- 1000
fit <- Segmentor3IsBack::Segmentor(one$logratio, model=2, Kmax=max.segments)
lik.df <- data.frame(lik=fit@likelihood, segments=1:max.segments)
times.list <- list()
for(n.segments in seq(10, max.segments, by=10)){
  some.lik <- lik.df[1:n.segments,]
  some.times <- microbenchmark::microbenchmark(
    R=pathR <- with(some.lik, modelSelectionR(lik, segments, segments)),
    C=pathC <- with(some.lik, modelSelectionC(lik, segments, segments)),
    times=5)
  times.list[[paste(n.segments)]] <- data.frame(n.segments, some.times)
}
times <- do.call(rbind, times.list)
## modelSelectionR and modelSelectionC should give identical results.
identical(pathR, pathC)
## However, modelSelectionC is much faster (linear time complexity)
## than modelSelectionR (quadratic time complexity).
library(ggplot2)
ggplot()+
  geom_point(aes(log10(n.segments), log10(time/1e9), color=expr), data=times)
ms <- modelSelection(some.lik, "lik", "segments")
with(ms, plot(segments, lik))
