N.peaks <- 2
Kmax <- N.peaks*2+1
N <- 50
peak.N <- 10
bkg.mean <- 100
library(data.table)
true.end.vec <- cumsum(c(N, rep(c(peak.N, N), times=N.peaks)))
true.vec <- true.end.vec[-length(true.end.vec)]+c(1,0)
error.dt.list <- list()
for(jump in seq(0, 50, l=10)){
  for(seed in 1:10){
    set.seed(seed)
    peak.mean <- bkg.mean+jump
    data.vec <- unlist(c(
      rpois(N, bkg.mean),
      lapply(1:N.peaks, function(i)c(rpois(peak.N, peak.mean), rpois(N, bkg.mean)))))
    chromEnd <- 1:length(data.vec)
    data.dt <- data.table(
      chromStart=chromEnd-1L,
      chromEnd,
      count=data.vec)
    fit <- PeakSegOptimal::PeakSegPDPAchrom(data.dt, as.integer(N.peaks))
    peaks <- subset(fit$segments, status=="peak" & peaks==N.peaks)
    unc <- Segmentor3IsBack::Segmentor(data.vec, Kmax=Kmax)
    end.vec <- unc@breaks[Kmax,]
    pred.mat <- cbind(
      constrained=as.integer(with(peaks, rbind(first, last))),
      unconstrained=end.vec[-length(end.vec)] + c(1,0))
    error.dt.list[[paste(seed, jump)]] <- data.table(
      seed, jump,
      error.bases=rowSums(abs(pred.mat-true.vec)),
      algo=colnames(pred.mat))
  }
}
error.dt <- do.call(rbind, error.dt.list)
stats.dt <- error.dt[, list(
  median=median(error.bases),
  mean=mean(error.bases),
  sd=sd(error.bases),
  qlo=quantile(error.bases, 0.25),
  qhi=quantile(error.bases, 0.75),
  N=.N
  ), by=list(jump, algo)]
library(ggplot2)
ggplot()+
  geom_ribbon(aes(
    jump, ymin=qlo, ymax=qhi, fill=algo), data=stats.dt, alpha=0.5)+
  geom_line(aes(
    jump, median, color=algo),
            data=stats.dt)

ggplot()+
  geom_ribbon(aes(
    jump, ymin=mean-sd/sqrt(N), ymax=mean+sd/sqrt(N), fill=algo), data=stats.dt, alpha=0.5)+
  geom_line(aes(
    jump, mean, color=algo),
            data=stats.dt)
