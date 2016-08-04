data("McGill0003_H3K4me3_chr1", package="cosegData")
library(coseg)
library(memtime)
library(data.table)
pen.info.list <- list()
for(N in 10^seq(4, 6, by=0.5)){
  some <- McGill0003_H3K4me3_chr1[1:N,]
  for(divisor in 10^seq(0, 3, by=0.5)){
    lambda <- N/divisor
    cat(sprintf("N=%f divisor=%f\n", N, divisor))
    for(rep.i in 1:2){
      info <- memtime({
        fpop <- PeakSegFPOPchrom(some, lambda)
      })
      pen.info.list[[paste(N, divisor, rep.i)]] <- data.table(
        fpop$loss, lambda, N, rep.i, seconds=info$time[["elapsed"]],
        kilobytes=info$memory["max.increase", "kilobytes"])
    }
  }
}
cosegData.timings <- do.call(rbind, pen.info.list)
save(cosegData.timings, file="cosegData.timings.RData")
