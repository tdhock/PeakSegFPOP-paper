source("packages.R")

load("whole.chrom.RData")
whole.chrom[, count := as.integer(coverage)]
whole.chrom[, zeros.after := c(chromStart[-1]-chromEnd[-.N], 0) ]

some <- whole.chrom
base.values <- some[, unique(sort(c(chromStart, chromEnd)))]
some.zeros <- data.table(
  chromStart=base.values[-length(base.values)],
  chromEnd=base.values[-1],
  count=0L)
setkey(some, chromStart)
setkey(some.zeros, chromStart)
some.zeros[some, count := as.integer(coverage)]

McGill0003_H3K4me3_chr1 <- data.frame(some.zeros)
##save(McGill0003_H3K4me3_chr1, file="~/R/cosegData/data/McGill0003_H3K4me3_chr1.RData", compress="xz")
prompt(McGill0003_H3K4me3_chr1, file="~/R/cosegData/man/McGill0003_H3K4me3_chr1.Rd")

N <- 1000000
one <- data.frame(some.zeros[1:N,])
seconds <- system.time({
  fpop <- PeakSegFPOPchrom(one, N)
})[["elapsed"]]
