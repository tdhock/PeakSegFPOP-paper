source("packages.R")

big.bedGraph <- "~/projects/chip-seq/data/McGill0003.1.brainTumor.H3K4me1.coverage.bedGraph"
big <- fread(big.bedGraph)
setnames(big, c("chrom", "chromStart", "chromEnd", "coverage"))
min.cov <- min(big$coverage)
stopifnot(0 < min.cov)
big[, coverage := coverage/min.cov]
whole.chrom <- big[chrom=="chr1",]

save(whole.chrom, file="whole.chrom.RData")
