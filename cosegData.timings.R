source("packages.R")

data("McGill0003_H3K4me3_chr1", package="cosegData")
data("hg19.gap", package="cosegData")

last <- 1e5
w <- 1e4
l <- McGill0003_H3K4me3_chr1[last,]
some <- McGill0003_H3K4me3_chr1[last+(-w):w,]
ggplot()+
  geom_vline(xintercept=l$chromStart/1e3)+
  geom_step(aes(chromStart/1e3, count),
          color="grey50",
            data=some)

## field	example	SQL type 	info 	description
## bin 	585	smallint(6) 	range 	Indexing field to speed chromosome range queries.
## chrom 	chr1	varchar(255) 	values 	Reference sequence chromosome or scaffold
## chromStart 	0	int(10) unsigned 	range 	start position in chromosome
## chromEnd 	10000	int(10) unsigned 	range 	end position in chromosome
## ix 	1	int(11) 	range 	index count of this fragment (obsolete/useless)
## n 	N	char(1) 	values 	'N' for gaps of known size, 'U' for gaps of unknown size
## size 	10000	int(10) unsigned 	range 	size of gap
## type 	telomere	varchar(255) 	values 	scaffold, contig, clone, fragment, etc.
## bridge 	no	varchar(255) 	values 	yes, no, mrna, bacEndPair, etc.
## colClasses <- c(
##   bin="NULL",
##   chrom="factor",
##   chromStart="integer",
##   chromEnd="integer",
##   ix="NULL",
##   n="NULL",
##   size="integer",
##   type="factor",
##   bridge="factor")
## hg19.gap <- read.table("~/Downloads/gap.txt.gz", colClasses=colClasses, col.names=names(colClasses))
## save(hg19.gap, file="~/R/cosegData/data/hg19.gap.RData")
## prompt(hg19.gap, file="~/R/cosegData/man/hg19.gap.Rd")

chr1.gap <- data.table(hg19.gap)[chrom=="chr1",]
chr1.cov <- data.table(McGill0003_H3K4me3_chr1)
setkey(chr1.cov, chromStart, chromEnd)
setkey(chr1.gap, chromStart, chromEnd)
over.dt <- foverlaps(chr1.cov, chr1.gap, nomatch=0L)

ggplot()+
  geom_segment(aes(chromStart/1e3, 0, xend=chromEnd/1e3, yend=0),
               data=chr1.gap)+
  geom_point(aes(i.chromStart/1e3, count),
             color="red",
             data=over.dt)+
  geom_point(aes(chromStart/1e3, 0),
             shape=1,
             data=chr1.gap)+
  geom_rect(aes(xmin=i.chromStart/1e3, xmax=i.chromEnd/1e3,
                ymin=0, ymax=count),
            data=over.dt)

non.zero <- over.dt[0 < count,]
ggplot()+
  geom_segment(aes(chromStart/1e3, 0, xend=chromEnd/1e3, yend=0),
               data=chr1.gap)+
  geom_point(aes(i.chromStart/1e3, count),
             color="red",
             data=non.zero)+
  geom_point(aes(chromStart/1e3, 0),
             shape=1,
             data=chr1.gap)

some.gaps <- chr1.gap[chromStart %in% non.zero$chromStart,]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("chromStart",scales="free")+
  geom_segment(aes(chromStart/1e3, 0, xend=chromEnd/1e3, yend=0),
               data=some.gaps)+
  geom_point(aes(i.chromStart/1e3, count),
             color="red",
             data=non.zero)+
  geom_point(aes(chromStart/1e3, 0),
             shape=1,
             data=some.gaps)+
  geom_rect(aes(xmin=i.chromStart/1e3, xmax=i.chromEnd/1e3,
                ymin=0, ymax=count),
            data=non.zero)

mid <- 121485434
w <- 1e6
some <- chr1.cov[!(chromStart < mid-w | mid+w < chromEnd),]
ggplot()+
  coord_cartesian(ylim=c(0, 100))+
  geom_step(aes(chromStart/1e3, count),
            color="grey50",
            data=some)+
  geom_segment(aes(chromStart/1e3, 0, xend=chromEnd/1e3, yend=0),
               size=2,
               data=chr1.gap[13:15,])

between.gaps <- chr1.gap[, data.table(problemStart=chromEnd[-.N], problemEnd=chromStart[-1])]
between.gaps[, bases := problemEnd-problemStart]
chr1.problems <- between.gaps[0 < bases,]
setkey(chr1.problems, problemStart, problemEnd)
chr1.problems.cov <- foverlaps(chr1.cov, chr1.problems, nomatch=0L)
chr1.problems.cov[chromStart < problemStart,]
chr1.problems.cov[problemEnd < chromEnd,]
problem.stats <- chr1.problems.cov[, list(
  sum.bases=sum(chromEnd-chromStart),
  data=.N
  ), by=.(problemStart, problemEnd)]
problem.stats[, bases := problemEnd-problemStart]
max(problem.stats$bases)
max(problem.stats$data)

all.between <-
  data.table(hg19.gap)[, data.table(
    problemStart=chromEnd[-.N],
    problemEnd=chromStart[-1]), by=chrom]
all.between[, max(problemEnd-problemStart)]

biggest <- chr1.problems.cov[(problemEnd-problemStart)==max(problemEnd-problemStart),]
biggest[chromStart < problemStart, chromStart := problemStart]
biggest[problemEnd < chromEnd, chromEnd := problemEnd]
biggest[chromEnd <= chromStart,]
biggest[, stopifnot(chromStart < chromEnd)]

## infer a gap track: all lines with zero counts which are longer than
## the biggest line with a positive count?
biggest.nonzero <- chr1.cov[0 < count, max(chromEnd-chromStart)]
maybe.gaps <- chr1.cov[count==0 & biggest.nonzero < chromEnd-chromStart, ]
maybe.gaps[, bases := chromEnd-chromStart]
other.problems <- maybe.gaps[, data.table(problemStart=chromEnd[-.N], problemEnd=chromStart[-1])]
other.problems[, stopifnot(problemStart < problemEnd)]
other.problems[, bases := problemEnd-problemStart]
## The inferred gap track heuristic looks like it is making too man
## segmentation problems, at least for these data.
ggplot()+
  geom_segment(aes(problemStart/1e3, 0,
                   xend=problemEnd/1e3, yend=0),
               size=1,
               data=chr1.problems)+
  geom_point(aes(problemStart/1e3, 0,
                   xend=problemEnd/1e3, yend=0),
               data=chr1.problems)+
  geom_segment(aes(problemStart/1e3, 1,
                   xend=problemEnd/1e3, yend=1),
               color="red",
               data=other.problems[38311 <= bases,])+
  geom_point(aes(problemStart/1e3, 1),
               color="red",
             data=other.problems[38311 <= bases,])

mid <- 121485434-5e7
w <- 1e5
some <- chr1.cov[!(chromStart < mid-w | mid+w < chromEnd),]
other.some <- other.problems[!(problemStart < mid-w | mid+w < problemEnd),]
some.gaps <- chr1.gap[!(chromStart < mid-w | mid+w < chromEnd),]
gg <- ggplot()+
  geom_step(aes(chromStart/1e3, count),
            color="grey50",
            data=some)+
  geom_segment(aes(problemStart/1e3, 0, xend=problemEnd/1e3, yend=0),
               color="red",
               data=other.some)+
  geom_point(aes(problemStart/1e3, 0, xend=problemEnd/1e3, yend=0),
               color="red",
             data=other.some)
if(nrow(some.gaps)){
  gg <- gg+  geom_segment(aes(chromStart/1e3, 0, xend=chromEnd/1e3, yend=0),
               size=2,
                          data=some.gaps)
}
print(gg)

N <- 10^6
divisor <- 10^3.5

pen.info.list <- list()
N.vec <- as.integer(10^seq(4, 6, by=0.5))
div.vec <- as.integer(10^seq(0, 3, by=0.5))
for(N in N.vec){
  some <- biggest[1:N,]
  for(divisor in div.vec){
    lambda <- N/divisor
    cat(sprintf("N=%d divisor=%f\n", N, divisor))
    for(rep.i in 1:2){
      info <- memtime({
        fpop <- PeakSegFPOPchrom(some, lambda)
      })
      pen.info.list[[paste(N, divisor, rep.i)]] <- data.table(
        algorithm="in.memory",
        fpop$loss, lambda, N, rep.i, seconds=info$time[["elapsed"]],
        memory.kilobytes=info$memory["max.increase", "kilobytes"],
        disk.kilobytes=0)
    }
  }
}

N <- nrow(biggest)
biggest$chrom <- "chr1"
for(N in c(N.vec, nrow(biggest))){
  some <- biggest[1:N,]
  some[, chromStart1 := chromStart + 1L]
  setkey(some, chrom, chromStart1, chromEnd)
  fwrite(some[, .(chrom, chromStart, chromEnd, count)], "memtest.bedGraph", col.names=FALSE, sep="\t")
  for(divisor in div.vec){
    lambda <- N/divisor
    cat(sprintf("N=%d divisor=%f\n", N, divisor))
    cmd <- paste("PeakSegFPOP memtest.bedGraph", lambda)
    for(rep.i in 1:2){
      unlink("tmp.db")
      info <- memtime({
        system(cmd)
      })
      bed.file <- paste0("memtest.bedGraph_penalty=", lambda, "_segments.bed")
      segs <- fread(bed.file)
      setnames(segs, c("chrom", "segStart", "segEnd", "status", "mean"))
      segs[, segStart1 := segStart+1L]
      setkey(segs, chrom, segStart1, segEnd)
      over.dt <- foverlaps(some, segs, nomatch=0L)
      stopifnot(nrow(over.dt)==nrow(some))
      ploss <- over.dt[, PoissonLoss(count, mean, chromEnd-chromStart)]
      peaks <- sum(segs$status=="peak")
      feasible <- all(diff(segs$mean) != 0)
      loss.dt <- data.table(
        segments=nrow(segs),
        peaks,
        penalized.loss=ploss+peaks*lambda,
        feasible)
      pen.info.list[[paste("on.disk", N, divisor, rep.i)]] <- data.table(
        algorithm="on.disk",
        loss.dt, lambda, N, rep.i, seconds=info$time[["elapsed"]],
        memory.kilobytes=info$memory["max.increase", "kilobytes"],
        disk.kilobytes=file.size("tmp.db")/1024)
    }
  }
}

cosegData.timings <- do.call(rbind, pen.info.list)
save(cosegData.timings, file="cosegData.timings.RData")
