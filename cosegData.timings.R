data("McGill0003_H3K4me3_chr1", package="cosegData")
data("hg19.gap", package="cosegData")
library(coseg)
library(memtime)
library(data.table)

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

N <- 1e5
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
