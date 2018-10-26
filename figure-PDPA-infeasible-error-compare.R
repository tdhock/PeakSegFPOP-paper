library(data.table)
library(ggplot2)

load("Segmentor.peaks.error.RData")
load("Segmentor.infeasible.error.RData")
load("PDPA.infeasible.error.RData")
load("PDPA.peaks.error.RData")
load("dp.peaks.error.RData")
load("dp.peaks.matrices.RData")

CDPA.error.list <- list()
for(chunk.name in names(dp.peaks.error)){
  chunk.error <- dp.peaks.error[[chunk.name]]
  for(cell.type in names(chunk.error)){
    type.dt <- data.table(
      chunk.name,
      chunk.error[[cell.type]])
    names(type.dt)[3] <- "peaks"
    CDPA.error.list[[paste(chunk.name, cell.type)]] <- type.dt
  }
}
CDPA.error <- do.call(rbind, CDPA.error.list)

baseline.vec <- c(
  macs.trained="MACS",
  hmcan.broad.trained="HMCanBroad")
baseline.error.list <- list()
for(set.name in names(dp.peaks.matrices)){
  set.chunks <- dp.peaks.matrices[[set.name]]
  for(chunk.name in names(set.chunks)){
    chunk.algos <- set.chunks[[chunk.name]]
    for(algo.name in names(baseline.vec)){
      algo <- baseline.vec[[algo.name]]
      algo.mat <- chunk.algos[[algo.name]]
      min.vec <- apply(algo.mat, 1, min)
      baseline.error.list[[paste(set.name, chunk.name, algo.name)]] <- data.table(
        algo,
        chunk.name,
        sample.id=names(min.vec),
        min.errors=min.vec)
    }
  }
}
baseline.error <- do.call(rbind, baseline.error.list)

CDPA.error[, segments := as.integer(paste(peaks))*2+1]
PDPA.peaks.error[, segments := as.integer(paste(peaks))*2+1]
Segmentor.peaks.error[, segments := as.integer(paste(peaks))*2+1]
col.name.vec <- names(PDPA.infeasible.error)
DT <- function(algo, dt){
  data.table(algo, dt[, ..col.name.vec])
}
all.error <- rbind(
  DT("CDPA", CDPA.error),
  DT("G.ig", PDPA.peaks.error),
  DT("G.jo", PDPA.infeasible.error),
  DT("S.ig", Segmentor.peaks.error),
  DT("S.jo", Segmentor.infeasible.error))

all.error[, table(chunk.name, algo)]

all.totals <- all.error[, list(
  total.fp=sum(fp),
  total.fn=sum(fn),
  total.errors=sum(fp+fn),
  possible.fp=sum(possible.fp),
  possible.fn=sum(possible.tp),
  labels=.N
  ), by=list(algo, chunk.name, sample.id, peaks, segments)]

all.totals[, table(chunk.name, algo)]

all.min <- rbind(baseline.error, all.totals[, list(
  min.errors=min(total.errors)
  ), by=list(algo, chunk.name, sample.id)])

all.min[, table(chunk.name, algo)]

min.wide <- dcast(all.min, chunk.name+sample.id~algo)

## as expected there are some (81) problems where the peak joining
## achieves fewer errors than ignoring the infeasible models with
## equality constraints.
min.wide[G.jo < G.ig]
min.wide[G.ig < G.jo]

## Somewhat strangely there are 35 problems when CDPA gets fewer label
## errors, and also 35 problems when GPDPA with join when infeasible
## gets fewer label errors.
min.wide[CDPA < G.jo]
min.wide[G.jo < CDPA]

## Even more strangely the distribution of differences is symmetrical!
min.wide[, table(CDPA-G.jo)]
min.wide[, set.name := sub("/.*", "", chunk.name)]
min.wide[, table(set.name, CDPA-G.jo)]
## > min.wide[, table(diff)]
## diff
##   -2   -1    0    1    2 
##    3   32 2682   32    3 
## >

min.wide[, list(
  count=.N),
  by=list(GPDPA.jo.ig=G.jo-G.ig)]

## ignore vs join
g <- function(x, Comparison){
  dt <- data.table(table(x), Comparison)
  names(dt)[1] <- "diff"
  dt
}
count.tall <- min.wide[, rbind(
  g(G.jo-G.ig, "Join-Ignore\nGPDPA"),
  g(S.jo-S.ig, "Join-Ignore\nPDPA"),
  g(G.jo-MACS, "GPDPAjoin-MACS\n"),
  g(G.jo-HMCanBroad, "GPDPAjoin-HMCanBroad\n"),
  g(G.jo-CDPA, "GPDPAjoin-CDPA\n"),
  g(G.jo-S.jo, "GPDPA-PDPA\n(Join post-processing)"),
  g(G.ig-S.ig, "GPDPA-PDPA\n(Ignore post-processing)"))]
count.tall[, diff.fac := factor(diff, (-10):2)]
(count.wide <- dcast(count.tall, Comparison ~diff.fac, value.var="N"))
## library(xtable)
## xt <- xtable(count.wide, align="lrrrrrrrrrr", digits=0)
## print(
##   xt,
##   file="PDPA-infeasible-error-compare.tex",
##   include.rownames=FALSE, floating=FALSE)
library(namedCapture)
pattern <- paste0(
  "(?<before>[^-]+)",
  "-",
  "(?<after>[a-zA-Z]+)")
levs <- rev(unique(count.tall$Comparison))
count.tall[, comp.fac := factor(Comparison, levs)]
count.tall[, comp.int := as.integer(comp.fac)]
count.tall[, diff.int := as.integer(diff)]
count.tall[, first := sub("\n.*", "", Comparison)]
count.tall[, last := sub(".*\n", "", Comparison)]
count.tall[, element := ifelse(
  first=="Join-Ignore",
  last,
  Comparison)]
count.tall[, panel := ifelse(
  first=="Join-Ignore",
  "Comparing
post-processing
methods", ifelse(
  first=="GPDPA-PDPA",
  "Comparing
constrained/
unconstrained", "Comparing
GPDPAjoin
with baselines"))]
better.dt <- count.tall[, {
  match.df <- str_match_named(paste(comp.fac), pattern)
  data.table(
    text=c(match.df[, "after"], match.df[, "before"]),
    total=c(sum(N[diff.int>0]), sum(N[diff.int<0])),
    diff.int=c(3.5, -8.5),
    hjust=c(0, 1))
}, by=list(comp.fac, panel, last)]
##((text=="Ignore") | (last=="(Ignore post-processing)"&text=="PDPA"))
better.dt[total==0, `:=`(
  diff.int=0.5,
  text="No difference")]
gg <- ggplot()+
  theme_bw()+
  theme(
    panel.margin=grid::unit(0, 'lines'),
    panel.grid.minor=element_blank())+
  facet_grid(panel ~ ., scales="free", space="free")+
  geom_vline(aes(
    xintercept=xint),
    data=data.table(xint=0),
    size=2)+
  geom_tile(aes(
    diff.int, comp.fac, fill=log10(N)),
    data=count.tall)+
  geom_text(aes(
    diff.int, comp.fac, label=N),
    data=count.tall)+
  scale_fill_gradient(
    "log10(
segmentation
problems)",
low="white", high="red")+
  scale_x_continuous(
    "Difference in minimum number of incorrect labels per segmentation problem",
    limits=c(-10, 5),
    breaks=unique(count.tall$diff.int))+
  scale_y_discrete("Comparison")+
  geom_text(aes(
    diff.int, comp.fac,
    label=ifelse(
      grepl("difference", text),
      ##text,
      "",
      paste0(text, " better
", total, " problems")),
    hjust=hjust),
    size=3,
    data=better.dt)
pdf("figure-PDPA-infeasible-error-compare.pdf", 12, 4)
print(gg)
dev.off()

