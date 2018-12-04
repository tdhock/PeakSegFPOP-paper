source("packages.R")

(objs <- load("Segmentor.peaks.error.RData"))
rm(Segmentor.model.list)
gc()
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
        min.errors=min.vec,
        labels=chunk.algos$regions)
    }
  }
}
baseline.error <- do.call(rbind, baseline.error.list)

CDPA.error[, segments := as.integer(paste(peaks))*2+1]
PDPA.peaks.error[, segments := as.integer(paste(peaks))*2+1]
Segmentor.peaks.error[, segments := as.integer(paste(peaks))*2+1]
col.name.vec <- c(
  "chunk.name", "sample.id", "peaks", "segments", "chromStart", 
  "chromEnd", "annotation", "tp", "possible.tp", "fp", "possible.fp", 
  "fp.status", "fn", "fn.status", "status")
DT <- function(algo, dt){
  data.table(algo, dt[, ..col.name.vec])
}
all.error <- rbind(
  DT("CDPA", CDPA.error),
  DT("G.ig", PDPA.peaks.error),
  DT("G.jo", PDPA.infeasible.error[rule=="join"]),
  DT("G.rm", PDPA.infeasible.error[rule=="remove"]),
  DT("S.ig", Segmentor.peaks.error),
  DT("S.jo", Segmentor.infeasible.error[rule=="join"]),
  DT("S.rm", Segmentor.infeasible.error[rule=="rm"]))

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
  min.errors=min(total.errors),
  labels=labels[1]
), by=list(algo, chunk.name, sample.id)])

all.min[, algo.fac := factor(algo)]
all.min[, row.fac := factor(paste(chunk.name, sample.id))]
one.comparison <- all.min[c("G.rm", "G.jo"), on=list(algo)]
## from ?family
## For the ‘binomial’ and ‘quasibinomial’ families the response can
## be specified in one of three ways:

##   1. As a factor: ‘success’ is interpreted as the factor not
##      having the first level (and hence usually of having the
##      second level).

##   2. As a numerical vector with values between ‘0’ and ‘1’,
##      interpreted as the proportion of successful cases (with the
##      total number of cases given by the ‘weights’).

##   3. As a two-column integer matrix: the first column gives the
##      number of successes and the second the number of failures.

## fit <- glm(
##   min.errors/labels ~ algo.fac + row.fac,
##   "binomial", one.comparison, labels)

## biglm::bigglm

all.min[, table(chunk.name, algo)]

labels.check <- all.min[, list(
  min.labels=min(labels),
  max.labels=max(labels)
), by=list(chunk.name, sample.id)]
labels.check[min.labels!=max.labels]

min.wide <- dcast(
  all.min,
  chunk.name+sample.id+labels~algo,
  value.var="min.errors")

## as expected there are some (81) problems where the peak joining
## achieves fewer errors than ignoring the infeasible models with
## equality constraints.
min.wide[G.jo < G.ig]
min.wide[G.ig < G.jo]

## 18 problems where join<rm, 51 problems where rm<join.
min.wide[G.jo < G.rm]
min.wide[G.rm < G.jo]

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

## 15 problems where CDPA<rm, 48 problems where rm<CDPA.
min.wide[CDPA < G.rm]
min.wide[G.rm < CDPA]
min.wide[, table(CDPA-G.rm)]

## ignore vs join
g <- function(pos, neg, Comparison){
  x <- pos-neg
  dt <- data.table(table(x), Comparison)
  names(dt)[1] <- "diff"
  dt
}
count.tall <- min.wide[, rbind(
  g(G.rm, G.jo, "Remove-Join\nGPDPA"),
  ##g(S.rm-S.jo, "Remove-Join\nPDPA"),
  g(G.jo, G.ig, "Join-Ignore\nGPDPA"),
  ##g(S.jo-S.ig, "Join-Ignore\nPDPA"),
  g(G.rm, MACS, "GPDPAremove-MACS\n"),
  g(G.rm, HMCanBroad, "GPDPAremove-HMCanBroad\n"),
  g(G.rm, CDPA, "GPDPAremove-CDPA\n"),
  g(G.rm, S.rm, "GPDPA-PDPA\n(Remove rule)"),
  g(G.jo, S.jo, "GPDPA-PDPA\n(Join rule)"),
  g(G.ig, S.ig, "GPDPA-PDPA\n(Ignore rule)"))]
count.tall[, dnum := as.integer(diff)]
count.tall[, diff.fac := factor(diff, seq(min(dnum), max(dnum), by=1))]
stopifnot(!is.na(count.tall$diff.fac))
(count.wide <- dcast(count.tall, Comparison ~diff.fac, value.var="N"))

l <- sum(min.wide$labels)
two.labels <- c(l, l)
g <- function(pos, neg, Comparison){
  pos.total <- sum(pos)
  neg.total <- sum(neg)
  mean.est <- (pos.total+neg.total)/2
  diff.total <- pos.total-neg.total
  abs.total <- abs(diff.total)
  t.res <- prop.test(
    c(pos.total, neg.total),
    two.labels)
  p.value <- skellam::pskellam(-abs.total, mean.est, lower.tail=TRUE)*2
  w.res <- wilcox.test(
    pos,
    neg,
    paired=TRUE)
  data.table(
    Comparison, p.skellam=p.value, abs.total, mean.est,
    p.wilcox=w.res$p.value,
    p.prop=t.res$p.value)
}
(pval.dt <- min.wide[, rbind(
  g(G.rm, S.rm, "GPDPA-PDPA\n(Remove rule)"),
  g(G.jo, S.jo, "GPDPA-PDPA\n(Join rule)"),
  g(G.ig, S.ig, "GPDPA-PDPA\n(Ignore rule)"),
  g(G.rm, MACS, "GPDPAremove-MACS\n"),
  g(G.rm, HMCanBroad, "GPDPAremove-HMCanBroad\n"),
  g(G.rm, CDPA, "GPDPAremove-CDPA\n"),
  ##g(G.jo, CDPA, "GPDPAjoin-CDPA\n"),
  g(G.rm, G.jo, "Remove-Join\nGPDPA"),
  ##g(S.rm-S.jo, "Remove-Join\nPDPA"),
  g(G.jo, G.ig, "Join-Ignore\nGPDPA")
  )])
  ##g(S.jo-S.ig, "Join-Ignore\nPDPA"),
  

pval.dt[order(p.wilcox), list(Comparison, p.wilcox, p.prop, p.skellam)]
pval.dt[, list(Comparison, p.wilcox, p.prop, p.skellam)]

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
  first %in% c("Join-Ignore", "Remove-Join"),
  last,
  Comparison)]
count.tall[, panel := ifelse(
  first %in% c("Join-Ignore", "Remove-Join"),
  "Comparing
post-processing
rules", ifelse(
  first=="GPDPA-PDPA",
  "Comparing
constrained/
unconstrained", "Comparing
GPDPAremove
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
    size=3,
    data=count.tall)+
  scale_fill_gradient(
    "log10(
segmentation
problems)",
low="white", high="red")+
  scale_x_continuous(
    "Difference in minimum number of incorrect labels per segmentation problem",
    limits=c(-11.75, 6.5),
    breaks=unique(count.tall$diff.int))+
  scale_y_discrete(
    "Comparison",
    labels=function(x){
      x <- sub("\n$", "", x)
      ifelse(grepl("\n", x), x, sub("-", "\n-", x))
    }
  )+
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
pdf("figure-PDPA-infeasible-error-compare.pdf", 9, 4.2)
print(gg)
dev.off()

