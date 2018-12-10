source("packages.R")

## In these data the oracle penalty never selects more than 5
## segments. Is that normal? The linear penalty selects up to 19
## segments.

chunk.name <- "H3K36me3_TDH_immune/1"
chunk.dir <- file.path(
  "../chip-seq-paper/chunks",
  chunk.name)
Segmentor.model.RData <- file.path(chunk.dir, "Segmentor.model.RData")
load(Segmentor.model.RData)
sample.id <- "McGill0004"
model.list <- Segmentor.model[[sample.id]]
counts.RData <- file.path(chunk.dir, "counts.RData")
load(counts.RData)
sample.counts <- data.table(counts)[sample.id, on=list(sample.id)]
##uncompressed.vec <- sample.counts[, rep(coverage, chromEnd-chromStart)]
##fit <- Segmentor(uncompressed.vec, Kmax=19, 
##Cr <- SelectModel(model.list, penalty="oracle", keep=TRUE)
weight.vec <- sample.counts[, chromEnd-chromStart]
total.bases <- sum(weight.vec)
identical(weight.vec, model.list$DataComp)

alice.oracle <- function
### An adapted version of Alice's Segmentor3IsBack::SelectModel code
### for penalty="oracle" -- this function can be used with any model
### (not just Segmentor S4 classes).
(n,
### number of base pairs to segment.
 Lik
### Numeric vector of Kmax -sum(dpois(x, mu, log=TRUE)) values.
 ){
  Kmax <- length(Lik)
  saut <- function(Lv, pen, Kseq, seuil = sqrt(n)/log(n), biggest = TRUE) {
    J = -Lv
    Kmax = length(J)
    k = 1
    kv = c()
    dv = c()
    pv = c()
    dmax = 1
    while(k < Kmax) {
      pk = (J[(k + 1):Kmax] - J[k])/(pen[k] - pen[(k + 1):Kmax])
      pm = max(pk)
      dm = which.max(pk)
      dv = c(dv, dm)
      kv = c(kv, k)
      pv = c(pv, pm)
      if (dm > dmax) {
        dmax = dm
        kmax = k
        pmax = pm
      }
      k = k + dm
    }
    if(biggest) {
      pv = c(pv, 0)
      kv = c(kv, Kmax)
      dv = diff(kv)
      dmax = max(dv)
      rt = max(dv)
      rt = which(dv == rt)
      pmax = pv[rt[length(rt)]]
      alpha = 2 * pmax
      km = kv[alpha >= pv]
      Kh = Kseq[km[1]]
      return(c(Kh, alpha))
    }
    else {
      paux <- pv[which(kv <= seuil)]
      alpha <- 2 * min(paux)
      km = kv[alpha >= pv]
      Kh = Kseq[km[1]]
      return(c(Kh, alpha))
    }
  }
  ## oracle penalty code.
  Kseq <- 1:Kmax
  in.square <- 1 + 4 * sqrt(1.1 + log(n/Kseq))
  pen <- Kseq * in.square * in.square
  from.saut <- saut(-Lik[Kseq], pen, Kseq, 
                    n/log(n), biggest = FALSE)
  list(crit=Lik + from.saut[2] * pen,
       segments=from.saut[1])
}

my.oracle <- function
### An adapted version of Alice's Segmentor3IsBack::SelectModel code
### for penalty="oracle" -- this function can be used with any model
### (not just Segmentor S4 classes).
(n,
### number of base pairs to segment.
 Lik,
### Numeric vector of Kmax -sum(dpois(x, mu, log=TRUE)) values.
 beta=10^seq(-2, 4, l=200)
### beta from equation (6) in arXiv:1301.2534.
 ){
  Kmax <- length(Lik)
  Kseq <- 1:Kmax
  in.square <- 1 + 4 * sqrt(1.1 + log(n/Kseq))
  pen <- Kseq * in.square * in.square
  sapply(beta, function(b){
    which.min(Lik + b * pen)
  })
}

alice.oracle(total.bases, model.list$Likelihood)

my.oracle(total.bases, model.list$Likelihood)

load("all.modelSelection.RData")
load("unsupervised.Segmentor.RData")

all.modelSelection[chunk.name== & segments==19 & algo=="PDPA", list(sample.id, min.lambda, max.lambda)]

oracle.Segmentor[["H3K36me3_TDH_immune/1"]]
