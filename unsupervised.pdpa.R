source("packages.R")

set.seed(1)
N <- 250 
x <- rpois(10*N, rep(c(8,1,5,3,16,33,2,12,7,1),each=N))
Kmax <- 40
res <- Segmentor(data=x, model=1, Kmax=Kmax)
Cr <- SelectModel(res, penalty='oracle', keep=TRUE)
Cr.mBIC <- SelectModel(res, penalty='mBIC', keep=TRUE)
plot(Cr$criterion, type="l")
best.k <- which.min(Cr$criterion)
points(best.k, Cr$criterion[best.k])
mu <- mean(x)
my.lik.1seg <- -sum(dpois(x, mu, log=TRUE))
stopifnot(all.equal(res@likelihood[1], my.lik.1seg))
## Segmentor likelihood computation is the same as dpois!
  
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

### From cghseg.
getmBIC <- function(K,lv,mu,CGHo){  
  M   = length(names(mu))
  N   = M*mu[[1]]$end[dim(mu[[1]])[1]]
  Ent = sum(unlist(lapply(mu,FUN = function(x){log(x$end-x$begin+1)})))

  if (CGHo["calling"]==FALSE){
    mBIC = ((N-K+1)/2)*(lv*(2/N)+1+log(2*pi))-0.5*Ent -(K-M)*log(N)+lgamma((N-K+1)/2)-((N-K+1)/2)*log(N)
  } else {
    P  = CGHo["nblevels"]
    Np = sapply(mu,FUN=function(x){
	tabulate(rep(x$levels,x$end-x$begin+1),P)
    })
    Np = apply(Np,1,sum)
    Ent  = sum(log(Np))
    mBIC = ((N-P+1)/2)*(lv*(2/N)+1+log(2*pi))-0.5*Ent-(K-M)*log(N)+lgamma((N-P+1)/2)-((N-P+1)/2)*log(N)
  }
  return(mBIC)
}

alice.mBIC <- function
### An adapted version of Alice's Segmentor3IsBack::SelectModel code
### for penalty="mBIC" -- this function can be used with any model
### (not just Segmentor S4 classes).
(end.mat,
### Kmax x Kmax lower diagonal matrix. Row K has the last base indices
### of the model with K segments.
 Lik
### Numeric vector of Kmax -sum(dpois(x, mu, log=TRUE)) values.
 ){
  stopifnot(length(Lik) == nrow(end.mat))
  stopifnot(length(Lik) == ncol(end.mat))
  n <- end.mat[1, 1]
  sizenr <- function(k) {
    ## This is copied directly from Alice's code and I suspect it is
    ## incorrect -- the 1 should be 0.
    bases.vec <- diff(c(1, end.mat[k, 1:k]))
    sum(log(bases.vec))
    ## the number of base pairs is used in the penalty computation
    ## ... this will change with our weighted problem!
  }
  ## mBIC penalty code.
  entropy.term <- sapply(1:Kmax, sizenr)
  crit <- Lik + 0.5 * entropy.term + (1:Kmax - 0.5) * log(n)
  K <- which.min(crit)
  list(crit=crit,
       segments=K)
}

incorrect.mBIC <- function(end.mat, Lik, weight){
  stopifnot(length(Lik) == nrow(end.mat))
  stopifnot(length(Lik) == ncol(end.mat))
  n <- sum(weight)
  sizenr <- function(k) {
    model.ends <- end.mat[k, 1:k]
    model.starts <- c(1, 1+model.ends[-length(model.ends)])
    bases.vec <-
      sapply(seq_along(model.starts), function(segment.i){
        from <- model.starts[[segment.i]]
        to <- model.ends[[segment.i]]
        sum(weight[from:to])
      })
    bases.vec[1] <- bases.vec[1] - 1 # THIS LINE IS THE ONLY DIFFERENCE!
    sum(log(bases.vec))
  }
  ## mBIC penalty code.
  entropy.term <- sapply(1:Kmax, sizenr)
  crit <- Lik + 0.5 * entropy.term + (1:Kmax - 0.5) * log(n)
  K <- which.min(crit)
  list(crit=crit,
       segments=K)
}

correct.mBIC <- function(end.mat, Lik, weight){
  stopifnot(length(Lik) == nrow(end.mat))
  stopifnot(length(Lik) == ncol(end.mat))
  n <- sum(weight)
  Kmax <- length(Lik)
  sizenr <- function(k) {
    model.ends <- end.mat[k, 1:k]
    if(all(!is.na(model.ends))){
      model.starts <- c(1, 1+model.ends[-length(model.ends)])
      bases.vec <-
        sapply(seq_along(model.starts), function(segment.i){
          from <- model.starts[[segment.i]]
          to <- model.ends[[segment.i]]
          sum(weight[from:to])
        })
      sum(log(bases.vec))
    }else{
      NA
    }
  }
  ## mBIC penalty code.
  entropy.term <- sapply(1:Kmax, sizenr)
  crit <- Lik + 0.5 * entropy.term + (1:Kmax - 0.5) * log(n)
  K <- which.min(crit)
  list(crit=crit,
       segments=K)
}

## My copy of Alice's oracle criterion computation is the same as
## Alice's code in Segmentor3IsBack.
crit.info <- alice.oracle(length(x), res@likelihood)
stopifnot(all.equal(as.numeric(crit.info$crit), as.numeric(Cr$criterion)))

## My copy of Alice's mBIC criterion computation is the same as
## Alice's code in Segmentor3IsBack.
mBIC.info <- alice.mBIC(res@breaks, res@likelihood)
stopifnot(all.equal(as.numeric(mBIC.info$crit), as.numeric(Cr.mBIC$criterion)))

PoissonLik <- function(count, bases, end.mat){
  Kmax <- nrow(end.mat)
  lik <- rep(NA, Kmax)
  loss <- rep(NA, Kmax)
  for(segments in 1:Kmax){
    seg.lik <- rep(NA, segments)
    seg.loss <- rep(NA, segments)
    ends <- end.mat[segments, 1:segments]
    if(all(!is.na(ends))){
      breaks <- ends[-length(ends)]
      starts <- c(1, breaks+1)
      for(segment.i in 1:segments){
        first <- starts[segment.i]
        last <- ends[segment.i]
        seg.data <- count[first:last]
        seg.bases <- bases[first:last]
        seg.mean <- sum(seg.data * seg.bases)/sum(seg.bases)
        loglik.vec <- dpois(seg.data, seg.mean, log=TRUE)
        seg.lik[segment.i] <- -sum(loglik.vec * seg.bases)
        seg.loss[segment.i] <- PoissonLoss(seg.data, seg.mean, seg.bases)
      }
      lik[segments] <- sum(seg.lik)
      loss[segments] <- sum(seg.loss)
    }
  }
  attr(lik, "loss") <- loss
  lik
}

unsupervised <- list()
oracle.segments <- list()
model.files <- Sys.glob("data/H*/*/PDPA.model.RData")
i.vec <- seq_along(model.files)
i.vec <- 1:2
for(model.file.i in i.vec){
  model.file <- model.files[[model.file.i]]
  chunk.name <- sub("data/", "", dirname(model.file))
  cat(sprintf("%4d / %4d %s\n", model.file.i, length(model.files), chunk.name))
  load(model.file)
  sample.list <- split(counts, counts$sample.id)
  segments.list <- list()
  oseg.list <- list()
  for(sample.id in names(sample.list)){
    fit <- PDPA.model[[sample.id]]
    lik.vec <- rep(NA, 19L)
    force.na <- rep(FALSE, 19L)
    for(n.segments in seq_along(lik.vec)){
      seg.mean.vec <- fit$mean.mat[n.segments, 1:n.segments]
      ends.vec <- c(fit$ends.mat[n.segments, 1:n.segments], fit$n.data)
      is.feasible <- all(diff(seg.mean.vec) != 0)
      lik.vec[[n.segments]] <- if(is.feasible && n.segments %% 2){
        data.mean.vec <- rep(seg.mean.vec, diff(ends.vec))
        data.lik.vec <- dpois(fit$count.vec, data.mean.vec, log=TRUE)
        last.lik <- -sum(data.lik.vec)
      }else{
        force.na[[n.segments]] <- TRUE
        last.lik
      }
    }
    n.bases <- sum(fit$weight.vec)
    sample.oracle <- alice.oracle(n.bases, lik.vec)
    oseg.list[[sample.id]] <- my.oracle(n.bases, lik.vec)
    segSeq <- seq(1, 19, by=2)
    write.na <- force.na[segSeq]
    penalty.mat <-  #peaks x penalties
      cbind(oracle=sample.oracle$crit[segSeq],
            none=lik.vec[segSeq])
    rownames(penalty.mat) <- segSeq
    penalty.mat[write.na, ] <- NA
    selected.rows <- apply(penalty.mat, 2, which.min)
    selected.segs <- segSeq[selected.rows]
    names(selected.segs) <- names(selected.rows)
    segments.list[[sample.id]] <- selected.segs
  }
  unsupervised[[chunk.name]] <- do.call(rbind, segments.list)
  oracle.segments[[chunk.name]] <- do.call(rbind, oseg.list)
}
test <- do.call(rbind, unsupervised)
stopifnot(!is.null(test))

save(unsupervised, oracle.segments, 
     file="unsupervised.pdpa.RData")
