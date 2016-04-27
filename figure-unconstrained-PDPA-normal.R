library(data.table)
library(ggplot2)
####Code to run Rigaill's method####
#y<-SimChange(10000,c(1,.4,.8,2,3,1.6,.3,.4,.2,1,4,.2,.8,1.2,1.4),.5)
#plot(y[[2]],typ="l")
#abline(v=y[[1]],col="red")
#Q2<-Rigaill(y[[2]],14);Q2

#y=rnorm(100)
#Rigaill(y,3)

####Functions for doing set operations on continuous sets, where a set is stored as a vector of endpoints####


####Union####
u1 <- function(A,B) 
{ 
  if(is.na(A[1])){
    return(B)
  }
  else if(is.na(B[1])){
    return(A)
  }
  else{
    Both<-c(A,B)
    m<-matrix(Both,ncol=2,byrow=T)
    m2<-m[order(m[,1]),]
    Both<-as.vector(t(m2)) #vector now of ordered pairs
    U<-c()
    U[1:2]<-Both[1:2]
    for(i in 2:((length(Both))/2)){
      if(U[length(U)-1]<=Both[2*i-1]&Both[2*i]<=U[length(U)]){
        #do nothing#
      }
      else if(U[length(U)-1]<=Both[2*i-1]&Both[2*i-1]<=U[length(U)]&U[length(U)]<=Both[2*i]){
        U[(length(U)-1):length(U)]<-c(U[length(U)-1],Both[2*i])
      }
      else{U[(length(U)+1):(length(U)+2)]<-Both[(2*i-1):(2*i)]}
    }
    return(U)}
} 

####Intersection####
in1<-function(A,B){
  if(is.na(A[1])|is.na(B[1])){return(NA)}
  else{
    U<-c()
    for(i in 1:(length(A)/2)){
      for(j in 1:(length(B)/2)){
        i1<-max(A[2*i-1],B[2*j-1])
        i2<-min(A[2*i],B[2*j])
        if(i1<=i2){U[(length(U)+1):(length(U)+2)]<-c(i1,i2)}
      }
    }
    if(is.null(U)){return(NA)}
    I<-u1(U,U[1:2])
    return(I)}
}

####Set Difference####
setdiff1<-function(A,B,M=10000){
  if(is.na(A[1])){return(NA)}
  else if(is.na(B[1])){return(A)}
  else{
    delta<-1/M
    m<-matrix(B,ncol=2,byrow=T)
    m2<-m[order(m[,1]),]
    B<-as.vector(t(m2)) #vector now of ordered pairs
    index<-rep(c(-delta,delta),length=length(B))
    BC<-c(-M,B+index,M)
    se<-in1(A,BC)
    return(se)}
}

####Function for running PDPA on Normal data with a change in mean and variance 1.5####
not.bold.lines.list <- list()
bold.lines.list <- list()
intervals.list <- list()
minima.list <- list()
Step <- function(s){
  factor(s, c("before", "unpruned", "pruned"))
}
plotFuns <- function(oloc.vec, step, timestep, n.segs){
  x <- seq(mu.min, mu.max, length=100)
  co <- seq_along(oloc.vec)
  for(kk in oloc.vec){
    if(kk==oloc.vec[1]){
      plot(x,a[kk]*x^2+b[kk]*x+c[kk],type="l",col=co[kk==oloc.vec],xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
    }else{
      lines(x,a[kk]*x^2+b[kk]*x+c[kk],col=co[kk==oloc.vec])
    }
    interval.vec <- Set[[kk]]
    n.intervals <- length(interval.vec)/2
    has.intervals <- !is.na(interval.vec[1])
    if(a[kk] != 0){
      mean.at.minimum <- -b[kk]/(2*a[kk])
      cost.at.minimum <- quad(mean.at.minimum, a[kk], b[kk], c[kk])
      minima.list[[paste(timestep, step, kk, n.segs)]] <<- data.table(
        timestep,
        step,
        kk, kk.fac=factor(kk),
        segments=n.segs,
        mean=mean.at.minimum,
        cost=cost.at.minimum)
    }
    not.bold.lines.list[[paste(timestep, step, kk, n.segs)]] <<- data.table(
      has.intervals,
      n.intervals,
      timestep,
      segments=n.segs,
      step=Step(step),
      mean=x,
      cost=a[kk]*x^2+b[kk]*x+c[kk],
      kk,
      kk.fac=factor(kk),
      is.min=FALSE)
    delta=.4 #interval boundary segment height
    delta2=.005 #interval boundary segment offset
    if(has.intervals){
      if(interval.vec[1]<mu.min){
        interval.vec[1]=mu.min
      }
      if(interval.vec[length(interval.vec)]>mu.max){
        interval.vec[length(interval.vec)]=mu.max
      }
      start.i <- seq(1, length(interval.vec), by=2)
      end.i <- seq(2, length(interval.vec), by=2)
      intervals.list[[paste(timestep, step, kk, n.segs)]] <<- data.table(
        timestep,
        segments=n.segs,
        step=Step(step),
        min=interval.vec[start.i],
        max=interval.vec[end.i],
        kk, kk.fac=factor(kk))
    }
    for(ii in 1:n.intervals) {
      lines(interval.vec[2*ii-c(1,0)]+c(delta2,-delta2),rep(min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1])),2),lwd=3,col=co[(kk)==oloc.vec])#interval
      lines(rep(interval.vec[2*ii-1],2)+delta2,c(min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1]))-delta,min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1]))+delta),lwd=3,col=co[(kk)==oloc.vec])#left interval boundary
      lines(rep(interval.vec[2*ii-0],2)-delta2,c(min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1]))-delta,min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1]))+delta),lwd=3,col=co[(kk)==oloc.vec])#right interval boundary
      ######boldbit
      if(is.na(interval.vec[2*ii-1])==FALSE){
        x2<-seq(interval.vec[2*ii-1],interval.vec[2*ii-0],length=50)
        lines(x2,a[kk]*x2^2+b[kk]*x2+c[kk],col=co[(kk)==oloc.vec],lwd=4)#bold part of function which attains the minimum.
        bold.lines.list[[paste(timestep, step, kk, ii, n.segs)]] <<- data.table(
          timestep,
          segments=n.segs,
          step=Step(step),
          mean=x2,
          cost=a[kk]*x2^2 + b[kk]*x2 + c[kk],
          is.min=TRUE,
          kk, kk.fac=factor(kk),
          ii)
      }#has an interval
    }#for(ii
  }##for(kk in oloc.vec
  legend("topleft", legend=paste(oloc.vec," "), col=co,pch=15,bg="white")
}
y<-c(rnorm(60,2,1.5),rnorm(60,0,1))
y <- c(1, 10, 14, 13)
##Rigaill(y,2,45,2,0.5,4.5)
maxseg <- 3
mu.min <- min(y)
mu.max <- max(y)

## My code for checking.
cum.y <- cumsum(y)
t.seq <- seq_along(y)
mean.vec <- cum.y/t.seq
cost.mat <- matrix(NA, maxseg, length(y))
cum.y.sq <- cumsum(y*y)
cost.mat[1,] <- t.seq*mean.vec*mean.vec + cum.y.sq - 2*mean.vec*cum.y

timestep <- 3
seg <- 2
n<-length(y)
S1<-0
for(i in 1:n){
  S1[i+1]<-S1[i]+y[i]
}
SS1<-0
for(i in 1:n){
  SS1[i+1]<-SS1[i]+(y[i])^2
}
C<-matrix(nrow=maxseg,ncol=n)
for (t in 1:(n)){
  mu<-S1[t+1]/t
  C[1,t] <- SS1[t+1] - 2*mu*S1[t+1] + t*mu^2
}

quad<-function(x,A,B,C){
  return(A*x^2+B*x+C)
}
tau<-matrix(nrow=maxseg,ncol=n)
output<-matrix(nrow=maxseg,ncol=maxseg+1)
output[1,1]<-C[1,n]
LOC<-list()
D<-c(min(y),max(y))
for (m in 2:maxseg){
  Set<-list()
  LOC[[m]] <- m-1
  a <- c()
  b <- c()
  c <- c()
  a[m-1] <- 0
  b[m-1] <- 0
  c[m-1] <- C[m-1, m-1]  #min cost of m-1 CPs, last one at m-1
  Set[[m-1]] <- D
  for (j in m:n){
    a[j] <- 0
    b[j] <- 0
    c[j] <- C[m-1, j]  #min cost of m-1 CPs, last one at j
    Set[[j]] <- D
    temp <- c()
    OLOC <- c(LOC[[m]],j)
    for (v in LOC[[m]]){
      ## This loss function is the full sum of squares (including the
      ## constant data-squared term).
      a[v] <- a[v] + 1
      b[v] <- b[v] - 2*y[j]
      c[v] <- c[v] + y[j]^2
      discriminant <- b[v]^2 - 4*a[v]*(c[v]-C[m-1,j])
      I <- if(discriminant < 0){
        NA
      }else{
        (-b[v]+c(-1,1)*sqrt(discriminant))/(2*a[v])
      }
      Set[[v]]<-in1(Set[[v]],I) #intersect function for continuous sets
      if (is.na(Set[[v]][1])){
        LOC[[m]]<-setdiff(LOC[[m]],v) #discrete set uses the setdiff function
      }
      Set[[j]]<-setdiff1(Set[[j]],I) #continuous set uses my setdiff1 function
    }
    if(!is.na(Set[[j]][1])){
      LOC[[m]]<-c(LOC[[m]],j)
    }
    ###find minimum of cost function
    min.vec <- c()
    for (v in 1:length(LOC[[m]])){
      vs <- LOC[[m]][v]
      interval.vec <- Set[[vs]]
      min.vec[v] <- NA
      if(a[vs]==0){
        min.vec[v] <- c[vs]
      }else{
        mean.at.minimum <- -b[vs]/(2*a[vs])
        cost.at.minimum <- quad(mean.at.minimum, a[vs], b[vs], c[vs])
        for(i in 1:(length(interval.vec)/2)){
          interval.min <- interval.vec[2*i-1]
          interval.max <- interval.vec[2*i]
          if(interval.min <= mean.at.minimum & mean.at.minimum <= interval.max){
            min.vec[v] <- cost.at.minimum
          }
        }
        if(is.na(min.vec[v])){
          min.vec[v] <- min(quad(interval.vec, a[vs], b[vs], c[vs]))
        }
      }
    }
    C[m,j] <- optimal.cost <- min(min.vec, na.rm=T)
    tau[m,j] <- last.segment.end <- LOC[[m]][(which(min.vec==C[m,j]))]
    ## Plot all functions before pruning.
    plotFuns(OLOC, "unpruned", j, m)
    ## Plot functions after pruning.
    has.intervals <- sapply(OLOC, function(it)!is.na(Set[[it]][1]))
    OLOCdash <- OLOC[has.intervals]
    plotFuns(OLOCdash, "before", j+1, m)
    plotFuns(OLOCdash, "pruned", j, m)
  }#for(j data point
  em<-m
  taustar<-c()
  taustar[m+1]<-n
  while(em >1){
    taustar[em]<-tau[em,taustar[em+1]]
    em<-em-1
  }
  output[m,1]<-C[m,n]
  output[m,2:(m+1)]<-taustar[-1]
}#for(m number of segments
minima <- do.call(rbind, minima.list)
intervals <- do.call(rbind, intervals.list)
not.bold.lines <- do.call(rbind, not.bold.lines.list)
bold.lines <- do.call(rbind, bold.lines.list)
cost.min <- not.bold.lines[, list(min=min(cost)), by=.(step, timestep, segments, kk)]

##dput(RColorBrewer::brewer.pal(Inf, "Set1"))
change.colors <-
  c("1"="#E41A1C", #red
    "2"="#377EB8", #blue
    "#4DAF4A", #green
    "#984EA3", #purple
    "#FF7F00", #orange
    "#FFFF33", #yellow
    "3"="#A65628",
    "4"="#F781BF", "#999999")
cost.limits <- max(cost.min$min)*c(-0.05, 1.05)
not.bold.unpruned <- not.bold.lines[step=="unpruned" & cost < cost.limits[2],]
bold.unpruned <- bold.lines[step=="unpruned",]
minima.unpruned <- minima[step=="unpruned",]
with.legend <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ timestep, labeller=label_both)+
  scale_color_manual(values=change.colors)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=0.25))+
  scale_linetype_manual(values=c("TRUE"="solid", "FALSE"="dashed"))+
  geom_line(aes(mean, cost, color=kk.fac, size=is.min,
                linetype=has.intervals,
                group=kk),
            data=not.bold.unpruned)+
  geom_line(aes(mean, cost, color=kk.fac, size=is.min,
                group=paste(kk, ii)),
            data=bold.unpruned)+
  geom_point(aes(mean, cost, color=kk.fac),
             shape=21,
             fill="white",
             data=minima.unpruned)+
  xlab("segment mean")+
  ylab("cost")
with.legend

with.legend <- ggplot()+
  theme_bw()+
  coord_cartesian(ylim=cost.limits)+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments+timestep ~ step, labeller=label_both)+
  scale_color_manual(values=change.colors)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=0.25))+
  geom_segment(aes(min, -Inf, xend=max, yend=-Inf, color=kk.fac),
               data=intervals,
               size=2)+
  geom_text(aes(min, -Inf, label="|"),
            color="white",
            size=10,
            data=intervals)+
  geom_text(aes(min, -Inf, label="|"),
            color="black",
            size=5,
            data=intervals)+
  geom_line(aes(mean, cost, color=kk.fac, size=is.min,
                group=kk),
            data=not.bold.lines)+
  geom_line(aes(mean, cost, color=kk.fac, size=is.min,
                group=paste(kk, ii)),
            data=bold.lines)+
  xlab("segment mean")+
  ylab("cost")
with.legend

pdf("figure-unconstrained-PDPA-normal.pdf")
print(with.legend)
dev.off()

