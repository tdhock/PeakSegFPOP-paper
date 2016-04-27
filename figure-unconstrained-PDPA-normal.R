source("packages.R")
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
show.segments.list <- list()
Step <- function(s){
  factor(s, c("before", "unpruned", "pruned"))
}
plotFuns <- function(oloc.vec, step, timestep, n.segs){
  x <- seq(mu.min, mu.max, length=100)
  co <- seq_along(oloc.vec)
  for(kk in oloc.vec){
    if(kk==oloc.vec[1]){
      plot(x,a[kk]*x^2+b[kk]*x+c[kk],type="l",col=co[kk==oloc.vec],
           xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
    }else{
      lines(x,a[kk]*x^2+b[kk]*x+c[kk],col=co[kk==oloc.vec])
    }
    interval.vec <- Set[[kk]]
    n.intervals <- length(interval.vec)/2
    has.intervals <- !is.na(interval.vec[1])
    if(a[kk] != 0){
      mean.at.minimum <- -b[kk]/(2*a[kk])
      cost.at.minimum <- quad(mean.at.minimum, a[kk], b[kk], c[kk])
      previous.segments <- optimal.segments.list[[paste(m-1, kk)]]
      segment.start <- kk+1
      new.segment <- data.table(
        segment.start,
        segment.end=timestep,
        mean=mean.at.minimum)
      show.segments.list[[paste(timestep, step, kk, n.segs)]] <<- data.table(
        has.intervals,
        timestep,
        step,
        kk, previous.segment.end=factor(kk),
        segments=n.segs,
        rbind(previous.segments, new.segment))
      minima.list[[paste(timestep, step, kk, n.segs)]] <<- data.table(
        has.intervals,
        timestep,
        step,
        kk, previous.segment.end=factor(kk, 1:n),
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
      previous.segment.end=factor(kk, 1:n),
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
        kk, previous.segment.end=factor(kk, 1:n))
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
          has.intervals=TRUE,
          segments=n.segs,
          step=Step(step),
          mean=x2,
          cost=a[kk]*x2^2 + b[kk]*x2 + c[kk],
          is.min=TRUE,
          kk, previous.segment.end=factor(kk, 1:n),
          ii)
      }#has an interval
    }#for(ii
  }##for(kk in oloc.vec
  legend("topleft", legend=paste(oloc.vec," "), col=co,pch=15,bg="white")
}
y <- c(rnorm(60,2,1.5),rnorm(60,0,1))
y <- c(1, 10, 14, 13, 5, 5, 10, 3, 10)
y <- c(1, 10, 14, 13)
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
n <- length(y)

optimal.segments.list <- list()
for(last.segment.end in seq_along(y)){
  optimal.segments.list[[paste(1, last.segment.end)]] <- data.table(
    segment.start=1,
    segment.end=last.segment.end,
    mean=mean.vec[[last.segment.end]])
}

quad<-function(x,A,B,C){
  return(A*x^2+B*x+C)
}
tau<-matrix(nrow=maxseg,ncol=n)
output<-matrix(nrow=maxseg,ncol=maxseg+1)
output[1,1]<-cost.mat[1,n]
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
  c[m-1] <- cost.mat[m-1, m-1]  #min cost of m-1 CPs, last one at m-1
  Set[[m-1]] <- D
  for (j in m:n){
    a[j] <- 0
    b[j] <- 0
    c[j] <- cost.mat[m-1, j]  #min cost of m-1 CPs, last one at j
    Set[[j]] <- D
    temp <- c()
    OLOC <- c(LOC[[m]],j)
    for (v in LOC[[m]]){
      ## This loss function is the full sum of squares (including the
      ## constant data-squared term).
      a[v] <- a[v] + 1
      b[v] <- b[v] - 2*y[j]
      c[v] <- c[v] + y[j]^2
      discriminant <- b[v]^2 - 4*a[v]*(c[v]-cost.mat[m-1,j])
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
    cost.mat[m,j] <- optimal.cost <- min(min.vec, na.rm=T)
    tau[m,j] <- segment.before.end <- LOC[[m]][(which(min.vec==cost.mat[m,j]))]
    optimal.segments.before <-
      optimal.segments.list[[paste(m-1, segment.before.end)]]
    segment.start <- segment.before.end+1
    new.segment <- data.table(
      segment.start,
      segment.end=j,
      mean=mean(y[segment.start:j]))
    optimal.segments.list[[paste(m, j)]] <- rbind(
      optimal.segments.before,
      new.segment)
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
  output[m,1]<-cost.mat[m,n]
  output[m,2:(m+1)]<-taustar[-1]
}#for(m number of segments
minima <- do.call(rbind, minima.list)
intervals <- do.call(rbind, intervals.list)
not.bold.lines <- do.call(rbind, not.bold.lines.list)
bold.lines <- do.call(rbind, bold.lines.list)
cost.min <- not.bold.lines[, list(min=min(cost)), by=.(step, timestep, segments, kk)]
show.segments <- do.call(rbind, show.segments.list)

##dput(RColorBrewer::brewer.pal(Inf, "Set1"))
change.colors <- if(length(y)==4){
  c("1"="#E41A1C", #red
    "2"="#377EB8", #blue
    "#4DAF4A", #green
    "#984EA3", #purple
    "#FF7F00", #orange
    "#FFFF33", #yellow
    "3"="#A65628",
    "4"="#F781BF", "#999999")
}else{
  c("#E41A1C", #red
    "#377EB8", #blue
    "#4DAF4A", #green
    "#984EA3", #purple
    "#FF7F00", #orange
    "#FFFF33", #yellow
    "#A65628",
    "#F781BF", "#999999")
}
cost.limits <- max(cost.min$min)*c(-0.05, 1.05)
not.bold.unpruned <- not.bold.lines[step=="unpruned" & cost < cost.limits[2],]
bold.unpruned <- bold.lines[step=="unpruned",]
minima.unpruned <- minima[step=="unpruned",]
segments.unpruned <- show.segments[step=="unpruned",]
with.legend <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ timestep, labeller=label_both)+
  scale_color_manual(values=change.colors)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=0.25))+
  scale_linetype_manual(values=c("TRUE"="solid", "FALSE"="dashed"))+
  geom_line(aes(mean, cost, color=previous.segment.end, size=is.min,
                linetype=has.intervals,
                group=kk),
            data=not.bold.unpruned)+
  geom_line(aes(mean, cost, color=previous.segment.end, size=is.min,
                group=paste(kk, ii)),
            data=bold.unpruned)+
  geom_point(aes(mean, cost, color=previous.segment.end),
             shape=21,
             fill="white",
             data=minima.unpruned)+
  xlab("segment mean")+
  ylab("cost")
with.legend

data.dt <- data.table(
  count=y,
  kk=seq_along(y),
  timestep=seq_along(y),
  previous.segment.end=factor(seq_along(y)))
interval.counts <-
  intervals[step=="pruned", list(intervals=.N), by=.(segments, timestep)]
addBoth <- function(dt, x.var, y.var){
  dt <- data.table(
    dt,
    x.var=factor(x.var, c("data", "cost")))
  if(!is.null(y.var)){
    dt$y.var <- factor(y.var, c("count", "intervals", "segments"))
  }
  dt
}
addX <- function(dt, x.var="cost", y.var="count"){
  addBoth(dt, x.var, y.var)
}
addY <- function(dt, y.var){
  addBoth(dt, "data", y.var)
}
dimnames(cost.mat) <- dimnames(tau) <- list("segments"=NULL, "timestep"=NULL)
cost.rects <- data.table(melt(cost.mat, value.name="cost"))[!is.na(cost),]
segment.rects <- data.table(
  segments=1:maxseg)
tau.text <- data.table(melt(tau, value.name="last.change"))[!is.na(last.change),]
tau.text[, previous.segment.end := factor(last.change, 1:n)]
viz <- list(
  title=paste(
    "Pruned Dynamic Programming Algorithm",
    "for optimal multiple change-point detection"),
  data=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=1000, height=500)+
    facet_grid(y.var ~ x.var, scales="free")+
    ylab("")+
    xlab("data sequence")+
    scale_color_manual("data", values=change.colors)+
    geom_point(aes(timestep, intervals, showSelected=segments),
               data=addY(interval.counts, "intervals"))+
    geom_segment(aes(segment.start-0.5, mean,
                     xend=segment.end+0.5, yend=mean,
                     showSelected=timestep,
                     showSelected2=segments,
                     showSelected3=previous.segment.end),
                 color="green",
                 data=addY(segments.unpruned, "count"))+
    geom_tallrect(aes(xmin=kk-0.5, xmax=kk+0.5, clickSelects=timestep),
                  alpha=0.5,
                  data=addX(data.dt, "data", NULL))+
    geom_point(aes(kk, count),
               data=addY(data.dt, "count"))+
    geom_point(aes(kk, data.dt$count[kk],
                   color=previous.segment.end,
                   showSelected=segments,
                   showSelected2=timestep,
                   clickSelects=previous.segment.end),
               size=5,
               data=addY(minima, "count"))+
    geom_tile(aes(timestep, segments, fill=log(cost+1)),
              data=addY(cost.rects, "segments"))+
    geom_widerect(aes(ymin=segments-0.5, ymax=segments+0.5,
                      clickSelects=segments),
                  fill=NA,
                  alpha=0.5,
                  data=addY(segment.rects, "segments"))+
    geom_text(aes(timestep, segments,
                  label=previous.segment.end,
                  clickSelects=previous.segment.end,
                  color=previous.segment.end),
              data=addY(tau.text, "segments"))+
    scale_size_manual(values=c("TRUE"=3, "FALSE"=1))+
    scale_linetype_manual(values=c("TRUE"="solid", "FALSE"="dashed"))+
    guides(color="none")+
    geom_path(aes(cost, mean, color=previous.segment.end, size=is.min,
                  linetype=has.intervals,
                  clickSelects=previous.segment.end,
                  showSelected=timestep,
                  showSelected2=segments,
                  group=paste(kk)),
              data=addX(not.bold.unpruned))+
    geom_path(aes(cost, mean, color=previous.segment.end, size=is.min,
                  linetype=has.intervals,
                  clickSelects=previous.segment.end,
                  showSelected=timestep,
                  showSelected2=segments,
                  group=paste(kk, ii)),
              data=addX(bold.unpruned))+
    geom_point(aes(cost, mean,
                   clickSelects=previous.segment.end,
                   showSelected=timestep,
                   showSelected2=segments,
                   showSelected3=has.intervals,
                   color=previous.segment.end),
               shape=21,
               fill="white",
               data=addX(minima.unpruned))
  )
animint2dir(viz, "figure-unconstrained-PDPA-normal")

with.legend <- ggplot()+
  theme_bw()+
  coord_cartesian(ylim=cost.limits)+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments+timestep ~ step, labeller=label_both)+
  scale_color_manual(values=change.colors)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=0.25))+
  geom_segment(aes(min, -Inf, xend=max, yend=-Inf, color=previous.segment.end),
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
  geom_line(aes(mean, cost, color=previous.segment.end, size=is.min,
                group=kk),
            data=not.bold.lines)+
  geom_line(aes(mean, cost, color=previous.segment.end, size=is.min,
                group=paste(kk, ii)),
            data=bold.lines)+
  xlab("segment mean")+
  ylab("cost")
with.legend

pdf("figure-unconstrained-PDPA-normal.pdf")
print(with.legend)
dev.off()

