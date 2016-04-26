## works_with_R("3.2.2",
##              "hadley/scales@2c3edf45de56d617444dc38e47e0404173817886",
##              "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
##              "tdhock/animint@b0a708d82b58a8ac266c70495500a1fb3dc2470b",
##              "tdhock/directlabels@7b4b08a5dd0ab86e0b90902b3a233903ddd42311",
##              data.table="1.9.6")
library(data.table)
library(animint)
library(directlabels)

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

####Function using Rigaill's pruning methods on optimal partitioning, for normal data with change in mean and variance of 1.5####
set.seed(124)
y<-c(rnorm(60,2,1.5),rnorm(60,0,1))
y <- c(1, 10, 14, 13)
mu.min <- min(y)
mu.max <- max(y)
pen <- log(length(y))
timestep <- 4
n<-length(y)
F<-c()
F[1] <- -pen
CP <- list()
CP[[1]] <-
  quad<-function(x,A,B,C){
    return(A*x^2+B*x+C)
  }
LOC<-c()
LOC[1]<-0
Set<-list()
D<-c(min(y),max(y))
Set[[1]]<-D
a<-c(0);b<-c(0);c<-c(F[1]+pen)
plot(0)
not.bold.lines.list <- list()
bold.lines.list <- list()
intervals.list <- list()
Step <- function(s){
  factor(s, c("before", "unpruned", "pruned"))
}
plotFuns <- function(oloc.vec, step, timestep){
  for(kk in oloc.vec){
    if(kk==oloc.vec[1]){
      plot(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],type="l",col=co[kk==oloc.vec],ylim=F[taustar+1]+c(0,20),xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
    }else{
      lines(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],col=co[kk==oloc.vec])
    }
    not.bold.lines.list[[paste(timestep, step, kk)]] <<- data.table(
      timestep,
      step=Step(step),
      mean=x,
      cost=a[kk+1]*x^2+b[kk+1]*x+c[kk+1],
      kk,
      kk.fac=factor(kk),
      is.min=FALSE)
    interval.vec <- Set[[kk+1]]
    if(length(interval.vec)>0){
      mm=length(interval.vec)/2
      if(is.na(interval.vec[1])==FALSE){
        if(interval.vec[1]<mu.min){interval.vec[1]=mu.min}
        if(interval.vec[length(interval.vec)]>mu.max){
          interval.vec[length(interval.vec)]=mu.max
        }
      }
      delta=.4 #interval boundary segment height
      delta2=.005 #interval boundary segment offset
      if(!is.na(interval.vec[1])){
        start.i <- seq(1, length(interval.vec), by=2)
        end.i <- seq(2, length(interval.vec), by=2)
        intervals.list[[paste(timestep, step, kk)]] <<- data.table(
          timestep,
          step=Step(step),
          min=interval.vec[start.i],
          max=interval.vec[end.i],
          kk, kk.fac=factor(kk))
      }
      for(ii in 1:mm) {
        lines(interval.vec[2*ii-c(1,0)]+c(delta2,-delta2),rep(min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1])),2),lwd=3,col=co[(kk)==oloc.vec])#interval
        lines(rep(interval.vec[2*ii-1],2)+delta2,c(min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1]))-delta,min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1]))+delta),lwd=3,col=co[(kk)==oloc.vec])#left interval boundary
        lines(rep(interval.vec[2*ii-0],2)-delta2,c(min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1]))-delta,min(c[oloc.vec+1]-0.25*b[oloc.vec+1]^2/(1e-3+a[oloc.vec+1]))+delta),lwd=3,col=co[(kk)==oloc.vec])#right interval boundary
        ######boldbit
        if(is.na(interval.vec[2*ii-1])==FALSE){
          x2<-seq(interval.vec[2*ii-1],interval.vec[2*ii-0],length=50)
          lines(x2,a[kk+1]*x2^2+b[kk+1]*x2+c[kk+1],col=co[(kk)==oloc.vec],lwd=4)#bold part of function which attains the minimum.
          bold.lines.list[[paste(timestep, step, kk, ii)]] <<- data.table(
            timestep,
            step=Step(step),
            mean=x2,
            cost=a[kk+1]*x2^2 + b[kk+1]*x2 + c[kk+1],
            is.min=TRUE,
            kk, kk.fac=factor(kk),
            ii)
        }#has an interval
      }#for(ii
    }#if(length(interval.vec
  }##for(kk in oloc.vec
  legend("topleft", legend=paste(oloc.vec," "), col=co,pch=15,bg="white")
}
for (taustar in 1:n){
  temp<-c()
  for (tau in 1:length(LOC)){
    vs<-LOC[tau]
    a[vs+1]<-a[vs+1]+1/3
    b[vs+1]<-b[vs+1]-(2/3)*y[taustar]
    c[vs+1]<-c[vs+1]+(1/2)*log(3*pi)+(1/3)*(y[taustar])^2
    temp[tau]<-NA
    if(a[vs+1]==0){
      temp[tau]<-c[vs+1]
    }else{
      for(i in 1:(length(Set[[vs+1]])/2)){
        if(Set[[vs+1]][2*i-1] <= -b[vs+1]/(2*a[vs+1])
           & -b[vs+1]/(2*a[vs+1]) <= Set[[vs+1]][2*i]){
          temp[tau]<-quad(-b[vs+1]/(2*a[vs+1]),a[vs+1],b[vs+1],c[vs+1])}
      }
      if(is.na(temp[tau])){
        temp[tau] <- min(quad(Set[[vs+1]], a[vs+1], b[vs+1], c[vs+1]))
      }
    }
  }
  F[taustar+1]<-min(temp,na.rm=T)
  taudash<-LOC[(which(temp==F[taustar+1]))]
  CP[[taustar+1]]<-c(CP[[taudash+1]],taudash)
  a[taustar+1]<-0;b[taustar+1]<-0;c[taustar+1]<-F[taustar+1]+pen
  Set[[taustar+1]]<-D
  OLOC=c(LOC,taustar)
  ####THE PRUNING BIT####
  for(tau in LOC){
    I <- if((b[tau+1]^2-4*a[tau+1]*(c[tau+1]-F[taustar+1]-pen))<0){
      NA
    }else{
      c(
      (-b[tau+1]-sqrt(b[tau+1]^2-4*a[tau+1]*(c[tau+1]-F[taustar+1]-pen)))/
        (2*a[tau+1]),
      (-b[tau+1]+sqrt(b[tau+1]^2-4*a[tau+1]*(c[tau+1]-F[taustar+1]-pen)))/
        (2*a[tau+1]))
    }
    Set[[tau+1]]<-in1(Set[[tau+1]],I) #intersect function for continuous sets
    if (is.na(Set[[tau+1]][1])){
      LOC<-setdiff(LOC,tau) #discrete set uses the setdiff function
    }
    Set[[taustar+1]]<-setdiff1(Set[[taustar+1]],I) #continuous set uses my setdiff1 function
  }
  if(!is.na(Set[[taustar+1]][1])){
    LOC<-c(LOC,taustar)
  }
  ####PRUNING BIT ENDS####
  ######PLOT timestep-1 ###########
  ##if(taustar==timestep-1){
  x=seq(mu.min,mu.max,length=100)
  co=1:length(OLOC)
  OLOCdash<-c()
  for(it in OLOC){
    if(is.na(Set[[it+1]][1])==FALSE){
      OLOCdash<-c(OLOCdash,it)
    }
  }
  plotFuns(OLOCdash, "before", taustar+1)
  ## for(kk in OLOCdash){
  ##   if(kk==OLOCdash[1]){
  ##     plot(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],type="l",col=co[kk==OLOCdash],ylim=F[taustar+1]+c(0,20),xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
  ##   }else{
  ##     lines(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],col=co[kk==OLOCdash])
  ##   }
  ##   not.bold.lines.list[[paste("before", kk)]] <- data.table(
  ##     step=Step("before"),
  ##     mean=x,
  ##     cost=a[kk+1]*x^2+b[kk+1]*x+c[kk+1],
  ##     kk,
  ##     kk.fac=factor(kk),
  ##     is.min=FALSE)
  ##   interval.vec <- Set[[kk+1]]
  ##   if(length(interval.vec)>0){
  ##     mm=length(interval.vec)/2
  ##     if(is.na(interval.vec[1])==FALSE){
  ##       if(interval.vec[1]<mu.min){interval.vec[1]=mu.min}
  ##       if(interval.vec[length(interval.vec)]>mu.max){
  ##         interval.vec[length(interval.vec)]=mu.max
  ##       }
  ##     }
  ##     delta=.4 #interval boundary segment height
  ##     delta2=.005 #interval boundary segment offset
  ##     if(!is.na(interval.vec[1])){
  ##       start.i <- seq(1, length(interval.vec), by=2)
  ##       end.i <- seq(2, length(interval.vec), by=2)
  ##       intervals.list[[paste("before", kk)]] <- data.table(
  ##         step=Step("before"),
  ##         min=interval.vec[start.i],
  ##         max=interval.vec[end.i],
  ##         kk, kk.fac=factor(kk))
  ##     }
  ##     for(ii in 1:mm) {
  ##       lines(interval.vec[2*ii-c(1,0)]+c(delta2,-delta2),rep(min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1])),2),lwd=3,col=co[(kk)==OLOCdash])#interval
  ##       lines(rep(interval.vec[2*ii-1],2)+delta2,c(min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1]))-delta,min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1]))+delta),lwd=3,col=co[(kk)==OLOCdash])#left interval boundary
  ##       lines(rep(interval.vec[2*ii-0],2)-delta2,c(min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1]))-delta,min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1]))+delta),lwd=3,col=co[(kk)==OLOCdash])#right interval boundary
  ##       ######boldbit
  ##       if(is.na(interval.vec[2*ii-1])==FALSE){
  ##         x2<-seq(interval.vec[2*ii-1],interval.vec[2*ii-0],length=50)
  ##         lines(x2,a[kk+1]*x2^2+b[kk+1]*x2+c[kk+1],col=co[(kk)==OLOCdash],lwd=4)#bold part of function which attains the minimum.
  ##         bold.lines.list[[paste("before", kk, ii)]] <- data.table(
  ##           step=Step("before"),
  ##           mean=x2,
  ##           cost=a[kk+1]*x2^2 + b[kk+1]*x2 + c[kk+1],
  ##           is.min=TRUE,
  ##           kk, kk.fac=factor(kk),
  ##           ii)
  ##       }
  ##       ######
  ##     }#for(ii
  ##   }#length(Set[kk+1
  ## }##for(kk in OLOCdash
  ## legend("topleft", legend=paste(OLOCdash," "), col=co,pch=15,bg="white")}
  ####PLOT END########
  ##}
  ######PLOT timestep mid###########
  ##if(taustar==timestep){
  x=seq(mu.min,mu.max,length=100)
  OLOCdash<-c(OLOCdash,taustar)
  co=1:length(OLOCdash)
  plotFuns(OLOC, "unpruned", taustar)
  ## for(kk in OLOC){
  ##   if(kk==OLOC[1]){
  ##     plot(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],type="l",col=co[kk==OLOC],ylim=F[taustar+1]+c(0,20),xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
  ##   }else{
  ##     lines(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],col=co[kk==OLOC])
  ##   }
  ##   not.bold.lines.list[[paste("unpruned", kk)]] <- data.table(
  ##     step=Step("unpruned"),
  ##     mean=x,
  ##     cost=a[kk+1]*x^2+b[kk+1]*x+c[kk+1],
  ##     kk,
  ##     kk.fac=factor(kk),
  ##     is.min=FALSE)
  ##   ## Set is a list of numeric vectors, one for each cost
  ##   ## function. if the vector has only one element which is
  ##   ## missing, then it has no intervals, and it will be
  ##   ## pruned. Otherwise each vector should be an even number of
  ##   ## elements, which indicate the starts and ends of the
  ##   ## intervals, e.g. c(start1,end1,start2,end2)
  ##   interval.vec <- Set[[kk+1]]
  ##   if(length(interval.vec)>0){
  ##     mm <- length(interval.vec)/2
  ##     if(is.na(interval.vec[1])==FALSE){
  ##       if(interval.vec[1] < mu.min){
  ##         interval.vec[1] <- mu.min
  ##       }
  ##       if(interval.vec[length(interval.vec)] > mu.max){
  ##         interval.vec[length(interval.vec)] <- mu.max
  ##       }
  ##     }
  ##     delta <- .4 #size of "fleck"
  ##     delta2 <- .005 #size of space between "flecks"
  ##     if(!is.na(interval.vec[1])){
  ##       start.i <- seq(1, length(interval.vec), by=2)
  ##       end.i <- seq(2, length(interval.vec), by=2)
  ##       intervals.list[[paste("unpruned", kk)]] <- data.table(
  ##         step=Step("unpruned"),
  ##         min=interval.vec[start.i],
  ##         max=interval.vec[end.i],
  ##         kk, kk.fac=factor(kk))
  ##     }
  ##     for(ii in 1:mm) {
  ##       lines(interval.vec[2*ii-c(1,0)]+c(delta2,-delta2),rep(min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1])),2),lwd=3,col=co[kk==OLOC])
  ##       lines(rep(interval.vec[2*ii-1],2)+delta2,c(min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1]))-delta,min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1]))+delta),lwd=3,col=co[kk==OLOC])
  ##       lines(rep(interval.vec[2*ii-0],2)-delta2,c(min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1]))-delta,min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1]))+delta),lwd=3,col=co[kk==OLOC])
  ##       ######boldbit
  ##       if(is.na(interval.vec[2*ii-1])==FALSE){
  ##         ## It is essential that x2 starts and ends at the min and
  ##         ## max of the interval (otherwise the min envelope will
  ##         ## not be completely bold).
  ##         x2 <- seq(interval.vec[2*ii-1], interval.vec[2*ii-0], length=50)
  ##         lines(x2,a[kk+1]*x2^2+b[kk+1]*x2+c[kk+1],col=co[kk==OLOC],lwd=4)
  ##         bold.lines.list[[paste("unpruned", kk, ii)]] <- data.table(
  ##           step=Step("unpruned"),
  ##           mean=x2,
  ##           cost=a[kk+1]*x2^2 + b[kk+1]*x2 + c[kk+1],
  ##           is.min=TRUE,
  ##           kk, kk.fac=factor(kk),
  ##           ii)
  ##       }
  ##       ######
  ##     }
  ##   }
  ## }
  ## legend("topleft", legend=paste(OLOC," "), col=co[OLOC%in%OLOCdash],pch=15,bg="white")
  ####PLOT END########
  ######PLOT timestep end##########
  OLOCdash2<-c()
  for(it in OLOC){
    if(is.na(Set[[it+1]][1])==FALSE){
      OLOCdash2<-c(OLOCdash2,it)}
  }
  plotFuns(OLOCdash2, "pruned", taustar)
  ## for(kk in OLOCdash2){
  ##   if(kk==OLOCdash2[1]){
  ##     plot(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],type="l",col=co[kk==OLOCdash],ylim=F[taustar+1]+c(0,20),xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
  ##   }else{
  ##     lines(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],col=co[kk==OLOCdash])
  ##   }
  ##   not.bold.lines.list[[paste("pruned", kk)]] <- data.table(
  ##     step=Step("pruned"),
  ##     mean=x,
  ##     cost=a[kk+1]*x^2+b[kk+1]*x+c[kk+1],
  ##     kk,
  ##     kk.fac=factor(kk),
  ##     is.min=FALSE)
  ##   interval.vec <- Set[[kk+1]]
  ##   if(length(interval.vec)>0){
  ##     mm=length(interval.vec)/2
  ##     if(is.na(interval.vec[1])==FALSE){
  ##       if(interval.vec[1]<mu.min){interval.vec[1]=mu.min}
  ##       if(interval.vec[length(interval.vec)]>mu.max){interval.vec[length(interval.vec)]=mu.max}}
  ##     delta=.4 #size of "fleck"
  ##     delta2=.005 #size of space between "flecks"
  ##     if(!is.na(interval.vec[1])){
  ##       start.i <- seq(1, length(interval.vec), by=2)
  ##       end.i <- seq(2, length(interval.vec), by=2)
  ##       intervals.list[[paste("pruned", kk)]] <- data.table(
  ##         step=Step("pruned"),
  ##         min=interval.vec[start.i],
  ##         max=interval.vec[end.i],
  ##         kk, kk.fac=factor(kk))
  ##     }
  ##     for(ii in 1:mm) {lines(interval.vec[2*ii-c(1,0)]+c(delta2,-delta2),rep(min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1])),2),lwd=3,col=co[kk==OLOCdash])
  ##       lines(rep(interval.vec[2*ii-1],2)+delta2,c(min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1]))-delta,min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1]))+delta),lwd=3,col=co[kk==OLOCdash])
  ##       lines(rep(interval.vec[2*ii-0],2)-delta2,c(min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1]))-delta,min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1]))+delta),lwd=3,col=co[kk==OLOCdash])
  ##       ######boldbit
  ##       if(is.na(interval.vec[2*ii-1])==FALSE){
  ##         x2<-seq(interval.vec[2*ii-1],interval.vec[2*ii-0],length=50)
  ##         lines(x2,a[kk+1]*x2^2+b[kk+1]*x2+c[kk+1],col=co[kk==OLOCdash],lwd=4)}
  ##         bold.lines.list[[paste("pruned", kk, ii)]] <- data.table(
  ##           step=Step("pruned"),
  ##           mean=x2,
  ##           cost=a[kk+1]*x2^2 + b[kk+1]*x2 + c[kk+1],
  ##           is.min=TRUE,
  ##           kk, kk.fac=factor(kk),
  ##           ii)
  ##       ######
  ##     }
  ##   }
  ## }
  ## legend("topleft", legend=paste(OLOCdash2," "), col=co[OLOCdash%in%OLOCdash2],pch=15,bg="white")
  ####PLOT END########
  #print(taustar)
  ##}#if(taustar==timestep)
} 
#return(c(F[n+1],CP[[n+1]]))
intervals <- do.call(rbind, intervals.list)
not.bold.lines <- do.call(rbind, not.bold.lines.list)
bold.lines <- do.call(rbind, bold.lines.list)
cost.min <- not.bold.lines[, list(min=min(cost)), by=.(step, timestep, kk)]

##dput(RColorBrewer::brewer.pal(Inf, "Set1"))
change.colors <-
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
    "#A65628", "#F781BF", "#999999")
with.legend <- ggplot()+
  theme_bw()+
  coord_cartesian(ylim=c(0, max(cost.min$min)+1))+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(timestep ~ step)+
  scale_color_manual(values=change.colors)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=0.5))+
  geom_line(aes(mean, cost, color=kk.fac, size=is.min,
                group=kk),
            data=not.bold.lines)+
  geom_segment(aes(min, -Inf, xend=max, yend=-Inf, color=kk.fac),
               data=intervals,
               size=2)+
  geom_text(aes(min, -Inf, color=kk.fac, label="|"),
            data=intervals)+
  geom_line(aes(mean, cost, color=kk.fac, size=is.min,
                group=paste(kk, ii)),
            data=bold.lines)

with.labels <- direct.label(with.legend, "last.qp")

pdf("figure-unconstrained-FPOP-normal.pdf")
print(with.labels)
dev.off()

####Code to run the function####
#y<-SimChange(1000,c(1,2,3,0,5,4,6,5,1),1.5)
#plot(y[[2]],typ="l")
#abline(v=y[[1]],col="red")

#OPR<-OptimalPartitioningRig(y[[2]],negloglike.norm.mean2,log(1000));OPR
#abline(v=OPR,col="blue")



