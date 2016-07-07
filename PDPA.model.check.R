source("packages.R")

PDPA.RData.vec <- Sys.glob("data/*/*/PDPA.model.RData")

for(file.i in seq_along(PDPA.RData.vec)){
  PDPA.RData <- PDPA.RData.vec[[file.i]]
  load(PDPA.RData)
  chunk.dir <- dirname(PDPA.RData)
  for(compare.base in c("dp.model.RData", "dp.model.reverse.RData")){
    dp.RData <- file.path(chunk.dir, compare.base)
    if(file.exists(dp.RData)){
      objs <- load(dp.RData)
      for(sample.id in names(dp.model)){
        dp.sample <- dp.model[[sample.id]]
        pdpa.sample <- PDPA.model[[sample.id]]
        n.data <- ncol(pdpa.sample$cost.mat)
        loss.dt <- data.table(dp.sample$error)
        ## Why are some models missing??
        if(is.numeric(loss.dt$error)){
          loss.dt[, dp := error]
          loss.dt[, pdpa := pdpa.sample$cost.mat[segments, n.data] ]
          loss.dt[, should.be.positive := dp - pdpa]
          bad.loss <- loss.dt[should.be.positive < -1e-9, ]
          if(nrow(bad.loss)){
            cat(sprintf("%s %s %s\n", PDPA.RData, dp.RData, sample.id))
            print(bad.loss)
            load(file.path(chunk.dir, "dp.timing.RData"))
            print(dp.timing[sample.id,])
            load(file.path(chunk.dir, "counts.RData"))
            counts.by.sample <- split(counts, counts$sample.id)
            one.sample <- counts.by.sample[[sample.id]]
            dp.segments <-
              data.table(dp.sample$segments)[segments %in% bad.loss$segments,]
            pdpa.segments.list <- list()
            for(n.segments in bad.loss$segments){
              break.vec <- pdpa.sample$ends.mat[n.segments, 2:n.segments]
              first <- c(1, break.vec+1)
              last <- c(break.vec, n.data)
              pdpa.segments.list[[paste(n.segments)]] <- data.table(
                mean=pdpa.sample$mean.mat[n.segments, 1:n.segments],
                first,
                last,
                chromStart=one.sample$chromStart[first],
                chromEnd=one.sample$chromEnd[last],
                status=rep(c("background", "peak"), l=n.segments),
                peaks=(n.segments-1)/2,
                segments=n.segments)
            }
            pdpa.segments <- do.call(rbind, pdpa.segments.list)
            rbind(
              pdpa.last=pdpa.segments$last,
              dp.last=dp.segments$last)
            data.vec <- one.sample$coverage
            weight.vec <- with(one.sample, chromEnd-chromStart)
            pdpa <- PeakSegPDPA(data.vec, weight.vec, 19L)
            fit <- cDPA(data.vec, weight.vec, 19L)
            all.loss <- data.table(
              dp=as.numeric(fit$loss),
              pdpa=as.numeric(pdpa.sample$cost.mat),
              segments=as.integer(row(fit$loss)),
              data=as.integer(col(fit$loss)))
            all.loss[, should.be.positive := dp - pdpa]
            all.loss[should.be.positive < -1e-8,]
            stop("dp model more likely than pdpa model")
          }
        }
      }
    }
  }
}

library(data.table)
compare.cost <- data.table(read.table(header=TRUE,text="
    Linear        Log   Constant   min_mean   max_mean     data_i
         0          0           901.495779    0.0000000000000000000000000000000000000000000000000000000    1.6448087431693989568515235077938996255397796630859375000 25
       732      -1204           296.635210    1.6448087431693989568515235077938996255397796630859375000    1.6864173339916479488209688497590832412242889404296875000 25
         0          0           901.874643    1.6864173339916479488209688497590832412242889404296875000    1.7140695915279879901760295979329384863376617431640625000 25
       661      -1133           379.414830    1.7140695915279879901760295979329384863376617431640625000    2.7167911989387127924544529378181323409080505371093750000 25
        91       -182           977.507459    2.7167911989387127924544529378181323409080505371093750000   33.0000000000000000000000000000000000000000000000000000000 25
"))
cost.model <- data.table(read.table(header=TRUE,text="
    Linear        Log   Constant   min_mean   max_mean     data_i
       317       -185           366.415053    0.0000000000000000000000000000000000000000000000000000000    0.1237290654904118630819098711981496307998895645141601562 22
       358       -226           275.666058    0.1237290654904118630819098711981496307998895645141601562    0.6718768670266760389964133537432644516229629516601562500 21
       359       -228           274.198821    0.6718768670266760389964133537432644516229629516601562500    1.8567468357541714318870162969687953591346740722656250000 20
       372       -267           274.195324    1.8567468357541714318870162969687953591346740722656250000    2.8252700540160855524618455092422664165496826171875000000 19
        94       -185           974.454875    2.8252700540160855524618455092422664165496826171875000000   33.0000000000000000000000000000000000000000000000000000000 24
"))
one.env <- data.table(read.table(header=TRUE, text="
    Linear        Log   Constant   min_mean   max_mean     data_i
         0          0           901.495779    0.0000000000000000000000000000000000000000000000000000000    0.0835962145110410170811832131221308372914791107177734375 25
       317       -185           366.415053    0.0835962145110410170811832131221308372914791107177734375    0.1237290654904118630819098711981496307998895645141601562 22
       358       -226           275.666058    0.1237290654904118630819098711981496307998895645141601562    0.6718768670266760389964133537432644516229629516601562500 21
       359       -228           274.198821    0.6718768670266760389964133537432644516229629516601562500    1.8567468357541714318870162969687953591346740722656250000 20
       372       -267           274.195324    1.8567468357541714318870162969687953591346740722656250000    2.8160722315714470376235567528055980801582336425781250000 19
        91       -182           977.507459    2.8160722315714470376235567528055980801582336425781250000   33.0000000000000000000000000000000000000000000000000000000 25
"))
pieces <- data.table(read.table(header=TRUE, text="
    Linear        Log   Constant   min_mean   max_mean     data_i
         0          0           901.495779    0.0000000000000000000000000000000000000000000000000000000    1.6448087431693989568515235077938996255397796630859375000 25
       317       -185           366.415053    0.0000000000000000000000000000000000000000000000000000000    0.1237290654904118630819098711981496307998895645141601562 22
"))
curve(ploss(pieces[1,]-pieces[2,], x), 0, 0.1)
ploss <- function(dt, x){
  ## need to make a new data table, otherwise ifelse may only get one
  ## element, and return only one element.
  new.dt <- data.table(dt, x)
  new.dt[, ifelse(Log==0, 0, Log*log(x)) + Linear*x + Constant]
}
getLines <- function(dt, n.data=100){
  line.list <- list()
  for(piece.i in 1:nrow(dt)){
    piece <- dt[piece.i,]
    mean.vec <- piece[, seq(min_mean, max_mean, l=n.data)]
    line.list[[piece.i]] <- data.table(
      piece.i,
      piece,
      mean=mean.vec,
      cost=ploss(piece, mean.vec))
  }
  do.call(rbind, line.list)
}
ggplot()+
  geom_line(aes(mean, cost, color="min"), data=getLines(one.env), size=3)+
  geom_line(aes(mean, cost, color="prev"), data=getLines(compare.cost), size=1)+
  geom_line(aes(mean, cost, color="model"), data=getLines(cost.model), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))

ggplot()+
  geom_line(aes(mean, cost, color="min"), data=getLines(one.env), size=3)+
  geom_line(aes(mean, cost, color="prev"), data=getLines(compare.cost), size=1)+
  geom_line(aes(mean, cost, color="model"), data=getLines(cost.model), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))+
  coord_cartesian(xlim=c(431,443))

foo
## discriminant=-0x1.c3a787f67885p-446 two_roots=1 Linear=-1404 Log=2
## Constant=0x1.3b213cef7cf4p+9
## smaller_mean=-0.000000 -0x0p+0
## larger_mean=0.447758 0x1.ca81279a572ecp-2
## 0.447758 in [0.000000,0.663577]
## not equal on the sides, 1 crossing point
##     Linear        Log   Constant   min_mean   max_mean     data_i
##          0          0 -62254.401644   0.000000   0.447758 542
##       1404         -2 -62884.661316   0.447758   0.663577 540
## ------


## exactly equal on entire interval
ggplot()+
  ##geom_line(aes(mean, cost, color="min"), data=getLines(one.env), size=3)+
  geom_line(aes(mean, cost, color="prev"), data=getLines(compare.cost[3,]), size=1)+
  geom_line(aes(mean, cost, color="model"), data=getLines(cost.model[7,]), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))

## model 1023 -10778 smaller
ggplot()+
  ##geom_line(aes(mean, cost, color="min"), data=getLines(one.env), size=3)+
  geom_line(aes(mean, cost, color="prev"), data=getLines(compare.cost[3,]), size=1)+
  geom_line(aes(mean, cost, color="model"), data=getLines(cost.model[8,]), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))
