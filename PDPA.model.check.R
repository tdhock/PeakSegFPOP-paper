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
         2        -12      -1253155.345559    0.0000000000000000000000000000000000000000000000000000000    0.0000000033308900973693182148977841685374284841358871745 27666
      3246      -1461      -1281439.863638    0.0000000033308900973693182148977841685374284841358871745    0.0000007708805774640572242392038207481341771654115291312 27666
      4213      -2056      -1289814.925143    0.0000007708805774640572242392038207481341771654115291312    0.0009941829803357931510943146236058964859694242477416992 27666
     21666     -23579      -1438633.458765    0.0009941829803357931510943146236058964859694242477416992    0.0069460097764935434597188823602209595264866948127746582 27666
     34993     -50321      -1571622.748318    0.0069460097764935434597188823602209595264866948127746582    0.2181121577155594737362775958899874240159988403320312500 27666
     35065     -50609      -1572077.003202    0.2181121577155594737362775958899874240159988403320312500    1.1120034748951577174835847472422756254673004150390625000 27666
     35278     -51895      -1572177.333912    1.1120034748951577174835847472422756254673004150390625000    1.4710301037473778595909834621124900877475738525390625000 27666
         0          0      -1540311.878931    1.4710301037473778595909834621124900877475738525390625000   58.0000000000000000000000000000000000000000000000000000000 27666
"))
cost.model <- data.table(read.table(header=TRUE,text="
    Linear        Log   Constant   min_mean   max_mean     data_i
       947       -679      -1266581.561075    0.0000000000000000000000000000000000000000000000000000000    0.0000000056006291739902859275065398936450933309316724262 27652
      2709      -1461      -1281439.863638    0.0000000056006291739902859275065398936450933309316724262    0.0000007708807582382837898237761507791443449377766228281 27636
      3857      -2056      -1289814.925143    0.0000007708807582382837898237761507791443449377766228281    0.0009941810500210428776929338923196155519690364599227905 27624
     21268     -23579      -1438633.458765    0.0009941810500210428776929338923196155519690364599227905    0.0069467304417468135266378048697788472054526209831237793 27304
     34993     -50321      -1571622.748318    0.0069467304417468135266378048697788472054526209831237793    0.1576360253831302626892352236609440296888351440429687500 27665
     34790     -50289      -1571531.629276    0.1576360253831302626892352236609440296888351440429687500    0.2197009389262863654490587350665009580552577972412109375 26939
     35065     -50609      -1572077.003202    0.2197009389262863654490587350665009580552577972412109375    1.1120034748951577174835847472422756254673004150390625000 27665
     35278     -51895      -1572177.333912    1.1120034748951577174835847472422756254673004150390625000    1.4679698138901497994623923659673891961574554443359375000 27665
        30       -152      -1540297.455754    1.4679698138901497994623923659673891961574554443359375000   11.9480290681335201696811054716818034648895263671875000000 27664
         2        -12      -1540310.190227   11.9480290681335201696811054716818034648895263671875000000   58.0000000000000000000000000000000000000000000000000000000 27665
"))
one.env <- data.table(read.table(header=TRUE, text="
    Linear        Log   Constant   min_mean   max_mean     data_i
         2        -12      -1253155.345559    0.0000000000000000000000000000000000000000000000000000000    0.0000000033308900973693182148977841685374284841358871745 27666
       947       -679      -1266581.561075    0.0000000033308900973693182148977841685374284841358871745    0.0000000056006291739902859275065398936450933309316724262 27652
      2709      -1461      -1281439.863638    0.0000000056006291739902859275065398936450933309316724262    0.0000007708807582382837898237761507791443449377766228281 27636
      3857      -2056      -1289814.925143    0.0000007708807582382837898237761507791443449377766228281    0.0009941810500210428776929338923196155519690364599227905 27624
     21268     -23579      -1438633.458765    0.0009941810500210428776929338923196155519690364599227905    0.0069467304417468135266378048697788472054526209831237793 27304
     34993     -50321      -1571622.748318    0.0069467304417468135266378048697788472054526209831237793    0.1576360253831302626892352236609440296888351440429687500 27666
     34790     -50289      -1571531.629276    0.1576360253831302626892352236609440296888351440429687500    0.2197009389262863654490587350665009580552577972412109375 26939
     35065     -50609      -1572077.003202    0.2197009389262863654490587350665009580552577972412109375    1.1120034748951577174835847472422756254673004150390625000 27666
     35278     -51895      -1572177.333912    1.1120034748951577174835847472422756254673004150390625000    1.4679698138901497994623923659673891961574554443359375000 27666
        30       -152      -1540297.455754    1.4679698138901497994623923659673891961574554443359375000   11.9480290681335201696811054716818034648895263671875000000 27664
         2        -12      -1540310.190227   11.9480290681335201696811054716818034648895263671875000000   15.6637381935888377881838096072897315025329589843750000000 27665
         0          0      -1540311.878931   15.6637381935888377881838096072897315025329589843750000000   58.0000000000000000000000000000000000000000000000000000000 27666
"))
pieces <- data.table(read.table(header=TRUE, text="
    Linear        Log   Constant   min_mean   max_mean     data_i
         0          0           901.495779    0.0000000000000000000000000000000000000000000000000000000    1.6448087431693989568515235077938996255397796630859375000 25
       317       -185           366.415053    0.0000000000000000000000000000000000000000000000000000000    0.1237290654904118630819098711981496307998895645141601562 22
"))
ploss <- function(dt, x){
  ## need to make a new data table, otherwise ifelse may only get one
  ## element, and return only one element.
  new.dt <- data.table(dt, x)
  new.dt[, ifelse(Log==0, 0, Log*log(x)) + Linear*x + Constant]
}
curve(ploss(pieces[1,]-pieces[2,], x), 0, 0.1)
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
  coord_cartesian(xlim=c(0, 1e-8))

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
