source("packages.R")

PDPA.RData.vec <- Sys.glob("data/H*/*/PDPA.model.RData")

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
"))
cost.model <- data.table(read.table(header=TRUE,text="
"))
one.env <- data.table(read.table(header=TRUE, text="
"))
compare.cost.log <- data.table(read.table(header=TRUE,text="
    Linear        Log   Constant min_log_mean max_log_mean     data_i
       156          0        -33637.990973                                                         -inf  -44.0229655582002976643707370385527610778808593750000000000 488
      1437       -198        -42354.538153  -44.0229655582002976643707370385527610778808593750000000000   -8.6199210595184911909427682985551655292510986328125000000 484
      4190      -3166        -67938.960704   -8.6199210595184911909427682985551655292510986328125000000   -3.2365403010716939924407142825657501816749572753906250000 440
      4193      -3175        -67968.207466   -3.2365403010716939924407142825657501816749572753906250000   -1.6763779035441404019479705311823636293411254882812500000 439
      4211      -3247        -68092.273579   -1.6763779035441404019479705311823636293411254882812500000   -0.7298660327086520238637490365363191813230514526367187500 438
      4254      -3462        -68269.919639   -0.7298660327086520238637490365363191813230514526367187500   -0.0625264321989785765154579166846815496683120727539062500 437
      4259      -3492        -68276.492373   -0.0625264321989785765154579166846815496683120727539062500    0.3480945451874953855408989511488471180200576782226562500 436
      4288      -3689        -68248.992367    0.3480945451874953855408989511488471180200576782226562500    1.4713604348927762366372462565777823328971862792968750000 432
      4289      -3699        -68238.633918    1.4713604348927762366372462565777823328971862792968750000    1.7150324780329941898315837534028105437755584716796875000 431
      4292      -3732        -68198.708414    1.7150324780329941898315837534028105437755584716796875000    2.3972943243266056079221471009077504277229309082031250000 430
      5323     -19670        -41324.818164    2.3972943243266056079221471009077504277229309082031250000    2.4550309613452787615983652358409017324447631835937500000 245
      5325     -19704        -41264.640700    2.4550309613452787615983652358409017324447631835937500000    2.5497462183816694825111426325747743248939514160156250000 244
       201        -45        -25783.153775    2.5497462183816694825111426325747743248939514160156250000    4.0943445622221004143170830502640455961227416992187500000 488
"))
cost.model.log <- data.table(read.table(header=TRUE,text="
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0        -33637.990973                                                         -inf    4.0943445622221004143170830502640455961227416992187500000 489
"))
one.env.log <- data.table(read.table(header=TRUE, text="
"))
ploss <- function(dt, x){
  ## need to make a new data table, otherwise ifelse may only get one
  ## element, and return only one element.
  new.dt <- data.table(dt, x)
  new.dt[, ifelse(Log==0, 0, Log*log(x)) + Linear*x + Constant]
}
LOG <- function(dt, n.data=1000){
  line.list <- list()
  for(piece.i in 1:nrow(dt)){
    piece <- dt[piece.i,]
    mean.vec <- piece[, seq(exp(min_log_mean), exp(max_log_mean), l=n.data)]
    line.list[[piece.i]] <- data.table(
      piece.i,
      piece,
      mean=mean.vec,
      cost=ploss(piece, mean.vec))
  }
  data.frame(do.call(rbind, line.list), x="log(mean)")
}
ID <- function(dt, n.data=1000){
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
  data.frame(do.call(rbind, line.list), x="mean")
}
ggplot()+
  geom_line(aes(log(mean), cost, color="min"), data=LOG(one.env.log), size=3)+
  geom_line(aes(log(mean), cost, color="prev"), data=LOG(compare.cost.log), size=1)+
  geom_line(aes(log(mean), cost, color="model"), data=LOG(cost.model.log), size=1)+
  ## geom_line(aes(mean, cost, color="min"), data=ID(one.env), size=3)+
  ## geom_line(aes(mean, cost, color="prev"), data=ID(compare.cost), size=1)+
  ## geom_line(aes(mean, cost, color="model"), data=ID(cost.model), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ x, scales="free")

ggplot()+
  geom_line(aes(mean, cost, color="prev"), data=LOG(compare.cost.log), size=1)+
  geom_line(aes(mean, cost, color="model"), data=LOG(cost.model.log), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ x, scales="free")

ggplot()+
  geom_line(aes(log(mean), cost, color="prev"), data=LOG(compare.cost.log), size=1)+
  geom_line(aes(log(mean), cost, color="model"), data=LOG(cost.model.log), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))+
  coord_cartesian(xlim=c(-44.1, -44), ylim=c(-33640,-33637.5), expand=FALSE)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ x, scales="free")

ggplot()+
  geom_line(aes(log(mean), cost, color="min"), data=LOG(one.env.log), size=3)+
  geom_line(aes(log(mean), cost, color="prev"), data=LOG(compare.cost.log), size=1)+
  geom_line(aes(log(mean), cost, color="model"), data=LOG(cost.model.log), size=1)+
  ## geom_line(aes(mean, cost, color="min"), data=ID(one.env), size=3)+
  ## geom_line(aes(mean, cost, color="prev"), data=ID(compare.cost), size=1)+
  ## geom_line(aes(mean, cost, color="model"), data=ID(cost.model), size=1)+
  geom_vline(xintercept=2.8481847238557489454535698314430192112922668457031250000)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  coord_cartesian(xlim=c(2.8475, 2.85), ylim=c(-2860000,-2850000), expand=FALSE)+
  facet_grid(. ~ x, scales="free")

ggplot()+
  geom_line(aes(mean, cost, color="min"), data=ID(one.env), size=3)+
  geom_line(aes(mean, cost, color="prev"), data=ID(compare.cost), size=1)+
  geom_line(aes(mean, cost, color="model"), data=ID(cost.model), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))+
  coord_cartesian(xlim=c(0,0.1), ylim=c(0, 2000))

foo

pieces <- data.table(read.table(header=TRUE, text="
    Linear        Log   Constant min_log_mean max_log_mean     data_i
       256          0            76.991564                                                         -inf    0.5523746266788631675836995782447047531604766845703125000 1
"))
curve(ploss(pieces[1,], x), -1, 0.5)

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
