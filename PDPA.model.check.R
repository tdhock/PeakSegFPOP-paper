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
         0          0           841.885175    0.0000000000000000000000000000000000000000000000000000000    1.5943837753510139432222558752982877194881439208984375000 24
       641      -1022           296.635210    1.5943837753510139432222558752982877194881439208984375000    1.6915632528339297202535362885100767016410827636718750000 24
       570       -951           379.414830    1.6915632528339297202535362885100767016410827636718750000    2.7263919513179435405447748053120449185371398925781250000 24
         3         -3           974.454875    2.7263919513179435405447748053120449185371398925781250000   33.0000000000000000000000000000000000000000000000000000000 24
"))
cost.model <- data.table(read.table(header=TRUE,text="
    Linear        Log   Constant   min_mean   max_mean     data_i
       226         -3           366.415053    0.0000000000000000000000000000000000000000000000000000000    0.1237290654904118630819098711981496307998895645141601562 22
       267        -44           275.666058    0.1237290654904118630819098711981496307998895645141601562    0.6718768670266760389964133537432644516229629516601562500 21
       268        -46           274.198821    0.6718768670266760389964133537432644516229629516601562500    1.8567468357541716539316212220001034438610076904296875000 20
       281        -85           274.195324    1.8567468357541716539316212220001034438610076904296875000    3.3402697894187407534616340853972360491752624511718750000 19
       321       -245           333.552785    3.3402697894187407534616340853972360491752624511718750000    3.4954719311684741001045040320605039596557617187500000000 18
       341       -329           368.766692    3.4954719311684741001045040320605039596557617187500000000    4.8760122516499011524615525559056550264358520507812500000 16
       362       -434           432.724846    4.8760122516499011524615525559056550264358520507812500000    5.1521700820567168932484491961076855659484863281250000000 15
       408       -671           584.267088    5.1521700820567168932484491961076855659484863281250000000    8.2042789713824166852873531752265989780426025390625000000 23
       226         -3           671.535757    8.2042789713824166852873531752265989780426025390625000000   33.0000000000000000000000000000000000000000000000000000000 23
"))
one.env <- data.table(read.table(header=TRUE, text="
    Linear        Log   Constant   min_mean   max_mean     data_i
       226         -3           366.415053    0.0000000000000000000000000000000000000000000000000000000    0.1237290654904118630819098711981496307998895645141601562 22
       267        -44           275.666058    0.1237290654904118630819098711981496307998895645141601562    0.6718768670266760389964133537432644516229629516601562500 21
       268        -46           274.198821    0.6718768670266760389964133537432644516229629516601562500    1.8567468357541716539316212220001034438610076904296875000 20
       281        -85           274.195324    1.8567468357541716539316212220001034438610076904296875000    2.8252700540160855524618455092422664165496826171875000000 19
         3         -3           974.454875    2.8252700540160855524618455092422664165496826171875000000   33.0000000000000000000000000000000000000000000000000000000 24
"))
compare.cost.log <- data.table(read.table(header=TRUE,text="
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0           841.885175                                                         -inf    0.4664873138429796450843412003450794145464897155761718750 24
       641      -1022           296.635210    0.4664873138429796450843412003450794145464897155761718750    0.5256531030614487454144523326249327510595321655273437500 24
       570       -951           379.414830    0.5256531030614487454144523326249327510595321655273437500    1.0029791055204002603318258479703217744827270507812500000 24
         3         -3           974.454875    1.0029791055204002603318258479703217744827270507812500000    3.4965075614664802294839773821877315640449523925781250000 24
"))
cost.model.log <- data.table(read.table(header=TRUE,text="
    Linear        Log   Constant min_log_mean max_log_mean     data_i
       226         -3           366.415053                                                         -inf   -2.0896610595980571467578101874096319079399108886718750000 22
       267        -44           275.666058   -2.0896610595980571467578101874096319079399108886718750000   -0.3976801888396973017059110588888870552182197570800781250 21
       268        -46           274.198821   -0.3976801888396973017059110588888870552182197570800781250    0.6188259433806467812999585476063657552003860473632812500 20
       281        -85           274.195324    0.6188259433806467812999585476063657552003860473632812500    1.2060515790015347015184943302301689982414245605468750000 19
       321       -245           333.552785    1.2060515790015347015184943302301689982414245605468750000    1.2514683969472388813670704621472395956516265869140625000 18
       341       -329           368.766692    1.2514683969472388813670704621472395956516265869140625000    1.5843277242594127063313180769910104572772979736328125000 16
       362       -434           432.724846    1.5843277242594127063313180769910104572772979736328125000    1.6394187446460362078681782804778777062892913818359375000 15
       408       -671           584.267088    1.6394187446460362078681782804778777062892913818359375000    2.1046558439448257438186828949255868792533874511718750000 23
       226         -3           671.535757    2.1046558439448257438186828949255868792533874511718750000    3.4965075614664802294839773821877315640449523925781250000 23
"))
one.env.log <- data.table(read.table(header=TRUE, text="
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0           841.885175                                                         -inf -158.4900405221358994367619743570685386657714843750000000000 24
       226         -3           366.415053 -158.4900405221358994367619743570685386657714843750000000000   -2.0896610595980571467578101874096319079399108886718750000 22
       267        -44           275.666058   -2.0896610595980571467578101874096319079399108886718750000   -0.3976801888396973017059110588888870552182197570800781250 21
       268        -46           274.198821   -0.3976801888396973017059110588888870552182197570800781250    0.6188259433806467812999585476063657552003860473632812500 20
       281        -85           274.195324    0.6188259433806467812999585476063657552003860473632812500    1.0386039543715330779605210409499704837799072265625000000 19
         3         -3           974.454875    1.0386039543715330779605210409499704837799072265625000000    3.4965075614664802294839773821877315640449523925781250000 24
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
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ x, scales="free")+
  geom_line(aes(log(mean), cost, color="min"), data=LOG(one.env.log), size=3)+
  geom_line(aes(log(mean), cost, color="prev"), data=LOG(compare.cost.log), size=1)+
  geom_line(aes(log(mean), cost, color="model"), data=LOG(cost.model.log), size=1)+
  geom_line(aes(mean, cost, color="min"), data=ID(one.env), size=3)+
  geom_line(aes(mean, cost, color="prev"), data=ID(compare.cost), size=1)+
  geom_line(aes(mean, cost, color="model"), data=ID(cost.model), size=1)+
  scale_color_manual(values=c(min="grey", prev="red", model="black"))

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
