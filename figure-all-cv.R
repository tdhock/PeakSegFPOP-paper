source("packages.R")

(tobjs <- load("test.error.RData"))
levs <- c(
  MACS="MACS(baseline)",
  HMCanBroad="HMCanBroad(baseline)",
  Segmentor="PDPA(unconstrained)",
  PeakSegDP="CDPA(previous best)",
  coseg.inf="GPDPA(proposed)")
test.error[, algorithm := levs[algorithm] ]
roc[, algorithm := levs[algorithm] ]
test.counts <- test.error[algorithm %in% levs, list(
  errors=sum(errors),
  labels=sum(regions),
  tp=sum(tp),
  possible.tp=sum(possible.tp),
  fp=sum(fp),
  possible.fp=sum(possible.fp)
  ), by=.(set.name, set.i, algorithm, train.type)]
test.counts[, TPR := tp/possible.tp]
test.counts[, FPR := fp/possible.fp]
possible.counts <- test.counts[algorithm==algorithm[1] & train.type=="supervised", {
  list(
    possible.fn=sum(possible.tp),
    possible.fp=sum(possible.fp)
    )
}, by=list(set.name)]
test.ranges <- test.counts[, list(
  min.labels=min(labels),
  max.labels=max(labels)
  ), by=.(set.name, set.i)]
test.ranges[, stopifnot(min.labels==max.labels)]
test.counts[, percent.accuracy := (1-errors/labels)*100]
test.mean <- test.counts[, list(
  mean.percent=mean(percent.accuracy)
), by=.(set.name, algorithm, train.type)]
roc.total <- roc[algorithm %in% levs, list(
  tp=sum(tp),
  possible.tp=sum(possible.tp),
  fp=sum(fp),
  possible.fp=sum(possible.fp)
), by=.(set.name, set.i, algorithm, param.i)]
roc.total[, TPR := tp/possible.tp]
roc.total[, FPR := fp/possible.fp]
roc.ord <- roc.total[order(param.i),]
roc.not.cvx <- roc.total[, {
  tpr <- TPR[which.max(FPR)]
  tpr <- 1
  list(
    FPR=c(1, FPR, 0, 1, 1),
    TPR=c(tpr, TPR, 0, 0, tpr)
    )
}, by=.(set.name, set.i, algorithm)]
auc <- roc.not.cvx[, list(
  auc=geometry::polyarea(FPR, TPR)
), by=.(set.name, set.i, algorithm)]

(objs <- load("all.cv.RData"))
algo.map <- levs
names(algo.map) <- sub("[(].*", "", levs)
roc.thresh[, algorithm := algo.map[algo] ]
roc.tall <- melt(
  roc.thresh[threshold=="predicted"],
  measure.vars=c("auc", "accuracy.percent"),
  id.vars=c("set.name", "fold.i", "algorithm"))
all.auc <- rbind(
  roc.thresh[threshold=="predicted", data.table(
    set.name, set.i=fold.i, algorithm, auc)],
  auc[grepl("baseline", algorithm)])
all.auc[, algo.fac := factor(algorithm, levs)]

mean.auc <- all.auc[, list(
  mean=mean(auc),
  sd=sd(auc),
  med=median(auc),
  q25=quantile(auc, 0.25),
  q75=quantile(auc, 0.75)
  ), by=.(set.name, algo.fac)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ set.name, scales="free")+
  geom_point(aes(
    med, algo.fac),
    shape=1,
    data=mean.auc)+
  geom_segment(aes(
    q25, algo.fac,
    xend=q75, yend=algo.fac),
    data=mean.auc)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ set.name, labeller=function(df){
    df$set.name <- gsub("_", "\n", df$set.name)
    df$set.name <- paste0(
      ifelse(
        grepl("H3K36me3", df$set.name),
        "Broad", "Sharp"),
      "\n",
      df$set.name)
    df
  })+
  geom_point(aes(
    auc, algo.fac),
    shape=1,
    data=all.auc)+
  scale_y_discrete(
    "Algorithm")+
  scale_x_continuous(
    "Test AUC (larger values indicate more accurate peak detection)",
    breaks=seq(0.6, 1, by=0.2))


ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(set.name ~ variable, scales="free")+
  geom_point(aes(
    value, algorithm),
    data=roc.tall)

pdf("figure-all-cv.pdf", 8, 2.2)
print(gg)
dev.off()
