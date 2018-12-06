source("packages.R")

load("all.modelSelection.RData")
load("dp.peaks.sets.RData")

fold.dt <- data.table(set.name=names(dp.peaks.sets))[, {
  data.table(fold.i=1:4)[, {
    data.table(chunk.name=dp.peaks.sets[[set.name]][[fold.i]])
  }, by=list(fold.i)]
}, by=list(set.name)]

## Compute auc of naive linear penalty.
selection.folds <- all.modelSelection[fold.dt, on=list(set.name, chunk.name)]
selection.folds[, fp := total.fp]
selection.folds[, fn := total.fn]
selection.folds[, errors := total.errors]
linear.auc <- selection.folds[, {
  pred.dt <- .SD[segments==1, data.table(
    chunk.name, sample.id, pred.log.lambda=1)]
  L <- ROChange(.SD, pred.dt, c("chunk.name", "sample.id"))
  with(L, data.table(thresholds, auc))
}, by=list(set.name, fold.i, algo)]

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
several <- "305-321"
algo.map <- levs
names(algo.map) <- sub("[(].*", "", levs)
both.roc <- rbind(
  roc.total[, data.table(
    parameters="1", set.name, fold.i=set.i, algorithm, FPR, TPR)],
  roc[, data.table(
    parameters=several, set.name, fold.i, algorithm=algo.map[algo], FPR, TPR)])
roc.thresh[, algorithm := algo.map[algo] ]
roc.tall <- melt(
  roc.thresh[threshold=="predicted"],
  measure.vars=c("auc", "accuracy.percent"),
  id.vars=c("set.name", "fold.i", "algorithm"))
all.auc <- rbind(
  linear.auc[threshold=="predicted", data.table(
    parameters="1", penalty="linear",
    set.name, set.i=fold.i, algorithm=algo.map[algo], auc)],
  roc.thresh[threshold=="predicted", data.table(
    parameters=several, penalty="linear",
    set.name, set.i=fold.i, algorithm, auc)],
  data.table(parameters="1", penalty="oracle", auc))
all.auc[, algo.fac := factor(algorithm, levs)]

select.dt <- data.table(
  set.name="H3K36me3_TDH_immune",
  algorithm="PDPA(unconstrained)")
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ fold.i)+
  geom_path(aes(
    FPR, TPR, color=parameters),
    data=both.roc[select.dt, on=list(set.name, algorithm)])+
  coord_equal(xlim=c(0,1), ylim=c(0,1))
all.auc[select.dt, on=list(set.name, algorithm)]

mean.auc <- all.auc[, list(
  mean=mean(auc),
  sd=sd(auc),
  med=median(auc),
  q25=quantile(auc, 0.25),
  q75=quantile(auc, 0.75)
  ), by=.(set.name, algo.fac)]

if(FALSE){
  
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

  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(set.name ~ variable, scales="free")+
    geom_point(aes(
      value, algorithm),
      data=roc.tall)

}

dpa.only <- all.auc[grepl("DPA", algorithm)]
(fit <- lm(auc ~ set.name + algorithm + parameters, dpa.only))
anova(fit)

dpa.wide <- dcast(
  dpa.only,
  set.name + set.i + algorithm ~ penalty + parameters,
  value.var="auc")
ggplot()+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")+
  theme_bw()+
  geom_point(aes(
    linear_1, `linear_305-321`),
    data=dpa.wide)
(t.res <- dpa.wide[, t.test(linear_1, `linear_305-321`, paired=TRUE)])
## t = -3.4981, df = 83, p-value = 0.0007558
dpa.wide[, hist(linear_1-`linear_305-321`)]

ggplot()+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")+
  theme_bw()+
  geom_point(aes(
    oracle_1, linear_1),
    data=dpa.wide)
## t = -1.7734, df = 83, p-value = 0.07984
dpa.wide[, t.test(oracle_1, linear_1, paired=TRUE)]

dpa.tall <- melt(
  dpa.wide,
  id.vars=c("set.name", "set.i", "algorithm", "linear_1"))
dpa.pvals <- dpa.tall[, {
  with(t.test(value, linear_1, paired=TRUE), data.table(p.value, estimate, df=length(value)-1))
}, by=list(algorithm, variable)]
ggplot()+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(algorithm ~ variable, labeller=function(df){
    if("variable" %in% names(df)){
      df$variable <- sub("_", " penalty\nparameters=", df$variable)
    }
    df
  })+
  geom_point(aes(
    value, linear_1),
    shape=1,
    data=dpa.tall)+
  geom_text(aes(
    0.6, 0.9, label=sprintf("p=%.3f", p.value)),
    data=dpa.pvals)+
  xlab("Test AUC of competitor")+
  ylab("Test AUC of baseline linear penalty with one parameter")

gg <- ggplot()+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(variable ~ algorithm, labeller=function(df){
    if("variable" %in% names(df)){
      df$variable <- paste(
        "Competitor:",
        sub("_", " penalty\nparameters=", df$variable))
    }
    df
  })+
  geom_point(aes(
    linear_1, value),
    shape=1,
    data=dpa.tall)+
  geom_text(aes(
    0.67, 0.9, label=sprintf("diff=%.05f
p=%.3f", estimate, p.value)),
    data=dpa.pvals)+
  ylab("Test AUC of competitor")+
  xlab("Test AUC of baseline linear penalty with one parameter")
pdf("figure-all-cv-learned-oracle-compare.pdf", 10, 5)
print(gg)
dev.off()

show.auc <- all.auc[
  ifelse(
    grepl("DPA", algorithm),
    parameters==several & penalty=="linear",
    TRUE)]
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
    shape=21,
    data=show.auc)+
  scale_y_discrete(
    "Algorithm")+
  scale_x_continuous(
    "Test AUC (larger values indicate more accurate peak detection)",
    breaks=seq(0.6, 1, by=0.2))
print(gg)

show.wide <- dcast(
  show.auc,
  set.name + set.i ~ algorithm,
  value.var="auc")
show.mean <- show.auc[, list(
  mean=mean(auc)
), by=list(set.name, algo=sub("[(].*", "", algorithm))]
show.mean.wide <- dcast(
  show.mean[algo %in% c("CDPA", "GPDPA")],
  set.name ~ algo,
  value.var="mean")
(show.pvals <- show.wide[, {
  with(t.test(
    `CDPA(previous best)`,
    `GPDPA(proposed)`,
    paired=TRUE),
    data.table(p.value, estimate))
}, by=list(set.name)][show.mean.wide, on=list(set.name)])

show.mean.wide <- dcast(
  show.mean[algo %in% c("PDPA", "GPDPA")],
  set.name ~ algo,
  value.var="mean")
(show.pvals <- show.wide[, {
  with(t.test(
    `PDPA(unconstrained)`,
    `GPDPA(proposed)`,
    paired=TRUE),
    data.table(p.value, estimate))
}, by=list(set.name)][show.mean.wide, on=list(set.name)])

pdf("figure-all-cv.pdf", 8, 2.2)
print(gg)
dev.off()
