source("packages.R")

load("test.error.RData")

test.counts <- test.error[, list(
  errors=sum(errors),
  labels=sum(regions),
  tp=sum(tp),
  possible.tp=sum(possible.tp),
  fp=sum(fp),
  possible.fp=sum(possible.fp)
  ), by=.(set.name, set.i, algorithm, train.type)]
test.counts[, TPR := tp/possible.tp]
test.counts[, FPR := fp/possible.fp]

test.ranges <- test.counts[, list(
  min.labels=min(labels),
  max.labels=max(labels)
  ), by=.(set.name, set.i)]
test.ranges[, stopifnot(min.labels==max.labels)]

test.counts[, percent.error := errors/labels*100]

test.mean <- test.counts[, list(
  mean.percent=mean(percent.error)
), by=.(set.name, algorithm, train.type)]

roc.total <- roc[, list(
  tp=sum(tp),
  possible.tp=sum(possible.tp),
  fp=sum(fp),
  possible.fp=sum(possible.fp)
), by=.(set.name, set.i, algorithm, param.name)]
roc.total[, TPR := tp/possible.tp]
roc.total[, FPR := fp/possible.fp]
roc.ord <- roc.total[order(as.numeric(paste(param.name))),]

algo.colors <-
  c(hmcan="#A6CEE3",
    HMCanBroad="#1F78B4",
    triform="#B2DF8A",
    PeakSegDP="#33A02C",
    "#FB9A99",
    coseg="black",
    PeakSegJoint="grey40",
    Segmentor="#E31A1C",
    MACS="#FDBF6F",
    macs.broad="#FF7F00",
    ccat.tf="#CAB2D6", #lite purple
    ccat.histone="#6A3D9A",#dark purple
    "#FFFF99", #yellow
    homer="#B15928") #brown

## TODO: compute AUC?

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  facet_grid(. ~ set.name, labeller=function(df){
    df$set.name <- gsub("_", "\n", df$set.name)
    df
  })+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=c(unsupervised="white", supervised="grey"))+
  geom_path(aes(FPR, TPR, color=algorithm),
            data=roc.ord)+
  geom_point(aes(FPR, TPR, color=algorithm, fill=train.type),
             shape=21,
             data=test.counts)

levs <- c("MACS", "HMCanBroad", "PeakSegDP", "Segmentor", "coseg")
test.mean[, algo.fac := factor(algorithm, levs)]
test.counts[, algo.fac := factor(algorithm, levs)]
set.best <- test.mean[, list(min.percent=min(mean.percent)), by=set.name]
dots <- ggplot()+
  geom_vline(aes(xintercept=min.percent), data=set.best)+
  geom_point(aes(mean.percent, algo.fac, color=train.type),
             alpha=0.3,
             size=4,
             data=test.mean)+
  geom_point(aes(percent.error, algo.fac, color=train.type),
             data=test.counts, pch=1)+
  scale_color_discrete("penalty / threshold training method")+
  facet_grid(. ~ set.name, labeller=function(df){
    df$set.name <- gsub("_", "\n", df$set.name)
    df
  }, scales="free_y", space="free_y")+
  scale_y_discrete("algorithm")+
  theme_bw()+
  guides(color=guide_legend())+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_x_continuous("percent incorrect peak region labels (test error)",
                     breaks=seq(0, 60, by=20))
pdf("figure-test-error-dots.pdf", h=2.5, w=8)
print(dots)
dev.off()
