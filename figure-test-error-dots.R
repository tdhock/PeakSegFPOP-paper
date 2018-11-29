source("packages.R")
library(animint)

load("test.error.RData")

levs <- c(
  MACS="MACS",
  HMCanBroad="HMCanBroad",
  Segmentor="PDPA",
  PeakSegDP="CDPA",
  ##coseg.inf="GPDPAinf",
  coseg="GPDPA")
levs <- c(
  MACS="MACS(popular baseline)",
  HMCanBroad="HMCanBroad(popular baseline)",
  Segmentor="PDPA(unconstrained baseline)",
  PeakSegDP="CDPA(previous best)",
  ##coseg.inf="GPDPAinf",
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

algo.colors <- c(
  hmcan="#A6CEE3",
  HMCanBroad="#1F78B4",
  triform="#B2DF8A",
  PeakSegDP="#33A02C",
  CDPA="#33A02C",
  "#FB9A99",
  coseg="black",
  GPDPA="black",
  GPDPAinf="grey",
  PeakSegJoint="grey40",
  "#E31A1C",
  MACS="#FDBF6F",
  macs.broad="#FF7F00",
  ccat.tf="#CAB2D6", #lite purple
  "#6A3D9A",#dark purple
  "#FFFF99", #yellow
  PDPA="#B15928", #brown
  Segmentor="#B15928") #brown

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  facet_grid(set.i ~ set.name, labeller=function(df){
    if("set.name" %in% names(df)){
      df$set.name <- gsub("_", "\n", df$set.name)
    }
    df
  })+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=c(unsupervised="white", supervised="black"))+
  geom_path(aes(FPR, TPR, color=algorithm, group=paste(set.i, algorithm)),
            data=roc.total)+
  geom_point(aes(FPR, TPR, color=algorithm, fill=train.type),
             shape=21,
             stroke=1,
             data=test.counts)

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
auc.wide <- dcast(auc, set.name + set.i ~ algorithm, value.var="auc")

ggplot()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(`CDPA(previous best)`, `GPDPA(proposed)`), data=auc.wide)+
  coord_equal()

err.wide <- dcast(test.counts[train.type=="supervised"], set.name+set.i ~ algorithm, value.var="percent.accuracy")
ggplot()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(`CDPA(previous best)`, `GPDPA(proposed)`), data=err.wide)+
  coord_equal()

ggplot()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(`PDPA(unconstrained baseline)`, `GPDPA(proposed)`), data=err.wide)+
  coord_equal()

ggplot()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(`PDPA(unconstrained baseline)`, `CDPA(previous best)`), data=err.wide)+
  coord_equal()

if(FALSE){
  auc.wide[, {
    L <- t.test(`CDPA(previous best)`, `GPDPA(proposed)`, paired=TRUE)
    res <- L[c("statistic", "p.value")]
    res$mean.cdpa <- mean(`CDPA(previous best)`)
    res$mean.gpdpa <- mean(`GPDPA(proposed)`)
    res
  }, by=list(set.name)][order(p.value)]
  auc.wide[, {
    L <- t.test(`PDPA(unconstrained baseline)`, `GPDPA(proposed)`, paired=TRUE)
    res <- L[c("statistic", "p.value")]
    res$mean.pdpa <- mean(`PDPA(unconstrained baseline)`)
    res$mean.gpdpa <- mean(`GPDPA(proposed)`)
    res
  }, by=list(set.name)][order(p.value)]
}

roc.cvx <- roc.not.cvx[, {
  fit <- chull(FPR, TPR)
  data.table(FPR, TPR)[fit,]
}, by=.(set.name, set.i, algorithm)]
mean.auc <- auc[, list(
  mean.auc=mean(auc),
  sd.auc=sd(auc),
  med.auc=median(auc),
  q25.auc=quantile(auc, 0.25),
  q75.auc=quantile(auc, 0.75)
  ), by=.(set.name, algorithm)]
## coseg is the best over all data sets, in terms of AUC.
mean.auc[, list(mean=mean(mean.auc)), by=algorithm][order(mean),]
best.auc <- mean.auc[, list(min.auc=max(mean.auc)), by=set.name]

some.algos <- function(dt)dt
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  facet_grid(set.i ~ set.name, labeller=function(df){
    if("set.name" %in% names(df)){
      df$set.name <- gsub("_", "\n", df$set.name)
      df$set.name <- paste0(
        ifelse(
          grepl("H3K36me3", df$set.name),
          "Broad", "Sharp"),
        "\n",
        df$set.name)
    }
    df
  })+
  geom_polygon(aes(FPR, TPR, color=algorithm),
               data=some.algos(roc.not.cvx),
               fill="grey")+
  ##scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=c(unsupervised="white", supervised="black"))+
  geom_path(aes(FPR, TPR, color=algorithm, group=paste(set.i, algorithm)),
            data=some.algos(roc.total))+
  coord_cartesian(xlim=c(0,0.5), ylim=c(0.5,1))+
  geom_point(aes(FPR, TPR, color=algorithm, fill=train.type),
             shape=21,
             stroke=1,
             data=some.algos(test.counts))

roc.show <- data.table(roc.not.cvx)
roc.show[FPR==1 & TPR==0, FPR := NA]
only <- function(dt)dt[set.name=="H3K4me3_TDH_immune" & set.i==2]
linetype.vec <- c(
  "error rate of predicted peaks\nat various pruning thresholds"="solid",
  "extrapolation to compute AUC"="dotted")
N.vec <- names(linetype.vec)
cfac <- function(i){
  factor(N.vec[i], N.vec)
}
both.roc <- rbind(
  data.table(curve=cfac(2), roc.show),
  data.table(curve=cfac(1), roc.total[, names(roc.show), with=FALSE]))
leg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position=c(0.8, 0.3))+
  coord_equal(xlim=c(0,1.2))+
  scale_x_continuous("FPR = False Positive Rate of predicted peaks relative to gold standard labels from biologist", breaks=seq(0, 1, by=0.2))+
  scale_y_continuous("TPR = True Positive Rate of predicted peaks relative to gold standard labels from biologist", breaks=seq(0, 1, by=0.2))+
  geom_path(aes(
    FPR, TPR,
    color=algorithm,
    linetype=curve,
    group=paste(curve, algorithm)),
            size=0.7,
            data=only(both.roc))+
  scale_linetype_manual(values=linetype.vec, guide=guide_legend(keyheight=2))+
  ## geom_path(aes(FPR, TPR, color=algorithm, group=paste(set.i, algorithm)),
  ##           data=only(roc.show))+
  ## geom_path(aes(FPR, TPR, color=algorithm, group=paste(set.i, algorithm)),
  ##           size=3,
  ##           data=only(roc.total))+
  geom_text(aes(FPR, TPR, label=label), data={
    rbind(
      data.table(FPR=1, TPR=0.6, label="Extrapolation to FPR=TPR=1\nto compute AUC"),
      data.table(FPR=0.3, TPR=0.6, label="No peak pruning\nmost predicted peaks"),
      data.table(FPR=0.4, TPR=0.4, label="All peaks pruned\nno peaks predicted\nFPR=TPR=0"))
  })+
  geom_segment(aes(x, y, xend=xend, yend=yend), arrow=grid::arrow(type="closed"), data={
    rbind(
      data.table(x=1, xend=1, y=0.7, yend=0.9),
      data.table(x=0.3, xend=0.3, y=0.7, yend=0.8),
      data.table(x=0.3, xend=0.1, y=0.3, yend=0.1))
  })
pdf("figure-test-error-dots-ROC-supp.pdf",8,8)
print(leg)
dev.off()
png("figure-test-error-dots-ROC-supp.png",800,800,res=100)
print(leg)
dev.off()
##direct.label(leg, "smart.grid")


test.mean[, algo.fac := factor(algorithm, levs)]
auc[, algo.fac := factor(algorithm, levs)]
mean.auc[, algo.fac := factor(algorithm, levs)]
test.counts[, algo.fac := factor(algorithm, levs)]
set.best <- test.mean[, list(min.percent=max(mean.percent)), by=set.name]
setkey(possible.counts, set.name)
dots <- ggplot()+
  geom_vline(aes(xintercept=min.percent), data=set.best)+
  geom_point(aes(mean.percent, algo.fac, color=train.type),
             alpha=0.3,
             size=4,
             data=test.mean)+
  geom_point(aes(percent.accuracy, algo.fac, color=train.type),
             data=test.counts, pch=1)+
  scale_color_discrete("penalty / threshold training method")+
  facet_grid(. ~ set.name, labeller=function(df){
    count.dt <- possible.counts[df$set.name]
    df$set.name <- paste0(
      gsub("_", "\n", df$set.name),
      "\n", count.dt$possible.fp, " possible fp",
      "\n", count.dt$possible.fn, " possible fn")
    df
  }, scales="free_y", space="free_y")+
  scale_y_discrete("algorithm")+
  theme_bw()+
  guides(color=guide_legend())+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_x_continuous("percent correct peak region labels (test accuracy)",
                     breaks=seq(40, 100, by=20))
dots
pdf("figure-test-error-dots-supp.pdf", 10, 3)
print(dots)
dev.off()

test.counts[, testSet := paste(set.name, "split", set.i)]
auc[, testSet := paste(set.name, "split", set.i)]
roc.total[, testSet := paste(set.name, "split", set.i)]
roc.not.cvx[, testSet := paste(set.name, "split", set.i)]
viz <- list(
  title="Peak detection test accuracy (4-fold cross-validation)",
  error=ggplot()+
    ggtitle("Click to select test set")+
    scale_fill_manual(values=c(unsupervised="white", supervised="black"))+
    geom_vline(aes(xintercept=min.percent),
               color="grey",
               data=set.best)+
    geom_point(aes(mean.percent, algo.fac, color=algorithm,
                  showSelectedalgo=algorithm,
                   showSelected=train.type),
               alpha=0.3,
               size=4,
               data=test.mean)+
    geom_point(aes(percent.accuracy, algo.fac, color=algorithm,
                  showSelectedalgo=algorithm,
                   clickSelects=testSet,
                   fill=train.type),
               stroke=1,
               size=3,
               data=test.counts,
               shape=21)+
    scale_color_manual(values=algo.colors)+
    facet_grid(. ~ set.name, labeller=function(df){
      df$set.name <- gsub("_", "\n", df$set.name)
      df
    })+
    guides(color="none")+
    scale_y_discrete("algorithm")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"),
          legend.position="top")+
    theme_animint(width=1100, height=170)+
    scale_x_continuous(
      "Test accuracy (percent correct peak region labels)",
      breaks=c(40, 60, 80, 100)),
  auc=ggplot()+
    geom_vline(aes(xintercept=min.auc),
               color="grey",
               data=best.auc)+
    scale_color_manual(values=algo.colors)+
    guides(color="none")+
    geom_point(aes(mean.auc, algo.fac,
                  showSelectedalgo=algorithm,
                   color=algorithm),
               alpha=0.3,
               size=4,
               data=mean.auc)+
    geom_point(aes(auc, algo.fac, color=algorithm,
                  showSelectedalgo=algorithm,
                   clickSelects=testSet),
               size=3,
               data=auc, pch=1)+
    facet_grid(. ~ set.name, labeller=function(df){
      df$set.name <- gsub("_", "\n", df$set.name)
      df
    })+
    scale_y_discrete("algorithm")+
    theme_bw()+
    theme_animint(width=1100, height=140)+
    theme(panel.margin=grid::unit(0, "cm"),
          legend.position="top")+
    scale_x_continuous(
      "Test AUC (larger values indicate more accurate peak detection)",
      breaks=c(0.6, 0.8, 1),
      labels=c("0.6", "0.8", "1")),
  rocZoom=ggplot()+
    ggtitle("Upper left corner of ROC space")+
    guides(color="none", fill="none")+
    scale_color_manual(values=algo.colors)+
    scale_fill_manual(values=c(unsupervised="white", supervised="black"))+
    geom_path(aes(FPR, TPR, color=algorithm,
                  showSelected=testSet,
                  showSelectedalgo=algorithm,
                  group=algorithm),
              alpha=0.5,
              data=roc.total)+
    coord_cartesian(xlim=c(0,0.5), ylim=c(0.5,1))+
    geom_point(aes(FPR, TPR, color=algorithm,
                  showSelectedalgo=algorithm,
                   showSelected=testSet,
                   showSelectedtype=train.type,
                   fill=train.type),
               shape=21,
               stroke=1,
               data=test.counts),
  roc=ggplot()+
    ggtitle("AUC computation (full ROC space)")+
    scale_color_manual(values=algo.colors)+
    scale_fill_manual(values=c(unsupervised="white", supervised="black"))+
    geom_polygon(aes(FPR, TPR, color=algorithm,
                     showSelected=testSet,
                     group=algorithm),
                 alpha=0.5,
                 fill=NA,
                 data=roc.not.cvx)+
    geom_point(aes(FPR, TPR, color=algorithm,
                   showSelected=testSet,
                   fill=train.type),
               shape=21,
               stroke=1,
               data=test.counts))
##print(viz$error)
##print(viz$auc)
##animint2dir(viz, "figure-test-error-dots")

dots <- ggplot()+
  geom_vline(aes(xintercept=min.auc),
             color="grey",
             data=best.auc)+
  geom_point(aes(mean.auc, algo.fac),
             shape=1,
             data=mean.auc)+
  geom_segment(aes(
    mean.auc+sd.auc, algo.fac, xend=mean.auc-sd.auc, yend=algo.fac),
               data=mean.auc)+
  facet_grid(. ~ set.name, labeller=function(df){
    df$set.name <- gsub("_", "\n", df$set.name)
    df
  }, scales="free_y", space="free_y")+
  scale_y_discrete("algorithm")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_x_continuous(
    "Test AUC (area under the Receiver Operating Characteristic curve)",
    breaks=c(0.6, 0.8, 1),
    labels=c("0.6", "0.8", "1"))
pdf("figure-test-error-mean.pdf", h=2, w=8)
print(dots)
dev.off()

dots <- ggplot()+
  geom_vline(aes(xintercept=min.auc),
             color="grey",
             data=best.auc)+
  geom_point(aes(mean.auc, algo.fac),
             alpha=0.3,
             size=3,
             data=mean.auc)+
  geom_point(aes(auc, algo.fac),
             data=auc, pch=1)+
  facet_grid(. ~ set.name, labeller=function(df){
    df$set.name <- gsub("_", "\n", df$set.name)
    df$set.name <- paste0(
      ifelse(
        grepl("H3K36me3", df$set.name),
        "Broad", "Sharp"),
      "\n",
      df$set.name)
    df
  }, scales="free")+
  scale_y_discrete("algorithm")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_x_continuous(
    "Test AUC (larger values indicate more accurate peak detection)",
    breaks=c(0.6, 0.8, 1),
    labels=c("0.6", "0.8", "1"))
dots <- ggplot()+
  geom_point(aes(med.auc, algo.fac),
             shape="|",
             size=3,
             data=mean.auc)+
  geom_segment(aes(q25.auc, algo.fac,
                   xend=q75.auc, yend=algo.fac),
               data=mean.auc)+
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
  scale_y_discrete("algorithm")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_x_continuous(
    "Test AUC (larger values indicate more accurate peak detection)",
    breaks=c(0.6, 0.8, 1),
    labels=c("0.6", "0.8", "1"))
pdf("figure-test-error-dots.pdf", h=2.2, w=8)
print(dots)
dev.off()

