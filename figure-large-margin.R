source("packages.R")

### Define what chunk and peaks to download.
set.name <- "H3K4me3_TDH_immune"
chunk.id <- 5

chunk.name <- paste0(set.name, "/", chunk.id)
regions.file <- file.path("data", set.name, chunk.id, "regions.RData")
load(regions.file)
counts.file <- file.path("data", set.name, chunk.id, "counts.RData")
load(counts.file)
model.file <- file.path("data", set.name, chunk.id, "dp.model.RData")
load(model.file)
region.list <- split(regions, regions$sample.id)
sample.ids <- c("McGill0322", "McGill0091", "McGill0002", "McGill0004")

chunk.peak.list <- list()
chunk.region.list <- list()
chunk.loss.list <- list()
for(sample.id in sample.ids){
  model.info <- dp.model[[sample.id]]
  chunk.loss.list[[sample.id]] <- data.table(
    sample.id, model.info$error)
  sample.regions <- region.list[[sample.id]]
  for(n.peaks in names(model.info$peaks)){
    peak.df <- model.info$peaks[[n.peaks]]
    region.df <- PeakErrorChrom(peak.df, sample.regions)
    if(nrow(peak.df)){
      chunk.peak.list[[paste(sample.id, n.peaks)]] <- data.table(
        sample.id, n.peaks, peak.df)
    }
    chunk.region.list[[paste(sample.id, n.peaks)]] <- data.table(
      sample.id, n.peaks, region.df)
  }
}
chunk.region <- do.call(rbind, chunk.region.list)
chunk.peak <- do.call(rbind, chunk.peak.list)
chunk.counts <- data.table(counts)[sample.id %in% sample.ids,]
chunk.loss <- do.call(rbind, chunk.loss.list)

error.counts <- chunk.region[, list(
  incorrect.labels=sum(fp+fn)
), by=.(sample.id, n.peaks)]

loss.list <- split(chunk.loss,
                   chunk.loss$sample.id,
                   drop=TRUE)
err.list <- split(error.counts, error.counts$sample.id, drop=TRUE)
signal.list <- split(chunk.counts,
                     chunk.counts$sample.id, drop=TRUE)
exact.dfs.list <- list()
intervals.list <- list()
for(sample.id in names(loss.list)){
  one <- loss.list[[sample.id]]
  exact.df <- with(one, exactModelSelection(error, peaks, peaks))
  err <- err.list[[sample.id]]
  setkey(err, n.peaks)
  signal <- signal.list[[sample.id]]
  exact.df$errors <-
    err[as.character(exact.df$model.complexity), ]$incorrect.labels
  indices <- with(exact.df, {
    largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
  })
  meta <- data.frame(sample.id, log.max.count=log(max(signal$coverage)))
  exact.dfs.list[[sample.id]] <- data.table(meta, exact.df)
  intervals.list[[sample.id]] <- data.table(
    meta,
    min.log.lambda=exact.df$min.log.lambda[indices$start],
    max.log.lambda=exact.df$max.log.lambda[indices$end])
}
exact.dfs <- do.call(rbind, exact.dfs.list)
intervals <- do.call(rbind, intervals.list)

what.peaks <- "peaks"
what.error <- "errors"
blank.df <- data.frame(x=10, what=what.peaks, sample.id="McGill0002")
exact.peaks <- data.frame(exact.dfs, what=what.peaks)
exact.error <- data.frame(exact.dfs, what=what.error)
funplot <- 
ggplot()+
  geom_segment(aes(min.log.lambda, model.complexity,
                   xend=max.log.lambda, yend=model.complexity),
               data=exact.peaks, size=2)+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=exact.error, size=2)+
  geom_blank(aes(x, 0), data=blank.df)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ sample.id, scales="free_y", space="free_y")+
  scale_x_continuous("$\\log \\lambda$", limits=c(7, 14),
                     breaks=seq(8, 12, by=2))+
  scale_y_continuous("", breaks=0:9)
print(funplot)

zero.peaks <- subset(exact.peaks, errors==0)
z.error <- subset(exact.error, errors==0)
tplot <- 
ggplot()+
  geom_segment(aes(min.log.lambda, model.complexity,
                   xend=max.log.lambda, yend=model.complexity),
               data=zero.peaks, size=3, color="green")+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=z.error, size=3, color="green")+
  geom_segment(aes(min.log.lambda, model.complexity,
                   xend=max.log.lambda, yend=model.complexity),
               data=exact.peaks, size=2)+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=exact.error, size=2)+
  geom_blank(aes(x, 0), data=blank.df)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ sample.id, scales="free", space="free")+
  scale_x_continuous("$\\log \\lambda$", limits=c(7, 14),
                     breaks=seq(8, 12, by=2))+
  scale_y_continuous("", breaks=0:9)
print(tplot)

what.intervals <- data.frame(intervals, what=what.error)
what.intervals$log.lambda <- what.intervals$max.log.lambda
left.intervals <- subset(what.intervals, is.finite(min.log.lambda))
left.intervals$log.lambda <- left.intervals$min.log.lambda
target.text <-
  rbind(data.frame(label="$\\overline L_i$", what.intervals, hjust=1),
        data.frame(label="$\\underline L_i$", left.intervals, hjust=0))
tsize <- 3
fun2 <- tplot+
  geom_text(aes(log.lambda, 1.5, label=label),
            data=target.text)+
  geom_point(aes(min.log.lambda, 0), shape=21, fill="white",
             data=left.intervals,
             size=tsize)+
  geom_point(aes(max.log.lambda, 0), shape=21, fill="black",
             data=what.intervals, size=tsize)

max.margin <- function
### Support vector interval regression for separable data. The idea is
### that we first normalize the feature matrix, giving normalized
### features x_i in R^p. Then we use a linear function f(x_i) = w'x_i
### + b to predict a log(lambda) that falls between all the log limits
### L_i^left and L_i^right and maximizes the margin. So the
### optimization problem is: max_{M,f} subject to, for all finite
### limits, L_i^right - f(x_i) >= M and f(x_i) - L_i^left >= M. Since
### we assume f is linear the problem becomes min_{w,b,M} -M subject
### to -M - w'x - b >= -L_i^right and -M + w'x + b >= L_i^left. We
### call M margin, w weights, b intercept.
(features,
### Matrix n x p of inputs: n signals, each with p features. We will
### scale these internally.
 limits,
### Matrix n x 2 of output lambda. Each row corresponds to the lower
### and upper bound of an interval on the log(lambda) which is optimal
### with respect to annotation error. Lower bound can be -Inf and
### upper bound can be Inf, which correspond to zero asymptotic
### cost. 
 verbose=0,
 ...
### ignored.
 ){
  ## reality checks.
  stopifnot(nrow(features)==nrow(limits))
  if(ncol(limits)!=2){
    cat("str(limits)=\n")
    str(limits)
    stop("limits should be a 2-column matrix")
  }
  stopifnot(is.matrix(features))
  
  ## check if there are any flat error curves, which have no limits.
  has.limits <- apply(is.finite(limits),1,any)
  ## we train the model on this subset.
  some.limits <- limits[has.limits,]
  some.features <- features[has.limits,,drop=FALSE]

  scaled <- scale(some.features)
  mu <- attr(scaled,"scaled:center")
  sigma <- attr(scaled,"scaled:scale")

  n <- nrow(scaled)
  p <- ncol(scaled)
  vars <- make.ids(margin=1,intercept=1,weights=p)
  constraints <- list(vars$margin*1 >= 0)
  for(i in 1:n){
    if(verbose >= 1)cat(sprintf("example constraints %5d / %5d",i,n))

    left <- some.limits[i,1]
    if(is.finite(left)){
      ivars <- with(vars,{
        intercept * 1 + sum(weights)*scaled[i,] + margin*-1
      })
      constraints <- c(constraints,list(ivars >= left))
    }

    right <- some.limits[i,2]
    if(is.finite(right)){
      ivars <- with(vars,{
        intercept * -1 + sum(weights)*scaled[i,]*-1 +margin*-1
      })
      constraints <- c(constraints,list(ivars >=  - right))
    }

    if(verbose >= 1)cat("\n")

  }
  const.info <- standard.form.constraints(constraints,vars)
  n.vars <- length(unlist(vars))
  Dvec <- rep(1e-10,n.vars)
  D <- diag(Dvec)
  d <- rep(0,n.vars)
  d[vars$margin] <- 1
  if(verbose >= 1)cat(sprintf("solving for %d variables and %d constraints... ",
              n.vars,length(constraints)))
  sol <- solve.QP(D,d,const.info$A,const.info$b0)
  if(verbose >= 1)cat("solved!\n")
  sol$mu <- mu
  sol$sigma <- sigma
  sol$scaled <- scaled
  sol$log.limits <- some.limits
  sol$features <- some.features
  sol$weights <- sol$solution[vars$weights]
  sol$intercept <- sol$solution[vars$intercept]
  sol$margin <- sol$solution[vars$margin]
  ## this function will be applied to new data before applying the
  ## model.
  sol$normalize <- function(X){
    mu.mat <- matrix(mu,nrow(X),ncol(X),byrow=TRUE)
    s.mat <- matrix(sigma,nrow(X),ncol(X),byrow=TRUE)
    (X-mu.mat)/s.mat
  }
  sol$f <- function(x){
    sum(x*sol$weights)+sol$intercept
  }
  sol$predict <- function(X){
    stopifnot(is.matrix(X))
    X.norm <- sol$normalize(X)
    weights.mat <- matrix(sol$weights,nrow(X),ncol(X),byrow=TRUE)
    L.hat <- rowSums(X.norm * weights.mat) + sol$intercept
    L.hat
  }
  sol$L.pred <- apply(scaled,1,sol$f)
  sol$lambda.pred <- sol$predict(features)
  sol
### List of solver results. For a feature matrix X with p columns, you
### can use list$predict(X) to get model estimates of log(lambda).
}

fit <- with(intervals, {
  max.margin(cbind(log.max.count), cbind(min.log.lambda, max.log.lambda))
})

zero.error <- subset(exact.dfs, errors==0)

intervals$predicted <- fit$L.pred
intervals$mid.log.lambda <- with(intervals, (max.log.lambda+min.log.lambda)/2)
prediction.in.mid <- with(intervals, abs(mid.log.lambda-predicted) < 1e-6)
## max margin line should be in the middle of one interval.
stopifnot(any(prediction.in.mid)) 
min.x <- min(intervals$log.max.count)
max.x <- max(intervals$log.max.count)
small.slope <- (9-11)/(min.x-max.x)
intervals$bad.pred <- small.slope*(intervals$log.max.count - max.x) + 11
intervals$right.margin <- with(intervals, max.log.lambda-predicted)
intervals$left.margin <- with(intervals, predicted-min.log.lambda)
min.log.count <- min(intervals$log.max.count)
max.log.count <- max(intervals$log.max.count)
seg.dt <- data.table(
  min.log.lambda=9, min.log.count,
  max.log.lambda=c(11, 11.5), max.log.count)
seg.dt[, slope := (min.log.lambda-max.log.lambda)/(min.log.count-max.log.count)]
seg.dt[, intercept := min.log.lambda-slope*min.log.count]
reg.dt <- rbind(
  seg.dt[, .(slope, intercept)],
  data.frame(slope=0, intercept=8.5))
count.grid <- c(3, 7)
penalty.grid.list <- list()
for(reg.i in 1:nrow(reg.df)){
  r <- reg.df[reg.i, ]
  penalty.grid.list[[reg.i]] <- data.table(
    count.grid, log.lambda=r$slope * count.grid + r$intercept, reg.i)
}
penalty.grid <- do.call(rbind, penalty.grid.list)

extreme <- subset(intervals, log.max.count %in% range(log.max.count))
extreme$max.log.lambda[1] <- 11
extreme$min.log.lambda[2] <- extreme$predicted[2]

## Plot the model that was selected by the large margin model.
setkey(intervals, sample.id)
exact.ids <- paste(exact.dfs$sample.id)
exact.dfs[, predicted := intervals[exact.ids,]$predicted ]
selected <- exact.dfs[min.log.lambda < predicted & predicted < max.log.lambda,]

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

sample.max.df <- chunk.counts[, list(
  count=max(coverage)
  ), by=sample.id]
sample.max <- sample.max.df$count
names(sample.max) <- as.character(sample.max.df$sample.id)

disp.counts <- chunk.counts[, {
  base <- seq(min(chromStart), max(chromEnd), l=500)
  data.table(base, count=approx(chromStart, coverage, base)$y)
}, by=sample.id]

selectedPlot <- 
ggplot()+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                data=chunk.region[n.peaks==0,],
                color="grey",
                alpha=0.5)+
  geom_line(aes(base/1e3, count),
            data=disp.counts, color="grey50")+
  ## geom_step(aes(chromStart/1e3, coverage),
  ##           data=chunk.counts, color="grey50")+
  ## geom_point(aes(chromStart/1e3, 0),
  ##            data=profile.list$peaks,
  ##            pch=1, size=2, color="deepskyblue")+
  ## geom_text(aes(118120, y.mid, label=label),
  ##           data=compare.labels, hjust=1, size=3)+
  ## geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
  ##               ymin=y.min, ymax=y.max,
  ##               linetype=status),
  ##           data=compare.regions,
  ##           fill=NA, color="black", size=0.5)+
  ## scale_linetype_manual("error type",
  ##                       values=c(correct=0,
  ##                         "false negative"=3,
  ##                         "false positive"=1))+
  ## geom_segment(aes(chromStart/1e3, y.mid,
  ##                  xend=chromEnd/1e3, yend=y.mid),
  ##              data=compare.peaks, size=1.5, color="deepskyblue")+
  ## geom_segment(aes(chromStart/1e3, 0,
  ##                  xend=chromEnd/1e3, yend=0),
  ##              data=profile.list$peaks, size=1.5, color="#6A3D9A")+
  coord_cartesian(xlim=c(118090, 118125), expand=FALSE)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free")+
  scale_y_continuous("aligned read coverage",
                     breaks=function(limits){
                       floor(limits[2])
                     })+
  xlab("position on chr11 (kilo base pairs)")+
  scale_fill_manual("label", values=ann.colors,
                    breaks=names(ann.colors))
print(selectedPlot)

p <- 
ggplot()+
  geom_segment(aes(min.log.lambda, log.max.count,
                   xend=max.log.lambda, yend=log.max.count),
               data=intervals, size=1.5)+
  geom_segment(aes(min.log.lambda, log.max.count,
                   xend=max.log.lambda, yend=log.max.count),
               data=extreme, color="red")+
  ## geom_line(aes(predicted, log.max.count), data=intervals, color="blue")+
  ## geom_line(aes(bad.pred, log.max.count), data=intervals, color="blue")+
  ## geom_vline(aes(xintercept=8.3), color="blue")+
  geom_line(aes(log.lambda, count.grid, group=reg.i),
            data=penalty.grid, color="blue")+
  coord_cartesian(ylim=c(4.2, 5.7))+
  ## geom_segment(aes(min.log.lambda, min.log.count,
  ##                  xend=max.log.lambda, yend=max.log.count),
  ##              data=seg.df, color="blue")+
  geom_text(aes(8.6, 5.3, label="1 error\nconstant"))+
  geom_text(aes(10.25, 5.3, label="0 errors\nsmall margin"))+
  geom_text(aes(11.5, 5.3, label="0 errors\nlarge margin"))+
  geom_point(aes(min.log.lambda, log.max.count), data=zero.error, size=5,pch="|")+
  geom_point(aes(max.log.lambda, log.max.count), data=zero.error, size=5, pch="|")+
  geom_text(aes((min.log.lambda + max.log.lambda)/2, log.max.count,
                label=sprintf("%d peak%s", model.complexity,
                  ifelse(model.complexity==1, "", "s")),
                hjust=ifelse(is.finite(min.log.lambda), 0.5, 0)),
            data=zero.error, vjust=-0.5, size=3)+
  geom_text(aes(max.log.lambda, log.max.count, label=sample.id),
            data=intervals, hjust=0)+
  ggtitle("max margin interval regression, margin in red")
print(p)
 
## Plot the max margin regression line.
zero.error[, `:=`(
  label=sprintf("%d peak%s", model.complexity,
                ifelse(model.complexity==1, "", "s")),
  label.penalty=(min.log.lambda + max.log.lambda)/2)]
text.df <- data.frame(zero.error)
rownames(text.df) <- with(text.df, paste(sample.id, model.complexity))
text.df["McGill0004 2", "label.penalty"] <- 8.75
text.df["McGill0091 1", "label.penalty"] <- 10
text.df["McGill0002 2", "label.penalty"] <- 12
tsize <- 2.5
modelsPlot <- 
  ggplot()+
  theme_bw()+
  geom_segment(aes(log.max.count, min.log.lambda, 
                   yend=max.log.lambda, xend=log.max.count),
               data=intervals, size=1.5)+
  geom_point(aes(log.max.count, min.log.lambda),
             data=subset(zero.error, is.finite(min.log.lambda)),
             size=tsize, pch=1)+
  geom_point(aes(log.max.count, max.log.lambda),
             data=zero.error,
             size=tsize, pch=1)+
  geom_point(aes(log.max.count, min.log.lambda),
             data=subset(intervals, is.finite(min.log.lambda)),
             size=tsize, shape=21, fill="white")+
  geom_point(aes(log.max.count, max.log.lambda),
             data=intervals,
             size=tsize, shape=21, fill="black")+
  geom_text(aes(log.max.count +
                ifelse(log.max.count==max(log.max.count), -1, 1)*0.03,
                label.penalty,
                label=label,
                vjust=ifelse(is.finite(min.log.lambda), 0.5, -0.5),
                hjust=ifelse(log.max.count==max(log.max.count), 1, 0)),
            data=text.df, size=3)+
  geom_text(aes(log.max.count, max.log.lambda, label=sample.id,
                hjust=ifelse(log.max.count==min(log.max.count), 0,
                  ifelse(log.max.count==max(log.max.count), 1, 0.5))),
            data=intervals, vjust=-0.5, size=3)+
  coord_cartesian(ylim=c(7.5, 13), xlim=c(3.6, 6.3))+
  scale_x_continuous("penalty $\\log\\lambda_i$", 
                     minor_breaks=NULL)+
  scale_y_continuous("feature $x_i = \\log\\max\\mathbf y_i$",
                     minor_breaks=NULL)+
  geom_segment(aes(log.max.count, min.log.lambda, 
                   yend=predicted, xend=log.max.count),
               data=intervals["McGill0002",], color="red")+
  geom_line(aes(count.grid, log.lambda, group=reg.i),
            data=penalty.grid, color="blue")+
  geom_text(aes(5.8, 12.5, label="0 errors\nlarge margin"),
            hjust=0, vjust=0, color="blue", size=3)+
  geom_text(aes(5.8, 11, label="0 errors\nsmall margin"),
            hjust=0, vjust=1, color="blue", size=3)+
  geom_text(aes(5.8, 9, label="1 error\nconstant"),
            hjust=0, color="blue", size=3)
print(modelsPlot)

ggplot()+
  geom_segment(aes(model.complexity, min.log.lambda, 
                   yend=max.log.lambda, xend=model.complexity),
               data=zero.peaks, size=3, color="green")+
  geom_segment(aes(errors, min.log.lambda, 
                   yend=max.log.lambda, xend=errors),
               data=z.error, size=3, color="green")+
  geom_segment(aes(model.complexity, min.log.lambda, 
                   yend=max.log.lambda, xend=model.complexity),
               data=exact.peaks, size=2)+
  geom_segment(aes(errors, min.log.lambda, 
                   yend=max.log.lambda, xend=errors),
               data=exact.error, size=2)+
  geom_blank(aes(0, x), data=blank.df)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sample.id ~ what, scales="free", space="free")+
  scale_y_continuous("$\\log \\lambda$", limits=c(7, 14),
                     breaks=seq(8, 12, by=2))+
  scale_x_continuous("", breaks=0:9)

min.log.pen <- min(exact.peaks$max.log.lambda)
max.log.pen <- max(exact.peaks$min.log.lambda)
notInf <- function(x)ifelse(x==Inf, max.log.pen+1, ifelse(x==-Inf, min.log.pen-1, x))
viz <- list(
  coverage=ggplot()+
    ggtitle("ChIP-seq data and peaks")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    ##theme_animint(width=1000)+
    facet_grid(sample.id ~ ., scales="free")+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    xlab("position on chr11 (kilo base pairs)")+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  data=chunk.region[n.peaks==0,],
                  color="grey",
                  alpha=0.5)+
    scale_linetype_manual("error type",
                          values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
    geom_tallrect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3,
      showSelected.variable=paste0(sample.id, "peaks"),
      showSelected.value=n.peaks,
      linetype=status),
      data=chunk.region,
      color="black",
      fill=NA,
      size=1.5)+
    geom_line(aes(base/1e3, count),
              data=disp.counts, color="grey50")+
    geom_segment(aes(
      chromStart/1e3, 0,
      xend=chromEnd/1e3, yend=0,
      key=paste(sample.id, chromStart, chromEnd),
      showSelected.variable=paste0(sample.id, "peaks"),
      showSelected.value=peaks),
      data=chunk.peak,
      color="deepskyblue",
      size=4)+
    geom_point(aes(
      chromStart/1e3, 0,
      key=paste(sample.id, chromStart),
      showSelected.variable=paste0(sample.id, "peaks"),
      showSelected.value=peaks),
      data=chunk.peak,
      color="black",
      fill="deepskyblue")+
    scale_fill_manual("label", values=ann.colors,
                      breaks=names(ann.colors)),
  penalty=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(. ~ what, scales="free")+
    geom_text(aes(
      feature, log.penalty,
      label=label, vjust=vjust, hjust=hjust),
      data=data.table(
        what="regression",
        feature=5.8,
        log.penalty=c(12.5, 11, 9),
        hjust=c(0.5,0,0),
        label=c(
          "0 errors\nlarge margin",
          "0 errors\nsmall margin",
          "1 error\nconstant"),
        vjust=c(0,1,0.5)),
      color="blue",
      size=3)+
    geom_segment(aes(
      log.max.count, min.log.lambda,
      clickSelects=sample.id,
      yend=max.log.lambda, xend=log.max.count),
      data=data.table(intervals, what="regression"),
      size=4,
      alpha=0.7,
      color="green")+
    geom_point(aes(
      log.max.count, min.log.lambda,
      clickSelects=sample.id),
      data=data.table(zero.error, what="regression")[is.finite(min.log.lambda),],
      size=tsize, pch=1)+
    geom_point(aes(
      log.max.count, max.log.lambda,
      clickSelects=sample.id),
      data=data.table(zero.error, what="regression"),
      size=tsize,
      pch=1)+
    geom_point(aes(
      log.max.count, min.log.lambda,
      clickSelects=sample.id),
      data=data.table(intervals, what="regression")[is.finite(min.log.lambda),],
      size=tsize,
      shape=21,
      fill="white")+
    geom_point(aes(
      log.max.count, max.log.lambda,
      clickSelects=sample.id),
      data=data.table(intervals, what="regression"),
      size=tsize,
      shape=21,
      fill="black")+
    ## geom_text(aes(
    ##   log.max.count +
    ##     ifelse(log.max.count==max(log.max.count), -1, 1)*0.03,
    ##   label.penalty,
    ##   label=label,
    ##   vjust=ifelse(is.finite(min.log.lambda), 0.5, -0.5),
    ##   hjust=ifelse(log.max.count==max(log.max.count), 1, 0)),
    ##   data=data.table(text.df, what="regression"))+
    geom_text(aes(
      log.max.count, max.log.lambda, label=sample.id,
      clickSelects=sample.id,
      hjust=ifelse(log.max.count==min(log.max.count), 0,
                   ifelse(log.max.count==max(log.max.count), 1, 0.5))),
      data=data.table(intervals, what="regression"),
      vjust=-0.5)+
    geom_segment(aes(
      log.max.count, min.log.lambda, 
      yend=predicted, xend=log.max.count),
      data=data.table(intervals["McGill0002",], what="regression"),
      color="red")+
    geom_line(aes(
      count.grid, log.lambda, group=reg.i),
      data=data.table(penalty.grid, what="regression"),
      size=1,
      color="blue")+
    geom_text(aes(
      x, y, label=label),
      data=data.table(
        what="regression",
        x=5, y=6,
        label="log.max.count"))+
    ##
    geom_segment(aes(model.complexity, min.log.lambda,
                     showSelected=sample.id,
                     yend=max.log.lambda, xend=model.complexity),
                 data=zero.peaks, size=4, color="green")+
    geom_segment(aes(errors, min.log.lambda, 
                     showSelected=sample.id,
                     yend=max.log.lambda, xend=errors),
                 data=z.error, size=4, color="green")+
    geom_segment(aes(
      model.complexity,
      notInf(min.log.lambda), 
      showSelected=sample.id,
      key=peaks,
      yend=notInf(max.log.lambda),
      xend=model.complexity),
      data=exact.peaks, size=2)+
    geom_segment(aes(
      errors, notInf(min.log.lambda), 
      showSelected=sample.id,
      key=peaks,
      yend=notInf(max.log.lambda), xend=errors),
      data=exact.error, size=2)+
    geom_widerect(aes(
      ymin=notInf(min.log.lambda), 
      showSelected=sample.id,
      key=peaks,
      clickSelects.variable=paste0(sample.id, "peaks"),
      clickSelects.value=peaks,
      ymax=notInf(max.log.lambda)),
      alpha=0.2,
      data=exact.peaks)+
    geom_widerect(aes(
      ymin=notInf(min.log.lambda), 
      showSelected=sample.id,
      key=peaks,
      clickSelects.variable=paste0(sample.id, "peaks"),
      clickSelects.value=peaks,
      ymax=notInf(max.log.lambda)),
      alpha=0.2,
      data=exact.error, size=2)+
    xlab("")+
    scale_y_continuous(breaks=0:9)
)
animint2dir(viz, "figure-large-margin")
