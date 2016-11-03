source("packages.R")

load("../PeakSeg-paper/dp.peaks.error.RData")
##load("../PeakSeg-paper/dp.peaks.RData")
load("PDPA.peaks.error.RData")
load("Segmentor.peaks.error.RData")
load("PDPA.peaks.RData")
load("macs.peaks.error.RData")
load("hmcan.broad.peaks.error.RData")
load("PDPA.cDPA.compare.RData")

macs.peaks.error$algorithm <- "macs"
hmcan.broad.peaks.error$algorithm <- "hmcan.broad"
pdpa <- data.table(PDPA.peaks.error)
names(pdpa)[3] <- "param.name"
pdpa$algorithm <- "coseg"
Seg <- data.table(Segmentor.peaks.error)
names(Seg)[3] <- "param.name"
Seg$algorithm <- "Segmentor"
error.regions.list <- list(
  macs=macs.peaks.error,
  hmcan.broad=hmcan.broad.peaks.error,
  coseg=pdpa,
  Segmentor=Seg)
for(chunk.name in names(dp.peaks.error)){
  one.chunk <- dp.peaks.error[[chunk.name]]
  error.regions.list[[chunk.name]] <- data.table(
    chunk.name, do.call(rbind, one.chunk),
    algorithm="PeakSegDP")
}
error.regions <- do.call(rbind, error.regions.list)

error.counts <- error.regions[, list(
  errors=sum(fp+fn),
  fn=sum(fn),
  possible.fn=sum(possible.tp),
  fp=sum(fp),
  possible.fp=sum(possible.fp),
  regions=.N
  ), by=.(algorithm, chunk.name, sample.id, param.name)]

min.train.error <-
  error.counts[, .SD[which.min(errors),], by=.(algorithm, chunk.name, sample.id)]

total.train.error <- min.train.error[, list(
  errors=sum(errors),
  fp=sum(fp),
  fn=sum(fn),
  problems=.N,
  labels=sum(regions),
  possible.fp=sum(possible.fp),
  possible.fn=sum(possible.fn)
  ), by=.(algorithm)]

## For each algorithm, how many problems had all ten models? how many
## models were computed and feasible in all?
feasible <- PDPA.cDPA.compare$all.loss[, list(
  PeakSegDP=sum(!is.na(dp.fwd)),
  coseg=sum(PDPA.feasible),
  Segmentor=sum(Seg.feasible),
  rows=.N
  ), by=.(chunk.name, sample.id)]
stopifnot(all(feasible$rows==10))
stopifnot(nrow(feasible)==2752)
molt <- melt(feasible, measure.vars=c("PeakSegDP", "coseg", "Segmentor"))
model.counts <- molt[, list(
  models=sum(value),
  ten=sum(value==10)
  ), by=.(variable)]

train.ord <- total.train.error[order(errors),]
total.row <- train.ord[1, .(errors=labels, fp=possible.fp, fn=possible.fn, models=27520, problems=2752)]
setkey(model.counts, variable)
only.algos <- cbind(
  train.ord[, .(errors, fp, fn)],
  model.counts[paste(train.ord$algorithm), .(models, problems=ten)])
no.names <- as.matrix(rbind(only.algos, total.row))
rownames(no.names) <- c(train.ord$algorithm, "possible")
library(xtable)
xt <- xtable(no.names, digits=0)
print(xt, file="table-min-train-error.tex", row.names=FALSE, floating=FALSE)

train.error.wide <-
  dcast(min.train.error, chunk.name + sample.id ~ algorithm, value.var="errors")
train.error.wide[PeakSegDP < coseg & coseg < Segmentor,]
## nice example to show the difference between algos.
##  8:  H3K4me3_TDH_immune/4 McGill0005         1         4     2
train.error.wide[coseg==0 & Segmentor==6,]
train.error.wide[coseg==5 & PeakSegDP==1,]
train.error.wide[coseg==3 & macs==0,]
train.error.wide[coseg==0 & macs==7,]

setkey(train.error.wide, chunk.name, sample.id)
show.dt <- train.error.wide[J(
  c("H3K4me3_PGP_immune/15", "H3K4me3_PGP_immune/14", "H3K4me3_TDH_immune/4",
    "H3K4me3_TDH_immune/3", "H3K36me3_TDH_immune/2"),
  c("McGill0079", "McGill0095", "McGill0005", "McGill0091", "McGill0009")),]

table(train.error.wide[, coseg-PeakSegDP])
table(train.error.wide[, coseg-Segmentor])

train.error.wide[, table(PeakSegDP, coseg)]
train.error.wide[, table(Segmentor, coseg)]
train.error.wide[, table(macs, coseg)]

abline.dt <- data.table(slope=1, intercept=0)
molt.list <- list()
for(other.name in c("PeakSegDP", "Segmentor", "macs")){
  table.args <- list()
  for(col.name in c(other.name, "coseg")){
    table.args[[col.name]] <- train.error.wide[[col.name]]
  }
  tab <- do.call(table, table.args)
  molt <- melt(tab, value.name="problems")
  molt$log10.problems <- log10(molt$problems)
  gg <- ggplot()+
    xlab(paste(
      "incorrect labels in best",
      other.name,
      "model"))+
    ylab(paste(
      "incorrect labels in best",
      "coseg",
      "model"))+
    theme_bw()+
    geom_abline(aes(slope=slope, intercept=intercept), data=abline.dt, color="grey")+
    coord_equal()+
    scale_fill_gradient(low="grey90", high=scales::muted("red"), na.value="white")+
    geom_tile(aes_string(x=other.name, y="coseg", fill="log10.problems"), data=molt)+
    geom_text(aes_string(x=other.name, y="coseg", label="problems"),
              data=subset(molt, 0 < problems))
  pdf(sprintf("figure-min-train-error-%s.pdf", other.name), h=5)
  print(gg)
  dev.off()
  names(molt)[1] <- "other.errors"
  molt.list[[other.name]] <- data.table(other.name, molt)
}
prob.counts <- do.call(rbind, molt.list)

some.counts <- subset(prob.counts, 0 < problems & other.name != "macs")
some.counts[, winner := ifelse(other.errors==coseg, "both", ifelse(other.errors<coseg, "other", "coseg"))]
some.counts[, list(problems=sum(problems)), by=.(other.name, winner)]
gg.some <- ggplot()+
  geom_abline(aes(slope=slope, intercept=intercept), data=abline.dt, color="grey")+
  xlab(paste(
    "incorrect labels in best",
    "competing",
    "model"))+
  ylab(paste(
    "incorrect labels in best",
    "coseg",
    "model"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ other.name)+
  coord_equal()+
  scale_fill_gradient(low="grey90", high=scales::muted("red"), na.value="white")+
  geom_tile(aes(other.errors, coseg, fill=log10.problems), data=some.counts)+
  geom_text(aes(other.errors, coseg, label=problems),
            size=3,
            data=some.counts)
print(gg.some)
pdf("figure-min-train-error-some.pdf", h=2.9)
print(gg.some)
dev.off()

nonzero.counts <- subset(prob.counts, 0 < problems)
gg.counts <- ggplot()+
  geom_abline(aes(slope=slope, intercept=intercept), data=abline.dt, color="grey")+
  xlab(paste(
    "incorrect labels in best",
    "competing",
    "model"))+
  ylab(paste(
    "incorrect labels in best",
    "coseg",
    "model"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ other.name)+
  coord_equal()+
  scale_fill_gradient(low="grey90", high=scales::muted("red"), na.value="white")+
  geom_tile(aes(other.errors, coseg, fill=log10.problems), data=nonzero.counts)+
  geom_text(aes(other.errors, coseg, label=problems),
            size=3,
            data=nonzero.counts)
print(gg.counts)

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
setkey(error.regions, chunk.name, sample.id)
setkey(error.counts, chunk.name, sample.id)
for(show.row.i in 1:nrow(show.dt)){
  show.row <- show.dt[show.row.i,]
  counts.file <- paste0("data/", show.row$chunk.name, "/counts.RData")
  load(counts.file)
  counts.by.sample <- split(counts, counts$sample.id)
  sample.counts <- data.table(counts.by.sample[[paste(show.row$sample.id)]])
  sample.regions <- error.regions[show.row]
  sample.error <- error.counts[show.row]
  print(dcast(sample.error, param.name ~ algorithm, value.var="errors"))
  load(sub("counts", "Segmentor.model", counts.file))
  load(sub("counts", "dp.model", counts.file))
  sample.peaks.raw <- list(
    coseg=PDPA.peaks[[paste(show.row$chunk.name)]][[paste(show.row$sample.id)]],
    Segmentor=Segmentor.model[[paste(show.row$sample.id)]]$peaks,
    PeakSegDP=dp.model[[paste(show.row$sample.id)]]$peaks)
  col.name.vec <- c("chromStart", "chromEnd", "peaks")
  sample.peaks.list <- list()
  for(algorithm in names(sample.peaks.raw)){
    algo.peaks <- sample.peaks.raw[[algorithm]]
    sample.peaks.list[[algorithm]] <-
      data.table(algorithm, do.call(rbind, algo.peaks)[, col.name.vec])
  }
  sample.peaks <- do.call(rbind, sample.peaks.list)

  gg.counts.prob <- gg.counts+
    geom_tile(aes(PeakSegDP, coseg),
              fill=NA,
              color="black",
              size=2,
              data=data.table(show.row, other.name="PeakSegDP"))+
    geom_tile(aes(macs, coseg),
              fill=NA,
              color="black",
              size=2,
              data=data.table(show.row, other.name="macs"))+
    geom_tile(aes(Segmentor, coseg),
              fill=NA,
              color="black",
              size=2,
              data=data.table(show.row, other.name="Segmentor"))
  pdf(sprintf("figure-min-train-error-problem%d.pdf", show.row.i), w=10)
  print(gg.counts.prob)
  dev.off()

  max.coverage <- max(sample.counts$coverage)
  first.chromStart <- min(sample.regions$chromStart)
  last.chromEnd <- max(sample.regions$chromEnd)
  y.key <- c(coseg=1, Segmentor=2, PeakSegDP=3, macs=4)*max.coverage*-0.15
  h <- abs(diff(y.key)[1]/3)

  for(peaks.str in paste(0:6)){
    png.name <- sprintf("figure-min-train-error-problem%d-%speaks.png", show.row.i, peaks.str)
    sample.gg <- ggplot()+
      theme_bw()+
      ggtitle(show.row[, paste(peaks.str, "peak models of problem", chunk.name, sample.id)])+
      xlab("position on chromosome (kb = kilo bases)")+
      ylab("aligned sequence reads")+
      scale_fill_manual(values=ann.colors)+
      geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
                    color="grey",
                    alpha=0.5,
                    data=sample.regions[algorithm=="coseg" & param.name==0,])+
      geom_step(aes(chromStart/1e3, coverage),
                color="grey50",
                data=sample.counts)+
      scale_linetype_manual("error type",
                            values=c(correct=0,
                              "false negative"=3,
                              "false positive"=1),
                            limits=c("correct", "false negative", "false positive"))+
      geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    linetype=status,
                    ymin=y.key[algorithm]-h, ymax=y.key[algorithm]+h),
                size=0.75,
                color="black",
                fill=NA,
                data=sample.regions[algorithm != "macs" & param.name==peaks.str,])
    these.peaks <- sample.peaks[peaks==peaks.str,]
    if(nrow(these.peaks)){
      sample.gg <- sample.gg+
      geom_segment(aes(chromStart/1e3, y.key[algorithm],
                       xend=chromEnd/1e3, yend=y.key[algorithm]),
                   size=2,
                   color="deepskyblue",
                   data=these.peaks[algorithm != "macs",])
    }
    sample.gg <- sample.gg+
      geom_text(aes(first.chromStart/1e3, y.key[algorithm], label=algorithm),
                hjust=1,
                size=3.5,
                data=sample.error[algorithm != "macs" & param.name==peaks.str,])+
      geom_text(aes(
        last.chromEnd/1e3, y.key[algorithm],
        label=paste0(" ", errors, " error", ifelse(errors==1, "", "s"))
        ),
                hjust=0,
                size=3.5,
                data=sample.error[algorithm != "macs" & param.name==peaks.str,])
    print(png.name)
    png(png.name, 9, 6, res=100, units="in")
    print(sample.gg)
    dev.off()
  }

  load(sub("counts", "peaks/macs.trained", counts.file))
  sample.error.min <- sample.error[, .SD[which.min(errors),], by=.(algorithm)]
  setkey(sample.error.min, algorithm, param.name)
  macs.one.param <- peaks[[paste(sample.error.min["macs"]$param.name)]]
  macs.one.sample <- subset(macs.one.param, paste(show.row$sample.id)==sample.id)
  sample.peaks[, param.name := peaks]
  setkey(sample.peaks, algorithm, param.name)
  not.macs <- sample.error.min[algorithm!="macs",]
  best.peaks <- rbind(
    data.table(algorithm="macs", macs.one.sample[, c("chromStart", "chromEnd")]),
    sample.peaks[not.macs, .(algorithm, chromStart, chromEnd)])
  setkey(sample.regions, algorithm, param.name)
  best.regions <- sample.regions[sample.error.min]

  png.name <- sprintf("figure-min-train-error-problem%d-best.png", show.row.i)
  sample.gg <- ggplot()+
    theme_bw()+
    ggtitle(show.row[, paste("best models for problem", chunk.name, sample.id)])+
    xlab("position on chromosome (kb = kilo bases)")+
    ylab("aligned sequence reads")+
    scale_fill_manual(values=ann.colors)+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
                  color="grey",
                  alpha=0.5,
                  data=sample.regions[algorithm=="coseg" & param.name==0,])+
    geom_step(aes(chromStart/1e3, coverage),
              color="grey50",
              data=sample.counts)+
    scale_linetype_manual("error type",
                          values=c(correct=0,
                            "false negative"=3,
                            "false positive"=1))+
    geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                  linetype=status,
                  ymin=y.key[algorithm]-h, ymax=y.key[algorithm]+h),
              size=0.75,
              color="black",
              fill=NA,
              data=best.regions)
  if(nrow(best.peaks)){
    sample.gg <- sample.gg+
      geom_segment(aes(chromStart/1e3, y.key[algorithm],
                       xend=chromEnd/1e3, yend=y.key[algorithm]),
                   size=2,
                   color="deepskyblue",
                   data=best.peaks)
  }
  sample.gg <- sample.gg+
    geom_text(aes(
      first.chromStart/1e3, y.key[algorithm],
      label=sprintf(
        "%s\n%s=%s",
        algorithm,
        ifelse(algorithm=="macs", "log(q)", "peaks"),
        paste(param.name)
        )),
              hjust=1,
              size=3.5,
              data=sample.error.min)+
    geom_text(aes(
      last.chromEnd/1e3, y.key[algorithm],
      label=paste0(" ", errors, " error", ifelse(errors==1, "", "s"))
      ),
              hjust=0,
              size=3.5,
              data=sample.error.min)
  print(png.name)
  png(png.name, 9, 6, res=100, units="in")
  print(sample.gg)
  dev.off()

  if(show.row.i==4){
    png.name <- sprintf("figure-min-train-error-problem%d-best-zoom.png", show.row.i)
    zoom.gg <- sample.gg+
      coord_cartesian(xlim=c(5440, 5465))
    print(png.name)
    png(png.name, 9, 6, res=100, units="in")
    print(zoom.gg)
    dev.off()
  }
}

train.error.wide[order(PeakSegDP-coseg),]
sort(train.error.wide[PeakSegDP<coseg, table(chunk.name)])
train.error.wide[PeakSegDP<coseg & grepl("TDH", chunk.name),]
error.counts[chunk.name=="H3K4me3_PGP_immune/15" & sample.id=="McGill0079",]
error.counts[chunk.name=="H3K4me3_PGP_immune/15" & algorithm=="coseg",]
error.counts[chunk.name=="H3K4me3_TDH_immune/1" & sample.id=="McGill0011",]

pdf("figure-min-train-error.pdf")
print(gg.counts)
dev.off()
