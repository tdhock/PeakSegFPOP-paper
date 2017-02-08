source("packages.R")

load("../PeakSeg-paper/dp.timings.RData")
load("PDPA.timings.RData")
load("Segmentor.timings.RData")
(objs <- load("dp.timings.reverse.RData"))
stopifnot(nrow(Segmentor.timings)==2752)
stopifnot(nrow(PDPA.timings)==2752)
stopifnot(nrow(dp.timings)==2752)
stopifnot(nrow(dp.timings.reverse)==2752)


algo <- function(algorithm, ...){
  data.table(algorithm, ...)
}
CDPA.name <- "PeakSegDP\nO(N^2)"
PDPA.name <- "Segmentor\nO(N log N)"
GPDPA.name <- "coseg\nO(N log N)"
CDPA.name <- "CDPA\nO(n^2)"
PDPA.name <- "PDPA\nO(n log n)"
GPDPA.name <- "GPDPA\nO(n log n)"
all.timings <- rbind(
  algo(CDPA.name, dp.timings),
  algo(PDPA.name, Segmentor.timings),
  ## algo("DP.fwd", dp.timings),
  ## algo("DP.rev", dp.timings.reverse),
  algo(GPDPA.name, PDPA.timings))
all.timings[, algo.tex := sub(
  "\nO[(]", "\n$O(", sub("[)]$", ")$", sub("log", "\\\\log", algorithm)))]

gg.quad <- ggplot()+
  ylab("hours of computation time")+
  xlab("data points to segment N")+
  geom_point(aes(data, seconds/60/60),
             shape=1,
             data=dp.timings)
pdf("figure-PDPA-timings-dp.pdf", h=4)
print(gg.quad)
dev.off()

totals <- all.timings[, list(seconds=sum(seconds)), by=algorithm]
totals[, minutes := seconds/60]
totals[, hours := minutes/60]
  
gg.linear <- ggplot()+
  ylab("hours")+
  geom_point(aes(data, seconds/60/60, color=algorithm),
             shape=1,
             data=all.timings)
print(gg.linear)

n.edge <- 41
edge.vec <-
  all.timings[, 10^seq(log10(min(data)), log10(max(data)+1), l=n.edge)]
window.dt <- data.table(min.data=edge.vec[-n.edge], max.data=edge.vec[-1])
setkey(window.dt, min.data, max.data)
all.timings[, data0 := data]
setkey(all.timings, data, data0)
over.dt <- foverlaps(all.timings, window.dt, nomatch=0L)
stopifnot(nrow(over.dt)==nrow(all.timings))
window.seconds <- over.dt[, list(
  median=median(seconds),
  problems=.N,
  min=min(seconds),
  max=max(seconds)
  ), by=.(algo.tex, min.data, max.data)]
stopifnot(0 < window.seconds$problems)
lab.df <- data.frame(
  seconds=c(1, 60, 60*60),
  label=c("1 second", "1 minute", "1 hour"))
gg.log <- ggplot()+
  theme_bw()+
  geom_hline(aes(yintercept=log10(seconds)),
             data=lab.df,
             color="grey")+
  geom_ribbon(aes(log10(max.data), ymin=log10(min), ymax=log10(max),
                  fill=algo.tex),
              data=window.seconds,
              alpha=0.5)+
  geom_line(aes(log10(max.data), log10(median),
                  color=algo.tex),
            data=window.seconds)+
  geom_text(aes(2.5, log10(seconds), label=label),
            data=lab.df,
            size=3,
            color="grey",
            vjust=-0.5)+
  ylab("log10(seconds)")
my.method <- list("last.points", dl.trans(x=x+0.1))
dl.log <- direct.label(gg.log, "last.polygons")+
  scale_x_continuous(
    "log10(data points to segment)")+
  coord_cartesian(xlim=c(2, 6.7), expand=FALSE)
tikz("figure-PDPA-timings-small.tex", 3.3, 2.5)
print(dl.log)
dev.off()

gg.log <- ggplot()+
  theme_bw()+
  geom_hline(aes(yintercept=seconds),
             data=lab.df,
             color="grey")+
  geom_ribbon(aes(max.data, ymin=min, ymax=max,
                  fill=algo.tex),
              data=window.seconds,
              alpha=0.5)+
  geom_line(aes(max.data, median,
                color=algo.tex),
            data=window.seconds)+
  geom_text(aes(10^2.5, seconds, label=label),
            data=lab.df,
            size=3,
            color="grey",
            vjust=-0.5)+
  scale_y_log10("seconds (log scale)")+
  scale_x_log10(
    "data points to segment (log scale)",
    breaks=c(1e3, 1e4, range(all.timings$data)))
my.method <- list("last.points", dl.trans(x=x+0.1))
dl.log <- direct.label(gg.log, "last.polygons")+
  coord_cartesian(xlim=c(min(all.timings$data), 3e6))
tikz("figure-PDPA-timings-log-log.tex", 3.3, 2.5)
print(dl.log)
dev.off()
pdf("figure-PDPA-timings-log-log.pdf", 3.3, 2.5)
print(dl.log)
dev.off()

gg.log <- ggplot()+
  geom_hline(aes(yintercept=log10(seconds)),
             data=lab.df,
             color="grey")+
  geom_text(aes(2.5, log10(seconds), label=label),
            data=lab.df,
            vjust=-0.5)+
  ggtitle("Timings on 2752 histone mark ChIP-seq data sets")+
  geom_point(aes(log10(data), log10(seconds), color=algorithm),
             shape=1,
             data=all.timings)
print(gg.log)

my.method <- list("last.points", dl.trans(x=x+0.1))
dl.log <- direct.label(gg.log, "last.polygons")+
  scale_x_continuous(
    "log10(data points to segment)",
    limits=c(min(log10(all.timings$data)), 6.2))

pdf("figure-PDPA-timings.pdf", 8, 5)
print(dl.log)
dev.off()
