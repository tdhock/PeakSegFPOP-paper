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
all.timings <- rbind(
  algo("PeakSegDP\nO(N^2)", dp.timings),
  algo("Segmentor\nO(N log N)", Segmentor.timings),
  ## algo("DP.fwd", dp.timings),
  ## algo("DP.rev", dp.timings.reverse),
  algo("coseg\nO(N log N)", PDPA.timings))
algo.code <- c(
  PeakSegDP="PeakSegDP\n$O(N^2)$",
  Segmentor="Segmentor\n$O(N \\log N)$",
  coseg="coseg\n$O(N \\log N)$")
all.timings[, algo.tex := algo.code[sub("\n.*", "", algorithm)]]

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
  xlab("log10(seconds)")+
  ## geom_point(aes(log10(data), log10(seconds), color=algo.tex),
  ##            shape=1,
  ##            data=all.timings)+
  geom_ribbon(aes(log10(max.data), ymin=log10(min), ymax=log10(max),
                  fill=algo.tex),
              data=window.seconds,
              alpha=0.5)+
  geom_line(aes(log10(max.data), log10(median),
                  color=algo.tex),
            data=window.seconds)+
  geom_text(aes(2.5, log10(seconds), label=label),
            data=lab.df,
            vjust=-0.5)+
  ylab("log10(seconds)")
my.method <- list("last.points", dl.trans(x=x+0.1))
dl.log <- direct.label(gg.log, "last.polygons")+
  scale_x_continuous(
    "log10(data points to segment)",
    limits=c(min(log10(all.timings$data)), 6.2))
tikz("figure-PDPA-timings-small.tex", 3.3, 2.5)
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
