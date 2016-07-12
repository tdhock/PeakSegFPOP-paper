source("packages.R")

load("../PeakSeg-paper/dp.timings.RData")
load("PDPA.timings.RData")
(objs <- load("dp.timings.reverse.RData"))
stopifnot(nrow(PDPA.timings)==2752)
stopifnot(nrow(dp.timings)==2752)
stopifnot(nrow(dp.timings.reverse)==2752)

algo <- function(algorithm, ...){
  data.frame(algorithm, ...)
}
all.timings <- rbind(
  algo("cDPA\nO(N^2)", dp.timings),
  ## algo("DP.fwd", dp.timings),
  ## algo("DP.rev", dp.timings.reverse),
  algo("cPDPA\nO(N log N)", PDPA.timings))
  
gg.linear <- ggplot()+
  ylab("hours")+
  geom_point(aes(data, seconds/60/60, color=algorithm),
             shape=1,
             data=all.timings)
print(gg.linear)

lab.df <- data.frame(
  seconds=c(1, 60, 60*60),
  label=c("1 second", "1 minute", "1 hour"))
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

pdf("figure-PDPA-timings.pdf", 5, 5)
print(dl.log)
dev.off()
