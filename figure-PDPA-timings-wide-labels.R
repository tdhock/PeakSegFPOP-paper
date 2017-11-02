source("packages.R")

options(
  tikzDocumentDeclaration=paste(
    "\\documentclass{beamer}",
    "\\usepackage{amsmath,amssymb,amsthm}"),
  tikzMetricsDictionary="tikzMetricsBeamer")

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
CDPA.name <- "CDPA = Constrained Dynamic Programming Algo
Up-down constrained, approximate solution, $O(d^2)$
H {\\it et al.} 2015, R pkg PeakSegDP"
PDPA.name <- "PDPA = Pruned Dynamic Programming Algorithm
Unconstrained, optimal solution, $O(d \\log d)$
Cleynen {\\it et al.} 2014, R pkg Segmentor3IsBack"
GPDPA.name <- "GPDPA = Generalized PDPA
Up-down constrained, optimal solution, $O(d \\log d)$
Proposed, R pkg PeakSegOptimal"
all.timings <- rbind(
  algo(CDPA.name, dp.timings),
  algo(PDPA.name, Segmentor.timings),
  ## algo("DP.fwd", dp.timings),
  ## algo("DP.rev", dp.timings.reverse),
  algo(GPDPA.name, PDPA.timings))
totex <- function(a){
  a
}
all.timings[, algo.tex := totex(algorithm)]

totals <- all.timings[, list(seconds=sum(seconds)), by=algorithm]
totals[, minutes := seconds/60]
totals[, hours := minutes/60]
  
gg.linear <- ggplot()+
  ylab("hours")+
  geom_point(aes(data, seconds/60/60, color=algorithm),
             shape=1,
             data=all.timings)

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

algo.colors <- c("#66C2A5", "#FC8D62",
                 "grey30")
#"#8DA0CB")
names(algo.colors) <- totex(c(PDPA.name, CDPA.name, GPDPA.name))
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
  scale_y_log10("seconds")+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=algo.colors)+
  guides(color="none", fill="none")+
  scale_x_log10(
    "$d$ = data points to segment",
    breaks=c(1e3, 1e4, range(all.timings$data)))
my.method <- list("last.points", dl.trans(x=x+0.1))

space <- 0.1
data.col <- "left"
na.col <- "right"
print.h <- function(d, ...){
  print(d$h)
  d
}

my.polygons <- list("last.points", "calc.boxes",
       function(d,...){
         for(xy in c("x", "y")){
           d[[sprintf("%s.%s", data.col, xy)]] <- d[[xy]]
           d[[sprintf("%s.%s", na.col, xy)]] <- NA
         }
         d$x <- d$x + space
         d
       },
                    print.h,
                    dl.trans(h=h*1.2),
                    print.h,
                    "calc.borders",
                    print.h,
                    qp.labels("y","bottom","top", make.tiebreaker("x","y"), ylimits),
                    print.h,
                    "calc.borders",
                    dl.trans(right=right+space),
                    "draw.polygons", print.h)
dl.log <- direct.label(gg.log, list(cex=0.6, "my.polygons"))+
  coord_cartesian(xlim=c(min(all.timings$data), 1e10))
print(dl.log)
tikz("figure-PDPA-timings-wide-labels.tex", 5, 2.2)
print(dl.log)
dev.off()
##system("pdflatex figure-PDPA-timings-wide-labels")
