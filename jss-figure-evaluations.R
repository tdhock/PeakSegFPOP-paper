source("jss-packages.R")

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)][log(bedGraph.lines) < penalty & penalty < bedGraph.lines & 1000 < bedGraph.lines]
bench.models[, gigabytes := megabytes/1024]

inf.evaluations <- readRDS("jss.evaluations.rds")
jss.evaluations <- inf.evaluations[others.penalty!=Inf]
jss.evaluations[peaks != loss.peaks, list(bedGraph.lines, peaks, loss.peaks)]

not.found <- inf.evaluations[bedGraph.lines==66031]
inf.pen <- 25000
not.found[, x := ifelse(others.penalty==Inf, inf.pen, others.penalty)]
not.found[, y := ifelse(others.penalty==Inf, others.total.loss-inf.pen*peaks, others.total.loss+(others.peaks-peaks)*others.penalty)]
biggest.it <- 4
it.text <- not.found[others.iteration <= biggest.it]
it.points <- not.found[biggest.it < others.iteration]
gg <- ggplot()+
  ggtitle(sprintf(
    "Numbers=iterations 1-%d\nPoints=iterations %d-%d",
    biggest.it,
    min(it.points$others.iteration),
    max(it.points$others.iteration)))+
  theme_bw()+
  xlab("Penalty $\\lambda$")+ 
  ylab("$G(\\lambda)=F(\\lambda)-\\lambda P^*$")+
  geom_abline(aes(
    slope=others.peaks-peaks,
    intercept=others.total.loss),
    color="grey50",
    data=not.found)+
  geom_text(aes(
    ifelse(others.penalty==Inf, inf.pen, others.penalty),
    ifelse(
      others.penalty==Inf,
      others.total.loss-inf.pen*peaks,
      others.total.loss+(others.peaks-peaks)*others.penalty),
    label=others.iteration),
    size=3,
    data=it.text)+
  geom_point(aes(
    others.penalty,
    others.total.loss+others.peaks*others.penalty-peaks*others.penalty),
    shape=1,
    data=it.points)
print(gg)

o <- 0.05
it.text[, vjust := c(0, 1, 0.5, 1.2, 0.9)]
it.text[, hjust := c(0-o, 1+o, 0-o, 0-o, 0-o)]
gg <- ggplot()+
  ggtitle("All 13 iterations")+
  theme_bw()+
  coord_cartesian(
    xlim=c(-1000, inf.pen),
    ylim=c(-1600000, -10000),
    expand=FALSE)+
  scale_x_continuous(
    "Penalty $\\lambda$",
    breaks=seq(0, 20000, by=5000))+ 
  ylab("$G(\\lambda)=F(\\lambda)-\\lambda P^*$")+
  geom_abline(aes(
    slope=others.peaks-peaks,
    intercept=others.total.loss),
    color="grey50",
    data=not.found)+
  geom_l(aes(
    x, y, hjust=hjust, vjust=vjust,
    label=sprintf("it=%d, $P=%d$", others.iteration, others.peaks)),
    size=3,
    alpha=0.7,
    color="red",
    data=it.text)+
  ## geom_text(aes(
  ##   x, y, hjust=hjust, vjust=vjust,
  ##   label=sprintf("it=%d, $P=%d$", others.iteration, others.peaks)),
  ##   size=3,
  ##   color="red",
  ##   data=it.text)+
  geom_point(aes(
    x, y),
    shape=1,
    color="red",
    data=not.found)
print(gg)
 
pdf("jss-figure-evaluations-concave.pdf")
print(gg)
dev.off()
tikz("jss-figure-evaluations-concave.tex", 3, 2)
print(gg)
dev.off()

biggest.it <- 8
it.text <- not.found[others.iteration <= biggest.it]
it.points <- not.found[biggest.it < others.iteration]
o <- 0.05
it.points[, hjust := c(0-o, 1+o, 0-o,      1+o, 0-o)]
it.points[, vjust := c(1+o, 0-o, 1+o,  1+o, 0-o)]
gg.zoom <- ggplot()+
  ggtitle(sprintf(
    "Zoom to iterations %d-%d",
    min(it.points$others.iteration),
    max(it.points$others.iteration)
  ))+
  theme_bw()+
  xlab("penalty")+
  ylab("concave function to maximize")+
  geom_abline(aes(
    slope=others.peaks-peaks,
    intercept=others.total.loss),
    color="grey50",
    data=it.points)+
  geom_text(aes(
    others.penalty,
    others.total.loss+others.peaks*others.penalty-peaks*others.penalty,
    vjust=vjust,
    hjust=hjust,
    label=sprintf(
      "iteration=%d
peaks=%d", others.iteration, others.peaks)),
    data=it.points)+
  geom_point(aes(
    others.penalty,
    others.total.loss+others.peaks*others.penalty-peaks*others.penalty),
    shape=1,
    data=it.points)
print(gg.zoom)

pdf("jss-figure-evaluations-concave-zoom.pdf")
print(gg.zoom)
dev.off()

gg.zoom <- ggplot()+
  ggtitle(sprintf(
    "Zoom to iterations %d-%d",
    min(it.points$others.iteration),
    max(it.points$others.iteration)
  ))+
  xlab("Penalty $\\lambda$")+ 
  scale_y_continuous(
    "$G(\\lambda)=F(\\lambda)-\\lambda P^*$",
    limits=c(NA, -71830)
  )+
  theme_bw()+
  geom_abline(aes(
    slope=others.peaks-peaks,
    intercept=others.total.loss),
    color="grey50",
    data=it.points)+
  geom_label(aes(
    others.penalty,
    others.total.loss+others.peaks*others.penalty-peaks*others.penalty,
    vjust=vjust,
    hjust=hjust,
    label=sprintf(
      "it=%d
$P=%d$", others.iteration, others.peaks)),
size=3,
alpha=0.5,
color="white",
    data=it.points)+
  geom_text(aes(
    others.penalty,
    others.total.loss+others.peaks*others.penalty-peaks*others.penalty,
    vjust=vjust,
    hjust=hjust,
    label=sprintf(
      "it=%d
$P=%d$", others.iteration, others.peaks)),
size=3,
    data=it.points)+
  geom_point(aes(
    others.penalty,
    others.total.loss+others.peaks*others.penalty-peaks*others.penalty),
    shape=1,
    data=it.points)

gg.zoom <- ggplot()+
  ggtitle(sprintf(
    "Zoom to iterations %d--%d",
    min(it.points$others.iteration),
    max(it.points$others.iteration)
  ))+
  xlab("Penalty $\\lambda$")+ 
  scale_y_continuous(
    "$G(\\lambda)=F(\\lambda)-\\lambda P^*$",
    limits=c(NA, -71830)
  )+
  theme_bw()+
  geom_abline(aes(
    slope=others.peaks-peaks,
    intercept=others.total.loss),
    color="grey50",
    data=it.points)+
  geom_l(aes(
    others.penalty,
    others.total.loss+others.peaks*others.penalty-peaks*others.penalty,
    vjust=vjust,
    hjust=hjust,
    label=sprintf(
      "it=%d, $P=%d$", others.iteration, others.peaks)),
    size=3,
    alpha=0.7,
    color="red",
    data=it.points)+
  ## geom_text(aes(
  ##   others.penalty,
  ##   others.total.loss+others.peaks*others.penalty-peaks*others.penalty,
  ##   vjust=vjust,
  ##   hjust=hjust,
  ##   label=sprintf(
  ##     "it=%d, $P=%d$", others.iteration, others.peaks)),
  ##   size=3,
  ##   color="red",
  ##   data=it.points)+
  geom_point(aes(
    others.penalty,
    others.total.loss+others.peaks*others.penalty-peaks*others.penalty),
    shape=1,
    color="red",
    data=it.points)

tikz("jss-figure-evaluations-concave-zoom.tex", 3, 2)
print(gg.zoom)
dev.off()

jss.evaluations[, others.minutes := others.seconds/60]
jss.evaluations[, others.gigabytes := others.megabytes/1024]
others.tall <- melt(
  jss.evaluations,
  measure.vars=c("others.minutes", "others.gigabytes"))
others.tall[, var := sub("others.", "", variable)]

prob.stats <- others.tall[, list(
  OP=.N,
  sum=sum(value),
  median=median(value),
  max=max(value),
  q95=quantile(value, 0.95),
  q05=quantile(value, 0.05)
  ), by=list(var, bedGraph.lines, segments=peaks*2+1, peaks)]

target.stats <- others.tall[, list(
  sum=sum(value),
  median=median(value),
  q95=quantile(value, 0.95),
  q05=quantile(value, 0.05)
  ), by=list(var, target.N)]

both.points <- rbind(
  prob.stats[var=="minutes", data.table(
    bedGraph.lines, value=sum, var, algorithm="find model")],
  prob.stats[var=="gigabytes", data.table(
    bedGraph.lines, value=max, var, algorithm="find model")],
  others.tall[, data.table(bedGraph.lines, value, var, algorithm="solve one")])
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(var ~ ., scales="free")+
  geom_point(aes(
    bedGraph.lines, value, color=algorithm),
    shape=1,
    data=both.points)+
  scale_x_log10("N = number of data to segment (log scale)")+
  scale_y_log10("(log scales)")

algo.key <- c(
  peaks="O(sqrt N) peaks\nin zero-error model",
  SN="Segment\nNeighborhood\nO(sqrt N)\niterations",
  OP="Optimal\nPartitioning\nO(log N)\niterations")
abbrev.colors <- c(
  "#E41A1C",#red
  "#377EB8",#blue
  OP="#4DAF4A",#green
  "#984EA3",#purple
  "#FF7F00",#orange
  "#FFFF33", #yellow
  "#A65628",#brown
  "#F781BF",#pink
  SN="#999999")#grey
op.color <- abbrev.colors[["OP"]]
algo.colors <- structure(abbrev.colors, names=algo.key[names(abbrev.colors)])
text.dt <- both.points[bedGraph.lines==max(bedGraph.lines), list(
  y=(max(value)+min(value))/2,
  hjust=0,
  vjust=1
), by=list(var, algorithm, bedGraph.lines)]
d <- function(bedGraph.lines, y, algorithm, var, hjust, vjust){
  data.table(bedGraph.lines, y, algorithm, var, hjust, vjust)
}
text.dt <- rbind(
  d(3e3, 3e3, "find model", "gigabytes", 0, 1),
  d(1e7, 0.002, "solve one", "gigabytes", 1, 0),
  d(3e3, 500, "find model", "minutes", 0, 1),
  d(1e7, 0.01, "solve one", "minutes", 1, 0))
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(var ~ ., scales="free")+
  geom_text(aes(
    bedGraph.lines, y, hjust=hjust, vjust=vjust, label={
      l <- ifelse(
        algorithm=="solve one",
        "Solve for one penalty\nO(N log N) time",
        "Find zero-error model\nwith O(sqrt N) peaks\nO(N(logN)^2) time")
      ifelse(var=="gigabytes", sub("time", "space", l), l)
    }),
    size=3,
    color=op.color,
    data=text.dt)+
  geom_ribbon(aes(
    target.N, ymin=q05, ymax=q95),
    alpha=0.5,
    fill=op.color,
    data=target.stats)+
  geom_line(aes(
    target.N, median),
    size=1,
    color=op.color,
    data=target.stats)+
  geom_point(aes(
    bedGraph.lines, value),
    shape=1,
    color=op.color,
    data=both.points[algorithm=="find model"])+
  scale_x_log10(
    "N = data to segment (log scale)"
  )+
  scale_y_log10(
    "Computational resources to find
zero-error model with O(sqrt N) peaks
via Optimal Partitioning (log scales)",
    labels=paste)
pdf("jss-figure-evaluations-computation.pdf", 3.5, 3)
print(gg)
dev.off() 

evals.dt <- prob.stats[var=="minutes"]
evals.dt[, SN := segments-1]
evals.tall <- melt(
  evals.dt,
  measure.vars=c(
    "SN",
    ##"peaks",
    "OP"
  ),
  variable.name="algo",
  value.name="evaluations")
evals.tall[, algorithm := algo.key[paste(algo)] ]
N.data <- 10^seq(4, 7, l=100)
fun.list <- list(
  N=identity,
  "log(N)"=log,
  "loglog(N)"=function(x)log(log(x)),
  "sqrt(N)"=sqrt)
ref.line.list <- list(
  OP=list(y=9, lines=c("log(N)", "sqrt(N)", "loglog(N)")),
  SN=list(y=60, lines=c("N", "log(N)", "sqrt(N)")))
ref.tall.list <- list()
for(ref.name in names(ref.line.list)){
  ref.info <- ref.line.list[[ref.name]]
  for(fun.name in ref.info$lines){
    fun <- fun.list[[fun.name]]
    first.y <- fun(min(N.data))
    ref.tall.list[[paste(fun.name, ref.name)]] <- data.table(
      N.data,
      ref.name,
      fun.name,
      value=fun(N.data)/first.y*ref.info$y)
  }
}
ref.tall <- do.call(rbind, ref.tall.list)
leg <- ggplot()+
  theme_bw()+
  scale_color_manual(values=algo.colors)+
  geom_point(aes(
    bedGraph.lines, evaluations, color=algorithm),
    data=evals.tall)+
  scale_x_log10(
    "N = data to segment (log scale)",
    limits=c(NA, 10^8.5),
    breaks=10^seq(4, 7, by=1)
  )+
  scale_y_log10(
    "Number of O(N log N) dynamic
programming iterations to find
model with O(sqrt N) peaks (log scale)",
    limits=c(NA, 2000)
  )
m <- list(
  cex=0.75,
  dl.trans(x=x+0.1),
  "last.points", "calc.boxes",
  "reduce.cex.lr",
  function(d,...){
    d$h <- d$h * 1.5
    d
  },
  "calc.borders",
  qp.labels("y","bottom","top", make.tiebreaker("x","y"), ylimits),
  "calc.borders")
dl <- direct.label(leg, m)
print(dl)

dl.ref <- dl+
  geom_line(aes(
    N.data, value, group=paste(ref.name, fun.name)),
    data=ref.tall)+
  geom_text(aes(
    N.data, value, label=fun.name),
    hjust=0,
    data=ref.tall[N.data==max(N.data)])
pdf("jss-figure-evaluations-ref.pdf", 3.5, 3)
print(dl.ref)
dev.off()

pdf("jss-figure-evaluations.pdf", 3.5, 3)
print(dl)
dev.off()

## BELOW TIKZ FIGURES

algo.key <- c(
  peaks="$O(\\sqrt N)$ peaks\nin zero-error model",
  SN="Segment\nNeighborhood\n$O(\\sqrt N)$\niterations",
  OP="Optimal\nPartitioning\n$O(\\log N)$\niterations")
abbrev.colors <- c(
  "#E41A1C",#red
  "#377EB8",#blue
  OP="#4DAF4A",#green
  "#984EA3",#purple
  "#FF7F00",#orange
  "#FFFF33", #yellow
  "#A65628",#brown
  "#F781BF",#pink
  SN="#999999")#grey
op.color <- abbrev.colors[["OP"]]
algo.colors <- structure(abbrev.colors, names=algo.key[names(abbrev.colors)])
text.dt <- both.points[bedGraph.lines==max(bedGraph.lines), list(
  y=(max(value)+min(value))/2,
  hjust=0,
  vjust=1
), by=list(var, algorithm, bedGraph.lines)]
d <- function(bedGraph.lines, y, algorithm, var, hjust, vjust){
  data.table(bedGraph.lines, y, algorithm, var, hjust, vjust)
}
text.dt <- rbind(
  d(3e3, 1e5, "find model", "gigabytes", 0, 1),
  d(1e7, 0.0005, "solve one", "gigabytes", 1, 0),
  d(3e3, 1e4, "find model", "minutes", 0, 1),
  d(1e7, 0.005, "solve one", "minutes", 1, 0))
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(var ~ ., scales="free")+
  geom_text(aes(
    bedGraph.lines, y, hjust=hjust, vjust=vjust, label={
      l <- ifelse(
        algorithm=="solve one",
        "Solve for one penalty\n$O(N \\log N)$ time",
        "Find zero-error model\nwith $O(\\sqrt N)$ peaks\n$O(N(\\log N)^2)$ time")
      ifelse(var=="gigabytes", sub("time", "space", l), l)
    }),
    size=2.5,
    color=op.color,
    data=text.dt)+
  geom_ribbon(aes(
    target.N, ymin=q05, ymax=q95),
    alpha=0.5,
    fill=op.color,
    data=target.stats)+
  geom_line(aes(
    target.N, median),
    size=1,
    color=op.color,
    data=target.stats)+
  geom_point(aes(
    bedGraph.lines, sum),
    shape=1,
    color=op.color,
    data=prob.stats)+
  scale_x_log10(
    "$N$ = data to segment (log scale)"
  )+
  scale_y_log10(
    "Computational resources to find
zero-error model with $O(\\sqrt N)$ peaks
via Optimal Partitioning (log scales)",
    labels=paste)

algo.key <- c(
  peaks="$O(\\sqrt N)$ peaks\nin zero-error model",
  SN="Segment\nNeighborhood\nGPDPA\n$O(\\sqrt N)$\niterations",
  OP="Optimal\nPartitioning\nGFPOP\n$O(\\log N)$\niterations")
abbrev.colors <- c(
  "#E41A1C",#red
  "#377EB8",#blue
  OP="#4DAF4A",#green
  "#984EA3",#purple
  "#FF7F00",#orange
  "#FFFF33", #yellow
  "#A65628",#brown
  "#F781BF",#pink
  SN="#999999")#grey
op.color <- abbrev.colors[["OP"]]
algo.colors <- structure(abbrev.colors, names=algo.key[names(abbrev.colors)])
text.dt <- both.points[bedGraph.lines==max(bedGraph.lines), list(
  y=(max(value)+min(value))/2,
  hjust=0,
  vjust=1
), by=list(var, algorithm, bedGraph.lines)]
d <- function(bedGraph.lines, y, algorithm, var, hjust, vjust){
  data.table(bedGraph.lines, y, algorithm, var, hjust, vjust)
}
text.dt <- rbind(
  d(3e3, 1e5, "find model", "gigabytes", 0, 1),
  d(1e7, 0.0005, "solve one", "gigabytes", 1, 0),
  d(3e3, 300, "find model", "minutes", 0, 1),
  d(1e7, 0.03, "solve one", "minutes", 1, 0))
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_text(aes(
    bedGraph.lines, y, hjust=hjust, vjust=vjust, label={
      l <- ifelse(
        algorithm=="solve one",
        "Solve for one penalty\n$O(N \\log N)$ time",
        "Sequential search for model\nwith $O(\\sqrt N)$ peaks\n$O(N(\\log N)^2)$ time")
      ifelse(var=="gigabytes", sub("time", "space", l), l)
    }),
    size=3,
    color=op.color,
    data=text.dt[var=="minutes"])+
  geom_ribbon(aes(
    target.N, ymin=q05, ymax=q95),
    alpha=0.5,
    fill=op.color,
    data=target.stats[var=="minutes"])+
  geom_line(aes(
    target.N, median),
    size=1,
    color=op.color,
    data=target.stats[var=="minutes"])+
  geom_point(aes(
    bedGraph.lines, sum),
    shape=1,
    color=op.color,
    data=prob.stats[var=="minutes"])+
  scale_x_log10(
    "$N$ = data to segment
(log scale)"
  )+
  scale_y_log10(
    "Time (minutes) to compute
     model with $O(\\sqrt N)$ peaks
     via GFPOP (log scale)",
labels=paste)
##print(gg)
tikz("jss-figure-evaluations-computation.tex", 3.1, 2.6)
print(gg)
dev.off() 

evals.dt <- prob.stats[var=="minutes"]
evals.dt[, SN := segments-1]
evals.tall <- melt(
  evals.dt,
  measure.vars=c(
    "SN",
    ##"peaks",
    "OP"
  ),
  variable.name="algo",
  value.name="evaluations")
evals.tall[, algorithm := algo.key[paste(algo)] ]
## N.data <- 10^seq(4, 7, l=100)
## fun.list <- list(
##   N=identity,
##   "log(N)"=log,
##   "loglog(N)"=function(x)log(log(x)),
##   "sqrt(N)"=sqrt)
## ref.line.list <- list(
##   OP=list(y=9, lines=c("log(N)", "sqrt(N)", "loglog(N)")),
##   SN=list(y=60, lines=c("N", "log(N)", "sqrt(N)")))
## ref.tall.list <- list()
## for(ref.name in names(ref.line.list)){
##   ref.info <- ref.line.list[[ref.name]]
##   for(fun.name in ref.info$lines){
##     fun <- fun.list[[fun.name]]
##     first.y <- fun(min(N.data))
##     ref.tall.list[[paste(fun.name, ref.name)]] <- data.table(
##       N.data,
##       ref.name,
##       fun.name,
##       value=fun(N.data)/first.y*ref.info$y)
##   }
## }
## ref.tall <- do.call(rbind, ref.tall.list)
leg <- ggplot()+
  theme_bw()+
  scale_color_manual(values=algo.colors)+
  geom_point(aes(
    bedGraph.lines, evaluations, color=algorithm),
    data=evals.tall)+
  scale_x_log10(
    "$N$ = data to segment
(log scale)",
    limits=c(NA, 10^8.5),
    breaks=10^seq(4, 7, by=1)
  )+
  scale_y_log10(
    "Number of $O(N \\log N)$
DP iterations to compute
model with $O(\\sqrt N)$
peaks (log scale)",
    limits=c(NA, 2000)
  )
m <- list(
  cex=0.8,
  dl.trans(x=x+0.1),
  "last.points", "calc.boxes",
  "reduce.cex.lr",
  function(d,...){
    d$h <- d$h +0.3
    d
  },
  "calc.borders",
  qp.labels("y","bottom","top", make.tiebreaker("x","y"), ylimits),
  "calc.borders")
dl <- direct.label(leg, m)
print(dl)

## dl.ref <- dl+
##   geom_line(aes(
##     N.data, value, group=paste(ref.name, fun.name)),
##     data=ref.tall)+
##   geom_text(aes(
##     N.data, value, label=fun.name),
##     hjust=0,
##     data=ref.tall[N.data==max(N.data)])
## pdf("jss-figure-evaluations-ref.pdf", 3.5, 3)
## print(dl.ref)
## dev.off()

tikz("jss-figure-evaluations.tex", 3.1, 2.6)
print(dl)
dev.off()
