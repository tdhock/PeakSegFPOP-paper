source("jss-packages.R")

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)][log(bedGraph.lines) < penalty & penalty < bedGraph.lines & 1000 < bedGraph.lines]
bench.models[, gigabytes := megabytes/1024]

jss.variable.peaks <- readRDS("jss.variable.peaks.rds")[others.penalty!=Inf]
jss.variable.peaks[, others.minutes := others.seconds/60]
jss.variable.peaks[, others.gigabytes := others.megabytes/1024]
others.tall <- melt(
  jss.variable.peaks,
  measure.vars=c("others.minutes", "others.gigabytes"))
others.tall[, var := sub("others.", "", variable)]

prob.stats <- others.tall[, list(
  OP=.N,
  sum=sum(value),
  median=median(value),
  max=max(value),
  q95=quantile(value, 0.95),
  q05=quantile(value, 0.05)
), by=list(
  var, target.N=as.integer(target.N),
  bedGraph.lines, segments=loss.peaks*2+1, peaks=loss.peaks)]

algo.key <- c(
  peaks="O(sqrt N) peaks\nin zero-error model",
  SN="Segment\nNeighborhood",
  OP="Optimal\nPartitioning")
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
sn.color <- abbrev.colors[["SN"]]
algo.colors <- structure(abbrev.colors, names=algo.key[names(abbrev.colors)])
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
evals.tall[, N.fac := paste("N =", format(bedGraph.lines, big.mark=","))]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ N.fac)+
  geom_point(aes(
    peaks, evaluations, color=algorithm),
    data=evals.tall[peaks <= 10])+
  coord_equal()+
  scale_x_continuous(breaks=seq(0, 10, by=2))+
  scale_y_continuous(breaks=seq(0, 20, by=2))
    
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ N.fac)+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(
    peaks, evaluations, color=algorithm),
    data=evals.tall)
    
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ target.N)+
  scale_x_log10()+
  scale_y_log10()+
  geom_line(aes(
    peaks, evaluations, group=bedGraph.lines),
    color=op.color,
    data=evals.tall[algo=="OP"])+
  geom_abline(
    slope=2, intercept=0,
    color=sn.color)
    
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ target.N)+
  geom_point(aes(
    peaks, evaluations, color=algorithm),
    data=evals.tall[algo=="OP" & peaks <= 10])+
  geom_abline(aes(
    slope=slope, intercept=intercept,
    color=algorithm),
    data=data.table(slope=2, intercept=0, algorithm="Segment\nNeighborhood"))+
  coord_equal()+
  scale_x_continuous(breaks=seq(0, 10, by=2))+
  scale_y_continuous("Number of DP iterations (log scale)", breaks=seq(0, 20, by=2))+
  scale_color_manual(values=algo.colors)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ N.fac)+
  scale_size_manual(values=c(
    "Optimal\nPartitioning"=1.5,
    "Segment\nNeighborhood"=1))+
  geom_line(aes(
    peaks, evaluations, color=algorithm, size=algorithm),
    data=evals.tall[peaks <= 10])+
  coord_equal()+
  scale_x_continuous(breaks=seq(0, 10, by=2))+
  scale_y_continuous("Number of DP iterations (log scale)", breaks=seq(0, 20, by=2))+
  scale_color_manual(values=algo.colors)

op.evals <- evals.tall[algo=="OP"]
first.peaks <- 3
N.peaks <- 10^seq(log10(first.peaks), 5, l=30)
fun.list <- list(
  ##N=identity,
  "$\\log(N)$"=log
  ##"loglog(N)"=function(x)log(log(x)),
  ##"sqrt(N)"=sqrt
)
first.dt <- op.evals[peaks==first.peaks, list(
  evaluations=mean(evaluations)
), by=list(target.N)]
ref.tall.list <- list()
for(first.i in 1:nrow(first.dt)){
  first <- first.dt[first.i]
  for(fun.name in names(fun.list)){
    fun <- fun.list[[fun.name]]
    first.y <- fun(first.peaks)
    ref.tall.list[[paste(fun.name, first.i)]] <- data.table(
      N.peaks,
      first,
      fun.name,
      value=fun(N.peaks)/first.y*first$evaluations)
  }
}
(ref.tall <- do.call(rbind, ref.tall.list))

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ target.N, labeller=function(df){
    df$target.N <- paste("N $\\approx$", df$target.N)
    df
  })+
  scale_x_log10("Number of peaks $P$ (log scale)")+
  scale_y_log10("Number of DP iterations (log scale)")+
  geom_point(aes(
    peaks, evaluations),
    color=op.color,
    data=op.evals)+
  geom_abline(
    slope=2, intercept=0,
    size=1,
    color=sn.color)

gg+
  geom_line(aes(
    N.peaks, value),
    data=ref.tall)


gg.zoom <- ggplot()+
  ggtitle("Zoom to $P \\leq 10$ peaks\n(linear scales)")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_x_continuous(
    "Number of peaks $P$",
    breaks=seq(0, 10, by=2))+
  scale_y_continuous("Number of DP iterations",
breaks=seq(0, 20, by=2))+
  geom_point(aes(
    peaks, evaluations),
    color=op.color,
    data=op.evals[target.N==max(target.N) & peaks <= 10])+
  geom_abline(
    slope=2, intercept=0,
    size=1,
    color=sn.color)+
  coord_equal()+
  scale_color_manual(values=abbrev.colors, guide=FALSE)+
  geom_text(aes(
    x, y, label=label, color=algo, hjust=hjust),
    size=3,
    data=rbind(
      data.table(x=4, y=14, algo="SN", label="SN\nfaster\nfor\n$P<5$", hjust=1),
      data.table(x=8, y=14, algo="OP", label="OP\nfaster\nfor\n$P>5$", hjust=0)),
    vjust=0.5)
tikz("jss-figure-variable-peaks-zoom.tex", 3, 3)
print(gg.zoom)
dev.off()

## only show figure for largest data set size.
gg <- ggplot()+
  ggtitle("All timings (log scales)")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_x_log10("Number of peaks $P$")+
  scale_y_log10("Number of DP iterations")+
  geom_point(aes(
    peaks, evaluations),
    color=op.color,
    data=op.evals[target.N==max(target.N)])+
  geom_abline(
    slope=2, intercept=0,
    size=1,
    color=sn.color)+
  scale_color_manual(values=abbrev.colors, guide=FALSE)+
  geom_text(aes(
    x, y, label=label, color=algo),
    size=3,
    data=rbind(
      data.table(x=10, y=25, algo="SN", label="Segment\nNeighborhood\nGPDPA\n$O(P)$ iterations"),
      data.table(x=100, y=11, algo="OP", label="Optimal\nPartitioning\nGFPOP\n$O(\\log P)$ iterations")),
    vjust=1,
    hjust=0)
gg

tikz("jss-figure-variable-peaks.tex", 3, 3)
print(gg)
dev.off()
