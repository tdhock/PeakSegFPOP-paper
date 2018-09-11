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
  scale_y_continuous("Number of $O(N \\log N)$
DP iterations (log scale)", breaks=seq(0, 20, by=2))+
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
  scale_y_continuous("Number of $O(N \\log N)$
DP iterations (log scale)", breaks=seq(0, 20, by=2))+
  scale_color_manual(values=algo.colors)

op.evals <- evals.tall[algo=="OP"]
bigO <- function(dt, x.var, y.var, fun.list, group.vec){
  dt[, {
    x <- .SD[[x.var]]
    y <- .SD[[y.var]]
    x.seq <- exp(seq(log(10), log(1e5), l=10))
    data.table(fun=names(fun.list))[, {
      f <- fun.list[[fun]]
      first.y <- f(min(x.seq))
      dt <- data.table(x=x.seq, y=f(x.seq)/first.y*10)
      dt[[x.var]] <- dt$x
      dt[[y.var]] <- dt$y
      dt
    }, by=list(fun)]
  }, by=group.vec]
}
(o <- bigO(
      op.evals,
      "peaks",
      "evaluations",
      list(
        "$O(\\log P)$"=log
        ),
      "target.N"))

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ target.N, labeller=function(df){
    df$target.N <- paste("N $\\approx$", df$target.N)
    df
  })+
  scale_x_log10("Number of peaks (log scale)")+
  scale_y_log10("Number of $O(N \\log N)$
DP iterations (log scale)")+
  geom_point(aes(
    peaks, evaluations),
    color=op.color,
    data=op.evals)+
  geom_abline(
    slope=2, intercept=0,
    color=sn.color)+
  geom_line(aes(
    peaks, evaluations, group=fun),
    data=o)
print(gg)
    
tikz("jss-figure-variable-peaks.tex", 6, 3)
print(gg)
dev.off()
