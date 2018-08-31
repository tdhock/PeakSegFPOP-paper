library(data.table)
library(ggplot2)
library(directlabels)

target.intervals.models <- fread("target.intervals.models.csv")
labeled_problems_features <- fread("labeled_problems_features.csv")
select.dt <- labeled_problems_features[, data.table(prob.dir)]
bench.models <- target.intervals.models[select.dt, on=list(prob.dir)][log(bedGraph.lines) < penalty & penalty < bedGraph.lines & 1000 < bedGraph.lines]
bench.models[, gigabytes := megabytes/1024]

jss.evaluations <- readRDS("jss.evaluations.rds")[others.penalty!=Inf]

others.tall <- melt(
  jss.evaluations,
  measure.vars=c("others.seconds", "others.megabytes"))
others.tall[, var := sub("others.", "", variable)]

prob.stats <- others.tall[, list(
  OP=.N,
  sum=sum(value),
  median=median(value),
  q95=quantile(value, 0.95),
  q05=quantile(value, 0.05)
  ), by=list(var, bedGraph.lines, segments, peaks)]

target.stats <- others.tall[, list(
  sum=sum(value),
  median=median(value),
  q95=quantile(value, 0.95),
  q05=quantile(value, 0.05)
  ), by=list(var, target.N)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(var ~ ., scales="free")+
  geom_ribbon(aes(
    target.N, ymin=q05, ymax=q95),
    alpha=0.5,
    data=target.stats)+
  geom_line(aes(
    target.N, median),
    size=1,
    data=target.stats)+
  geom_point(aes(
    bedGraph.lines, sum),
    shape=1,
    data=prob.stats)+
  scale_x_log10("N = number of data to segment (log scale)")+
  scale_y_log10("(log scales)")

evals.dt <- prob.stats[var=="seconds"]
evals.dt[, SN := segments-1]
evals.tall <- melt(
  evals.dt,
  measure.vars=c("SN", "OP", "peaks"),
  variable.name="algo",
  value.name="evaluations")
algo.key <- c(
  peaks="O(sqrt N) peaks\nin zero-error model",
  SN="Segment Neighborhood\nO(sqrt N) DP iterations",
  OP="Optimal Partitioning\nO(log N) DP iterations")
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
  geom_point(aes(
    bedGraph.lines, evaluations, color=algorithm),
    shape=1,
    data=evals.tall)+
  scale_x_log10(
    "N = number of data to segment (log scale)",
    limits=c(NA, 1e8)
  )+
  scale_y_log10("Number of O(N log N) dynamic programming iterations")
dl <- direct.label(leg, list("last.qp", dl.trans(x=x+0.1)))
print(dl)

dl.ref <- dl+
  geom_line(aes(
    N.data, value, group=paste(ref.name, fun.name)),
    data=ref.tall)+
  geom_text(aes(
    N.data, value, label=fun.name),
    hjust=0,
    data=ref.tall[N.data==max(N.data)])
pdf("jss-figure-evaluations-ref.pdf", 3, 3)
print(dl.ref)
dev.off()

pdf("jss-figure-evaluations.pdf", 3, 3)
print(dl)
dev.off()
