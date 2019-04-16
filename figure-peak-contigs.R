library(data.table)
peak.contigs <- fread("peak.contigs.csv")


peak.contigs[, bin.bases := 10^round(log10(bases))]
bin.stats <- peak.contigs[, list(
  median=as.double(median(peaks)),
  q25=quantile(peaks, 0.25),
  q75=quantile(peaks, 0.75)
), by=list(model, experiment, bin.bases)]

x.vec <- unique(bin.stats$bin.bases)

## ???
## lm.dt <- data.table(fun.name=c("identity", "sqrt"))[, {
##   fun <- get(fun.name)
##   peak.contigs[, {
##     fit <- lm(fun(peaks) ~ bases)
##     data.table(bases=x.vec, peaks=predict(fit, data.table(bases=x.vec)))
##   }, by=list(experiment, model)]
## }, by=list(fun.name)]
    

xref <- 1e5
bin.stats[, min.peaks := median[bin.bases==xref], by=list(model, experiment)]
bin.stats[, sqrt := sqrt(bin.bases)/sqrt(xref)*min.peaks]
bin.stats[, linear := bin.bases/xref*min.peaks]
ref.dt <- melt(bin.stats, measure.vars=c("sqrt", "linear"))
library(ggplot2)
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ model)+
  geom_point(aes(
    bases, peaks),
    shape=21,
    data=peak.contigs)+
  geom_line(aes(
    bin.bases, value, color=variable),
    data=ref.dt)+
  scale_x_log10(limits=c(1e4, 1e9))+
  scale_y_log10()
directlabels::direct.label(gg, "last.polygons")

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ model)+
  geom_ribbon(aes(
    bin.bases, ymin=q25, ymax=q75),
    alpha=0.5,
    data=bin.stats)+
  geom_line(aes(
    bin.bases, median),
    data=bin.stats)+
  geom_line(aes(
    bin.bases, value, color=variable),
    data=ref.dt)+
  scale_x_log10(
    "Bases per contig in the human genome hg19",
    limits=c(1e4, 1e9))+
  scale_y_log10(
    "Peaks predicted per contig")
dl <- directlabels::direct.label(gg, "last.polygons")

png("figure-peak-contigs.png", 600, 300, res=100)
print(dl)
dev.off()
