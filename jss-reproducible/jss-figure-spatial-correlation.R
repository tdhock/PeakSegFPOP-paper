library(PeakSegDisk)
library(data.table)
library(ggplot2)

## Load data set with one row for every genomic region with a
## unique aligned read, and compute mean read size in bases.
data(ChIPreads, envir=environment())
experiments <- ChIPreads[, .(
  mean.bases=mean(chromEnd-chromStart),
  median.bases=median(chromEnd-chromStart),
  chromStart=min(chromStart)
), by=list(experiment)]

## Compute data set with two representations of these aligned
## reads: count each read at each aligned base in the read, or
## just the end/last base of the read.
end.counts <- ChIPreads[, list(
  count=.N #ignores dup reads, sum(count) would not.
), by=list(experiment, chrom, chromEnd)]
aligned.dt <- rbind(
  ChIPreads[, .(
    bases.counted="each", experiment, chrom,
    chromStart, chromEnd,
    count=1)], #ignore duplicate reads.
  end.counts[, .(
    bases.counted="end", experiment, chrom,
    chromStart=chromEnd-1L, chromEnd,
    count)])

## Compute count profile for each base in these genomic regions.
seq.dt <- aligned.dt[, {
  event.dt <- rbind(
    data.table(count, pos=chromStart+1L),
    data.table(count=-count, pos=chromEnd+1L))
  edge.vec <- event.dt[, {
    as.integer(seq(min(pos), max(pos), l=100))
  }]
  event.bins <- rbind(
    event.dt,
    data.table(count=0L, pos=edge.vec))
  total.dt <- event.bins[, .(
    count=sum(count)
  ), by=list(pos)][order(pos)]
  total.dt[, cum := cumsum(count)]
  total.dt[, bin.i := cumsum(pos %in% edge.vec)]
  ## it is somewhat confusing because total.dt pos is the first base
  ## with cum, and cum goes all the way up to but not including the
  ## pos of the next row.
  total.dt[, data.table(
    chromStart=pos[-.N]-1L,
    chromEnd=pos[-1]-1L,
    count=cum[-.N],
    bin.i=bin.i[-.N])]
}, by=list(bases.counted, experiment, chrom)]
gg.data <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(
    bases.counted ~ experiment,
    scales="free",
    labeller=label_both)+
  geom_step(aes(
    chromStart/1e3, count, color=data.type),
    data=data.table(seq.dt, data.type="exact"))+
  scale_color_manual(values=c(
    exact="black",
    bins="red",
    model="deepskyblue"
  ))+
  scale_x_continuous("Position on hg19 chrom (kb = kilo bases)")
if(interactive())print(gg.data)

## Compute mean profile in bins.
bin.dt <- seq.dt[, {
  bases <- chromEnd - chromStart
  data.table(
    binStart=min(chromStart),
    binEnd=max(chromEnd),
    mean.count=sum(count*bases)/sum(bases),
    bases=sum(bases)
  )}, by=list(bases.counted, experiment, bin.i)]
gg.bins <- gg.data+
  geom_step(aes(
    binStart/1e3, mean.count, color=data.type),
    alpha=0.75,
    size=0.7,
    data=data.table(bin.dt, data.type="bins"))+
  scale_y_log10("Aligned DNA sequence reads (log scale)")
if(interactive())print(gg.bins)

## Compute optimal segmentation model with 2 peaks.
segs.dt <- seq.dt[, {
  data.dir <- file.path("figure-spatial-correlation", bases.counted, experiment)
  dir.create(data.dir, showWarnings=FALSE, recursive=TRUE)
  coverage.bedGraph <- file.path(data.dir, "coverage.bedGraph")
  fwrite(
    .SD[, .(chrom, chromStart, chromEnd, count)],
    coverage.bedGraph,
    sep="\t",
    quote=FALSE,
    col.names=FALSE)
  fit <- PeakSegDisk::sequentialSearch_dir(data.dir, 2L, verbose=1)
  data.table(fit$segments, data.type="model")
}, by=list(bases.counted, experiment)]
changes.dt <- segs.dt[, {
  .SD[-1]
}, by=list(bases.counted, experiment, data.type)]
gg.model <- gg.bins+
  geom_segment(aes(
    chromStart/1e3, mean,
    xend=chromEnd/1e3, yend=mean,
    color=data.type),
    data=segs.dt)+
  geom_vline(aes(
    xintercept=chromEnd/1e3,
    color=data.type),
    data=changes.dt)
if(interactive())print(gg.model)

## Compute difference between peak positions of two models.
peaks.dt <- segs.dt[status=="peak"]
peaks.dt[, peak.i := rep(1:2, l=.N)]
peak.pos.tall <- melt(
  peaks.dt,
  measure.vars=c("chromStart", "chromEnd"))
peak.pos.wide <- dcast(
  peak.pos.tall,
  experiment + variable + peak.i ~ bases.counted)
peak.pos.wide[, diff.bases := abs(each-end)]

read.size.panel <- "each"
bases.max.dt <- seq.dt[, .(max.count=max(count)), by=list(bases.counted)]
read.size.y <- bases.max.dt[
  read.size.panel, max.count, on=list(bases.counted)]
read.size.y <- Inf
diff.panel <- "end"
diff.y <- bases.max.dt[
  diff.panel, max.count, on=list(bases.counted)]
diff.y <- Inf
diff.vjust <- 1.1
text.size <- 3
gg.text <- gg.model+
  geom_text(aes(
    chromStart/1e3, read.size.y, label=sprintf(
      "Median read size:\n%.0f bases",
      median.bases)),
    hjust=0,
    size=text.size,
    vjust=1.1,
    data=data.table(experiments, bases.counted=read.size.panel))+
  geom_text(aes(
    end/1e3, diff.y,
    label=diff.bases,
    color=data.type),
    size=text.size,
    data=data.table(
      bases.counted=diff.panel,
      data.type="model",
      peak.pos.wide),
    vjust=diff.vjust,
    hjust=0)+
  geom_text(aes(
    chromStart/1e3, diff.y,
    label="Difference=",
    color=data.type,
  ),
  size=text.size,
  hjust=0,
  vjust=diff.vjust,
  data=data.table(
    data.type="model",
    bases.counted=diff.panel,
    experiments["H3K36me3", on=list(experiment)]))
png("jss-figure-spatial-correlation.png", 8, 3.2, units="in", res=300)
print(gg.text)
dev.off()
