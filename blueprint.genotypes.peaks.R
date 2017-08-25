genotype.mat <- readRDS("blueprint.genotypes.rds")

if(!file.exists("blueprint_peaks_matrix.tsv")){
  download.file(
    "http://hubs.hpc.mcgill.ca/~thocking/PeakSegFPOP-labels/H3K27ac-H3K4me3_TDHAM_BP/peaks_matrix.tsv",
    "blueprint_peaks_matrix.tsv")
}

if(!file.exists("blueprint_peaks_summary.tsv")){
  download.file(
    "http://hubs.hpc.mcgill.ca/~thocking/PeakSegFPOP-labels/H3K27ac-H3K4me3_TDHAM_BP/peaks_summary.tsv",
    "blueprint_peaks_summary.tsv")
}

peaks.dt <- fread("blueprint_peaks_matrix.tsv")

summary.dt <- fread("blueprint_peaks_summary.tsv")

pattern <- paste0(
  "^Mono",
  "[0-9]+",
  "_",
  "H3K27ac",
  "/",
  "(?<person>",
  paste(colnames(genotype.mat), collapse="|"),
  ")",
  "_",
  "NCMLS")
grep(pattern, names(peaks.dt), value=TRUE, perl=TRUE)
match.mat <- str_match_named(names(peaks.dt), pattern)

has.genotype <- !is.na(match.mat)

peaks.mat <- as.matrix(peaks.dt[, has.genotype, with=FALSE])
rownames(peaks.mat) <- peaks.dt$peak
rbind(colnames(peaks.mat), match.mat[has.genotype])
colnames(peaks.mat) <- match.mat[has.genotype]

summary.dt[, genotyped.up := rowSums(peaks.mat)]
summary.dt[, genotyped.down := ncol(peaks.mat)-genotyped.up]

big.peaks.tall <- data.table(melt(peaks.mat))
setnames(big.peaks.tall, c("peak.name", "person", "has.peak"))
peaks.tall <- big.peaks.tall[has.peak==1]

zoom <- 4
summary.dt[, zoomStart := peakStart - peakBases*zoom]
summary.dt[, zoomEnd := peakEnd + peakBases*zoom]
center.peaks <- summary.dt[, data.table(
  center.peak.name=peak.name, chrom, zoomStart, zoomEnd)]
other.peaks <- summary.dt[, data.table(
  other.peak.name=peak.name, chrom, peakStart, peakEnd)]
setkey(center.peaks, chrom, zoomStart, zoomEnd)
setkey(other.peaks, chrom, peakStart, peakEnd)
over.dt <- foverlaps(center.peaks, other.peaks, nomatch=0L)[center.peak.name != other.peak.name]
over.dt[center.peak.name=="chr17:29634022-29641640"]
over.dt[other.peak.name=="chr17:29634022-29641640"]

over.dt[, list(peaks.in.zoom=.N), by=list(center.peak.name)][order(peaks.in.zoom)]
over.dt[, list(peaks.in.zoom=.N), by=list(zoom.peak.name)][order(peaks.in.zoom)]

big.and.other <- big.peaks.tall[over.dt, on=list(peak.name=center.peak.name)]
setkey(big.peaks.tall, person, peak.name)
setkey(big.and.other, person, other.peak.name)
big.and.other[, other.has.peak := big.peaks.tall[big.and.other, has.peak]]
others.up.counts <- big.and.other[, list(others.up=sum(other.has.peak)), by=list(peak.name, has.peak, person)]
people.counts <- others.up.counts[, list(
  people.up=sum(has.peak),
  people.other.up=sum(1 < others.up),
  people=.N
  ), by=list(peak.name)][order(people.other.up, people.up)]

data.frame(people.counts[people.other.up==40])

down.but.nearby.up <- others.up.counts[has.peak==0 & 0<others.up]

down.but.nearby.up[, list(
  count=.N
  ), by=list(peak.name)][count==38]


