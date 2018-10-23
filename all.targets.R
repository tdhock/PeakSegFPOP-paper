source("packages.R")

load("all.modelSelection.RData")
all.modelSelection[, errors := total.errors]

one <- all.modelSelection[chunk.name=="H3K36me3_AM_immune/10" & sample.id=="McGill0322"]
targetIntervals(one, "algo")

all.min <- all.modelSelection[, list(
  min.errors=min(total.errors)
  ), by=list(chunk.name, sample.id, algo)]

wide.min <- dcast(all.min, chunk.name+sample.id~algo)

## PDPA never has fewer min errors than GPDPA.
wide.min[GPDPA<PDPA]
wide.min[, table(GPDPA,PDPA)]

## CDPA has fewer errors for 38 problems, GPDPA has fewer for 31
## problems...
wide.min[, table(GPDPA,CDPA)]
wide.min[GPDPA<CDPA]
wide.min[CDPA<GPDPA]

all.targets <- targetIntervals(all.modelSelection, c("chunk.name", "sample.id", "algo"))

save(all.targets, file="all.targets.RData")
