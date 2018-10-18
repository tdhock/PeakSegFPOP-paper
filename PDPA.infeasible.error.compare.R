source("packages.R")

load("PDPA.infeasible.error.RData")
load("PDPA.peaks.error.RData")
load("dp.peaks.error.RData")

CDPA.error.list <- list()
for(chunk.name in names(dp.peaks.error)){
  chunk.error <- dp.peaks.error[[chunk.name]]
  for(cell.type in names(chunk.error)){
    type.dt <- data.table(
      chunk.name,
      chunk.error[[cell.type]])
    names(type.dt)[3] <- "peaks"
    CDPA.error.list[[paste(chunk.name, cell.type)]] <- type.dt
  }
}
CDPA.error <- do.call(rbind, CDPA.error.list) 

all.error <- rbind(
  data.table(algo="CDPA", CDPA.error),
  data.table(algo="PDPA.ignore", PDPA.peaks.error),
  data.table(algo="PDPA.join", PDPA.infeasible.error))

all.error[, table(chunk.name, algo)]

all.totals <- all.error[, list(
  total.fp=sum(fp),
  total.fn=sum(fn),
  total.errors=sum(fp+fn),
  possible.fp=sum(possible.fp),
  possible.fn=sum(possible.tp),
  labels=.N
  ), by=list(algo, chunk.name, sample.id, peaks)]

all.totals[, table(chunk.name, algo)]

all.min <- all.totals[, list(
  min.errors=min(total.errors)
  ), by=list(algo, chunk.name, sample.id)]

all.min[, table(chunk.name, algo)]

min.wide <- dcast(all.min, chunk.name+sample.id~algo)

## as expected there are some (81) problems where the peak joining
## achieves fewer errors than ignoring the infeasible models with
## equality constraints.
min.wide[PDPA.join < PDPA.ignore]
min.wide[PDPA.ignore < PDPA.join]

## Somewhat strangely there are 35 problems when CDPA gets fewer label
## errors, and also 35 problems when PDPA with join when infeasible
## gets fewer label errors.
min.wide[CDPA < PDPA.join]
min.wide[PDPA.join < CDPA]

## Even more strangely the distribution of differences is symmetrical!
min.wide[, table(CDPA-PDPA.join)]

## Not so when ignoring infeasible models -- CDPA is better.
min.wide[, table(CDPA-PDPA.ignore)]

## TODO take a look at the models where CDPA is better.
