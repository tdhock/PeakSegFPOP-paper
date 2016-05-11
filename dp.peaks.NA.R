source("packages.R")

load("dp.peaks.matrices.RData")

model.info <- list()
for(set.name in names(dp.peaks.matrices)){
  chunk.list <- dp.peaks.matrices[[set.name]]
  for(chunk.name in names(chunk.list)){
    chunk.info <- chunk.list[[chunk.name]]
    regions <- chunk.info$regions
    na.mat <- is.na(chunk.info$PeakSeg)
    na.models <- rowSums(na.mat)
    min.err <- apply(chunk.info$PeakSeg, 1, min, na.rm=TRUE)
    max.peaks <- apply(chunk.info$PeakSeg, 1, paste, collapse=" ")
    model.info[[paste(set.name, chunk.name)]] <- 
      data.table(set.name, chunk.name,
                 sample.id=names(regions), regions,
                 missing=na.models,
                 min.err, max.peaks)
  }
}
models.tall <- do.call(rbind, model.info)
models.tall[, list(
  regions=sum(regions),
  missing=sum(missing),
  problems=.N,
  all10=sum(missing == 0))]
## 14 problems with less than 10 models... there are 3 for which there
## is some error.
models.tall[missing != 0 & min.err != 0,]
dp.peaks.matrices$H3K4me3_XJ_immune[["H3K4me3_XJ_immune/2"]]$PeakSeg["McGill0028", ]
14/2752*100 # 0.5% of problems for which we did not compute all 10 models.
3/2752*100 # 0.1% of problems for which this actually caused an error.
models.tall[missing != 0 & !grepl("NA$", max.peaks),]
5/2752*100 # 0.2% could compute 9 peaks model but not a smaller model.
## TODO: plot these data/models!

dp.peaks.NA <- models.tall[missing != 0,]
save(dp.peaks.NA, models.tall, file="dp.peaks.NA.RData")
