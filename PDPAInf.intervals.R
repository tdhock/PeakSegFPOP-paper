source("packages.R")

library(data.table)
library(coseg)

PDPA.intervals.list <- list()
for(algo.type in c("PDPA", "PDPAInf")){
  model.file.vec <- Sys.glob(paste0("data/H*/*/", algo.type, ".model.RData"))
  file.i <- 42
  sample.i <- 24
  file.i.vec <- seq_along(model.file.vec)
  for(file.i in file.i.vec){
    model.file <- model.file.vec[[file.i]]
    chunk.id.dir <- dirname(paste(model.file))
    chunk.id <- basename(chunk.id.dir)
    set.dir <- dirname(chunk.id.dir)
    set.name <- basename(set.dir)
    cat(sprintf("%4d / %4d %s\n", file.i, length(model.file.vec), model.file))
    if(!file.exists(model.file)){
      stop("run PDPA.timings.R first")
    }else{
      obj <- load(model.file) #counts
      model.list <- get(obj)
      for(sample.id in names(model.list)){
        model <- model.list[[sample.id]]
        PDPA.intervals.list[[paste(algo.type, model.file, sample.id)]] <- data.table(
          algo.type, set.name, chunk.id, sample.id,
          mean.intervals=model$mean.intervals,
          max.intervals=model$max.intervals)
      }
    }
  }
}
PDPAInf.intervals <- do.call(rbind, PDPA.intervals.list)

save(PDPAInf.intervals, file="PDPAInf.intervals.RData")
