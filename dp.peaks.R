source("packages.R")

files <- Sys.glob("../chip-seq-paper/chunks/H*/*/dp.model.RData")
pattern <- paste0(
  "chunks/",
  "(?<set_name>.+?)",
  "/",
  "(?<chunk_id>[0-9]+)")
(matched <- str_match_named(files, pattern))

dp.peaks <- list()
dp.loss.list <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  regions.str <- paste0(set.name, "/", chunk.id)
  f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), f))
  load(f)
  for(sample.id in names(dp.model)){
    sample.list <- dp.model[[sample.id]]
    dp.loss.list[[paste(regions.str, sample.id)]] <- data.table(
      set.name, chunk.id, chunk.name=regions.str, sample.id,
      sample.list$error)
    dp.peaks[[regions.str]][[sample.id]] <- sample.list$peaks
  }
}
dp.loss <- do.call(rbind, dp.loss.list)

save(dp.loss, dp.peaks, file="dp.peaks.RData")
