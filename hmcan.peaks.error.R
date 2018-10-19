source("packages.R")

files <- Sys.glob("../chip-seq-paper/chunks/H*/*/peaks/hmcan.broad.trained.RData")
## Parse the first occurance of pattern from each of several strings
## using (named) capturing regular expressions, returning a matrix
## (with column names).
str_match_perl <- function(string,pattern){
  stopifnot(is.character(string))
  stopifnot(is.character(pattern))
  stopifnot(length(pattern)==1)
  parsed <- regexpr(pattern,string,perl=TRUE)
  captured.text <- substr(string,parsed,parsed+attr(parsed,"match.length")-1)
  captured.text[captured.text==""] <- NA
  captured.groups <- do.call(rbind,lapply(seq_along(string),function(i){
    st <- attr(parsed,"capture.start")[i,]
    if(is.na(parsed[i]) || parsed[i]==-1)return(rep(NA,length(st)))
    substring(string[i],st,st+attr(parsed,"capture.length")[i,]-1)
  }))
  result <- cbind(captured.text,captured.groups)
  colnames(result) <- c("",attr(parsed,"capture.names"))
  result
}
pattern <-
  paste0("../chip-seq-paper/chunks/",
         "(?<set_name>.+?)",
         "/",
         "(?<chunk_id>[0-9]+)")
matched <- str_match_perl(files, pattern)

hmcan.broad.peaks.error.list <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  chunk.name <- paste0(set.name, "/", chunk.id)
  f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), f))
  load(f)
  regions.RData <- sprintf("../chip-seq-paper/chunks/%s/regions.RData", chunk.name)
  load(regions.RData)
  regions.by.sample <- split(regions, regions$sample.id, drop=TRUE)
  for(param.name in names(peaks)){
    peaks.df <- peaks[[param.name]]
    peaks.by.sample <- split(peaks.df, peaks.df$sample.id)
    for(sample.id in names(regions.by.sample)){
      sample.regions <- regions.by.sample[[sample.id]]
      sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
        peaks.by.sample[[sample.id]]
      }else{
        Peaks()
      }
      err.df <- PeakErrorChrom(sample.peaks, sample.regions)
      hmcan.broad.peaks.error.list[[paste(chunk.name, sample.id, param.name)]] <- 
        data.table(chunk.name, sample.id, param.name, err.df)
    }
  }
}
hmcan.broad.peaks.error <- do.call(rbind, hmcan.broad.peaks.error.list)
save(hmcan.broad.peaks.error, file="hmcan.broad.peaks.error.RData")
