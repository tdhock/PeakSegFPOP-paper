source("packages.R")

files <- Sys.glob("../chip-seq-paper/chunks/H*/*/Segmentor.model.RData")
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

Segmentor.peaks.error.list <- list()
Segmentor.peaks.list <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  chunk.name <- paste0(set.name, "/", chunk.id)
  f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), f))
  load(f)
  Segmentor.model.list[[chunk.name]] <- Segmentor.model
  regions.RData <- sprintf("../chip-seq-paper/chunks/%s/regions.RData", chunk.name)
  load(regions.RData)
  regions.by.sample <- split(regions, regions$sample.id, drop=TRUE)
  for(sample.id in names(Segmentor.model)){
    sample.regions <- regions.by.sample[[sample.id]]
    sample.peaks <- Segmentor.model[[sample.id]]$peaks
    for(peaks.str in names(sample.peaks)){
      peak.df <- sample.peaks[[peaks.str]]
      if(is.null(peak.df)){
        cat(paste(chunk.name, sample.id, peaks.str, "peaks\n"))
        peak.df <- Peaks()
      }
      err.df <- PeakErrorChrom(peak.df, sample.regions)
      Segmentor.peaks.error.list[[paste(chunk.name, sample.id, peaks.str)]] <- 
        data.table(chunk.name, sample.id, peaks=peaks.str, err.df)
    }
  }
}
Segmentor.peaks.error <- do.call(rbind, Segmentor.peaks.error.list)
save(Segmentor.model.list, Segmentor.peaks.error, file="Segmentor.peaks.error.RData")
