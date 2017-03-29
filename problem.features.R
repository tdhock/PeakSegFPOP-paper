source("packages.R")

files <- Sys.glob("data/H*/*/counts.RData")
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
  paste0("data/",
         "(?<set_name>.+?)",
         "/",
         "(?<chunk_id>[0-9]+)")
matched <- str_match_perl(files, pattern)
problem.features <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  regions.str <- paste0(set.name, "/", chunk.id)
  count.f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), count.f))
  load(count.f)
  count.list <- split(counts, counts$sample.id)
  sample.ids <- names(count.list)
  chunk.mat <-
    matrix(NA, length(sample.ids), 14*4,
           dimnames=list(sample.id=sample.ids, feature=NULL))
  for(sample.id in sample.ids){
    count.df <- count.list[[sample.id]]
    bases <- with(count.df, chromEnd-chromStart)
    long <- rep(count.df$coverage, bases)
    n.bases <- sum(bases)
    n.data <- nrow(count.df)
    count.df <- count.list[[sample.id]]
    feature.vec <-
      c(unweighted.quartile=quantile(count.df$coverage),
        weighted.quartile=quantile(long),
        unweighted.mean=mean(count.df$coverage),
        weighted.mean=mean(long),
        bases=n.bases,
        ## sqrt=under.sqrt,
        ## square=in.square,
        ## cleynen=cleynen,
        ## p.err,
        data=n.data)        
    log.features <-
      c(feature.vec,
        `log+1`=log(feature.vec+1),
        log=log(feature.vec),
        log.log=log(log(feature.vec)))
    chunk.mat[sample.id, ] <- log.features
  }
  colnames(chunk.mat) <- names(log.features)
  colnames(chunk.mat)[colMeans(is.finite(chunk.mat)) == 1]
  problem.features[[regions.str]] <- chunk.mat
}

save(problem.features, file="problem.features.RData")

