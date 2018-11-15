source("jss-packages.R")
job.i <- 1

job.i <- commandArgs(trailingOnly=TRUE)

all.models <- fread("jss.bench.models.several.csv")
job.models <- all.models[job==job.i]
for(i in 1:nrow(job.models)){
  row.i <- job.models[i, row]
  out.csv <- file.path("jss.bench.models.rules", paste0(row.i, ".csv"))
  if(file.exists(out.csv)){
    cat(out.csv, "exists, skipping\n", sep=" ")
  }else{
    cmd <- paste("Rscript jss.bench.models.one.R", row.i)
    status <- system(cmd)
    if(status!=0){
      stop(status)
    }
  }
}
