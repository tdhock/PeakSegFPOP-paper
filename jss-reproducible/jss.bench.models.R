source("jss-packages.R")

if(!file.exists("jss.bench.models.csv")){
  download.file(
    "https://raw.githubusercontent.com/tdhock/PeakSegFPOP-paper/master/jss.bench.models.csv.gz",
    "jss.bench.models.csv.gz")
  system("gunzip jss.bench.models.csv.gz")
}
system("sed -i 's/total.cost/total.loss/' jss.bench.models.csv")
