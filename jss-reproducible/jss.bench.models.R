source("jss-packages.R")

if(!file.exists("jss.bench.models.csv")){
  download.file(
    "https://raw.githubusercontent.com/tdhock/PeakSegFPOP-paper/master/jss.bench.models.csv",
    "jss.bench.models.csv")
}
system("sed -i 's/total.cost/total.loss/' jss.bench.models.csv")
