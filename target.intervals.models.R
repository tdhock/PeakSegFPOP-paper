if(!file.exists("labeled_problems_features.csv")){
  download.file(
    "https://raw.githubusercontent.com/tdhock/feature-learning-benchmark/master/labeled_problems_features.csv",
    "labeled_problems_features.csv")
}
if(!file.exists("target.intervals.models.csv")){
  download.file(
    "http://members.cbio.ensmp.fr/~thocking/data/target.intervals.models.csv",
    "target.intervals.models.csv")
}
system("touch target.intervals.models.csv labeled_problems_features.csv")
