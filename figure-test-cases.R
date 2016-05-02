source("packages.R")

test.case.list <- list(
  undefined3=list(
    data.vec=c(1, 2, 3),
    ends=rbind(
  undefined4=list(
    data.vec=c(1, 10, 14, 13),
    ends=rbind(
  undefined4rev=list(
    data.vec=c(13, 14, 10, 1),
    ends=rbind(
  ok=list(
    data.vec=c(1, 2, 10, 14, 3, 5),
    ends=rbind(

for(test.name in names(test.case.list)){
  test.case <- test.case.list[[test.name]]
  cDPA.result <- cDPA(test.case$data.vec, maxSegments=3)
  print(test.name)
  dput(cDPA.result$ends)
}
