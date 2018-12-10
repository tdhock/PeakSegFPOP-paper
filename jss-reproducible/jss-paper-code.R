### R code from vignette source '/home/tdhock/projects/PeakSegFPOP-paper/jss-paper.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(prompt="R> ")


###################################################
### code chunk number 2: loadData
###################################################

library("PeakSegDisk")
data(Mono27ac, package="PeakSegDisk")
Mono27ac$coverage



###################################################
### code chunk number 3: saveData
###################################################

data.dir <- file.path("Mono27ac", "chr11:60000-580000")
dir.create(data.dir, showWarnings=FALSE, recursive=TRUE)
write.table(Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")



###################################################
### code chunk number 4: problemPeakSegFPOP
###################################################

fit <- PeakSegDisk::problem.PeakSegFPOP(data.dir, "10000")



###################################################
### code chunk number 5: fitLoss
###################################################

fit$loss



###################################################
### code chunk number 6: plotModel
###################################################

library("ggplot2")
gg <- ggplot()+theme_bw()+
  geom_step(aes(chromStart, count), color="grey50", data=Mono27ac$coverage)+
  geom_segment(aes(chromStart, mean, xend=chromEnd, yend=mean),
    color="green", size=1, data=fit$segments)+
  coord_cartesian(xlim=c(2e5, 3e5))
print(gg)



###################################################
### code chunk number 7: seq-search
###################################################

fit <- PeakSegDisk::problem.sequentialSearch(data.dir, 17L)



###################################################
### code chunk number 8: jss-paper.Rnw:1128-1131
###################################################

fit$others[, list(iteration, under, over, penalty, peaks, total.loss)]



