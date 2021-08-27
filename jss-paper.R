### R code from vignette source 'C:/Users/th798/projects/PeakSegFPOP-paper/jss-paper.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(prompt="R> ", datatable.print.trunc.cols=FALSE)


###################################################
### code chunk number 2: loadData
###################################################

library("PeakSegDisk")
data("Mono27ac", package="PeakSegDisk")
Mono27ac$coverage


###################################################
### code chunk number 3: saveData
###################################################

data.dir <- file.path("Mono27ac", "chr11-60000-580000")
dir.create(data.dir, showWarnings = FALSE, recursive = TRUE)
PeakSegDisk::writeBedGraph(
  Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"))


###################################################
### code chunk number 4: problemPeakSegFPOP
###################################################

fit <- PeakSegDisk::PeakSegFPOP_dir(data.dir, 10000.5)


###################################################
### code chunk number 5: fitLoss
###################################################

summary(fit)


###################################################
### code chunk number 6: plotModel
###################################################

library("ggplot2")
gg <- ggplot() + theme_bw() + coord_cartesian(xlim = c(2e5, 3e5)) +
  geom_step(aes(chromStart, count), color = "grey50",
    data = Mono27ac$coverage) +
  geom_segment(aes(chromStart, mean, xend = chromEnd, yend = mean),
    color = "green", size = 1, data = fit$segments) +
  xlab("Position on chromosome (bases)")
print(gg)


###################################################
### code chunk number 7: plot2
###################################################

fit <- PeakSegDisk::PeakSegFPOP_df(Mono27ac$coverage, 999.9)
class(fit)
gg <- plot(fit)
print(gg)


###################################################
### code chunk number 8: zoom
###################################################

print(gg + ggplot2::coord_cartesian(xlim = c(205000, 210000)))


###################################################
### code chunk number 9: seq-search
###################################################

fit <- PeakSegDisk::sequentialSearch_dir(data.dir, 17L)


###################################################
### code chunk number 10: jss-paper.Rnw:1247-1249
###################################################

fit$others[, list(iteration, under, over, penalty, peaks, total.loss)]


