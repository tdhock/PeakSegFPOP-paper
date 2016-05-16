source("packages.R")

load("dp.timings.RData")

PDPA.model.RData.vec <- Sys.glob("data/H*/*/PDPA.model/*.RData")
PDPA.timings.list <- list()
for(PDPA.model.RData.i in seq_along(PDPA.model.RData.vec)){
  PDPA.model.RData <- PDPA.model.RData.vec[[PDPA.model.RData.i]]
  load(PDPA.model.RData)
  PDPA.timings.list[[PDPA.model.RData]] <- data.table(
    result$timing,
    result$model$intervals[, list(
      mean.intervals=mean(intervals))])
}
PDPA.timings <- do.call(rbind, PDPA.timings.list)

gg <- ggplot()+
  ylab("hours")+
  geom_point(aes(data, seconds/60/60, color=mean.intervals),
             shape=1,
             data=PDPA.timings)

ggplot()+
  ylab("hours")+
  geom_point(aes(data, mean.intervals),
             shape=1,
             data=PDPA.timings)

pdf("figure-PDPA-timings.pdf")
print(gg)
dev.off()
