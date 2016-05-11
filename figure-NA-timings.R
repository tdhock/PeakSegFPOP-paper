source("packages.R")

load("dp.timings.RData")
load("dp.peaks.NA.RData")

(forward.timings <- data.table(dp.timings))
forward.timings[, chunk.name := paste0(set.name, "/", chunk.id)]

setkey(models.tall, set.name, chunk.name, sample.id)
setkey(forward.timings, set.name, chunk.name, sample.id)
join.dt <- models.tall[forward.timings]
(missing.dt <- join.dt[0 < missing,])
ggplot()+
  ylab("minutes")+
  geom_point(aes(data, seconds/60),
             shape=1,
             data=join.dt)+
  geom_point(aes(data, seconds/60),
             shape=1,
             color="red",
             data=missing.dt)

max.data <- max(missing.dt$data)
forward.small <- join.dt[data <= max.data,]
gg <- ggplot()+
  geom_text(aes(data, seconds, label=missing),
            shape=1,
            data=forward.small)+
  geom_text(aes(data, seconds, label=missing),
            shape=1,
            color="red",
            data=missing.dt)

## TODO: same figure for backward direction.

pdf("figure-NA-timings.pdf")
print(gg)
dev.off()
