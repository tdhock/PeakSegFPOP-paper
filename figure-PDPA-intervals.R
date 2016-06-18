library(ggplot2)
library(data.table)

load("PDPA.intervals.RData")

intervals.dt <- data.table(PDPA.intervals)[segments<data & 1<segments,]
max.dt <- intervals.dt[, .SD[intervals==max(intervals),], by=segments]
text.dt <- max.dt[, .SD[1,], by=segments]

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ .)+
  geom_line(aes(data, intervals), data=intervals.dt)+
  geom_point(aes(data, intervals), data=max.dt, color="red")+
  geom_text(aes(0, intervals, label=intervals),
            color="red",
            data=text.dt,
            hjust=1.1,
            vjust=1)+
  scale_y_continuous(
    "intervals stored by the constrained optimal segmentation algorithm",
    breaks=seq(0, 40, by=20))

png("figure-PDPA-intervals.png", w=7, h=6, units="in", res=200)
print(gg)
dev.off()
