source("packages.R")

load("all.cv.RData")

roc.tall <- melt(
  roc.thresh[threshold=="predicted"],
  measure.vars=c("auc", "accuracy.percent"),
  id.vars=c("set.name", "fold.i", "algo", "validation.crit"))
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(variable ~ set.name, scales="free")+
  geom_point(aes(
    algo, value, color=validation.crit),
    shape=1,
    data=roc.tall)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(set.name ~ variable, scales="free")+
  geom_point(aes(
    value, algo),
    data=roc.tall)

save(roc, roc.thresh, file="all.cv.RData")
