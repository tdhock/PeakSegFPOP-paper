source("packages.R")

load("figure-constrained-PDPA-normal-real.RData")

cost.lines[, data.i.fac := factor(data.i)]
cost.lines[, minimization := paste(
  total.segments, "segments up to data point", timestep)]
left.of.intervals <-
  cost.lines[, .SD[1,],
             by=.(total.segments, timestep, minimization, cost.type, piece.i)]
between.intervals <- left.of.intervals[min.mean != min(data.vec),]

type.code <- c(
  model="min cost\nin 3 segments\nup to data $t-1$\n$C_{3,t-1}$",
  compare="cost\nof change\nafter $t-1$\n$C^{\\geq}_{2,t-1}$")
cfac <- function(cost.type)factor(cost.type, c("model", "compare", "minimum"))
cost.lines[, cost.type.fac := cfac(cost.type)]
ti <- 5:7
ti <- 35:36
gg.pruning <- ggplot()+
  ##coord_cartesian(xlim=c(-0.2, 0), ylim=c(0, 0.3))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(timestep ~ total.segments, scales="free",
             labeller=function(var, val){
               if(var %in% c("total.segments", "timestep")){
                 paste(var, "=", val)
               }else{
                 paste(val)
               }
             })+
  geom_line(aes(mean, cost),
            color="grey",
            size=2,
            data=envelope[total.segments==3 & timestep %in% ti,])+
  geom_line(aes(mean, cost, group=paste(data.i.fac, piece.i),
                color=data.i.fac),
            data=cost.lines[total.segments==3 & timestep %in% ti,])+
  ## geom_point(aes(min.cost.mean, min.cost, color=data.i.fac),
  ##            data=minima[total.segments==3 & timestep %in% ti,])+
  coord_cartesian(ylim=c(0.1,0.4))
gg.pruning <- ggplot()+
  ##coord_cartesian(xlim=c(-0.2, 0), ylim=c(0, 0.3))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ timestep, 
             labeller=function(var, val){
               if(var %in% c("total.segments", "timestep")){
                 paste(var, "=", val)
               }else{
                 paste(val)
               }
             })+
  geom_line(aes(mean, cost),
            color="grey",
            size=2,
            data=envelope[total.segments==3 & timestep %in% ti,])+
  geom_point(aes(mean, cost,
                 key=mean,
                 showSelected=total.segments, showSelected2=timestep),
             shape=1,
             data=between.intervals[total.segments==3 & timestep %in% ti,])+
  scale_color_discrete("previous\nchange")+
  geom_line(aes(mean, cost, group=paste(data.i.fac, piece.i),
                color=data.i.fac),
            data=cost.lines[total.segments==3 & timestep %in% ti,])+
  ## geom_point(aes(min.cost.mean, min.cost, color=data.i.fac),
  ##            data=minima[total.segments==3 & timestep %in% ti,])+
  coord_cartesian(ylim=c(0.12,0.4), xlim=c(-0.2, 0.8))
l <- function(mean, cost, timestep, cost.type, label){
  data.table(
    mean, cost, timestep, cost.type,
    cost.type.fac=cfac(cost.type), label)
}
label.dt <- rbind(
  l(0.5,0.2,35,"compare","$C^{\\geq}_{2,34}$"),
  l(0,0.3,35,"model","$C_{3,34}$"),
  l(0.2,0.35,36,"model","$C_{3,35}$"),
  l(0.25,0.15,36,"compare","$C^{\\geq}_{2,35}$"))
gg.pruning <- ggplot()+
  coord_cartesian(xlim=c(-0.2, 0), ylim=c(0, 0.3))+
  scale_color_manual(
    "cost type",
    ##breaks=c("model", "compare", "minimum"),
    values=c(
    model="#E41A1C", compare="#377EB8",
    "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
    "#A65628", "#F781BF",
    minimum="#999999"), labels=type.code)+
  scale_size_manual(
    "cost type",
    ##breaks=c("model", "compare", "minimum"),
    values=c(
    model=1,
    compare=1,
    minimum=3), labels=type.code)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ timestep, 
             labeller=function(var, val){
               paste("pruning at data $t=", val, "$")
             })+
  geom_line(aes(mean, cost,
                color=cfac("minimum"),
                size=cfac("minimum")),
            data=envelope[total.segments==3 & timestep %in% ti,])+
  geom_point(aes(mean, cost,
                 key=mean,
                 showSelected=total.segments, showSelected2=timestep),
             shape=1,
             data=between.intervals[total.segments==3 & timestep %in% ti,])+
  geom_line(aes(mean, cost, 
                color=cost.type.fac,
                size=cost.type.fac),
            data=cost.lines[total.segments==3 & timestep %in% ti,])+
  ## geom_point(aes(min.cost.mean, min.cost, color=data.i.fac),
  ##            data=minima[total.segments==3 & timestep %in% ti,])+
  xlab("mean $\\mu$")+
  ylab("cost value $C(\\mu)$")+
  coord_cartesian(ylim=c(0.12,0.4), xlim=c(-0.19, 0.79))+
  guides(color=guide_legend(keyheight=4))+
  geom_text(aes(mean, cost, label=label, color=cost.type.fac),
            data=label.dt)
tikz("figure-2-min-envelope.tex", 6, 2.5)
print(gg.pruning)
dev.off()

