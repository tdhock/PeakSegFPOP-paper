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
  minimum="min $M_{3,t}$",
  add="$\\ell_t(u_3)=$\n$\\ell(y_t,u_3)$\ncost of data $t$",
  model="$C_{3,t}=$cost\nto data $t$",
  compare="$C^{\\geq}_{2,t}=$cost of\nnon-increasing\nchange after $t$")
cfac <- function(cost.type){
  factor(cost.type, c("model", "compare", "minimum", "add"))
}
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
  geom_point(aes(mean, cost),
             shape=1,
             data=between.intervals[total.segments==3 & timestep %in% ti,])+
  scale_color_discrete("previous\nchange")+
  geom_line(aes(mean, cost, group=paste(data.i.fac, piece.i),
                color=data.i.fac),
            data=cost.lines[total.segments==3 & timestep %in% ti,])+
  ## geom_point(aes(min.cost.mean, min.cost, color=data.i.fac),
  ##            data=minima[total.segments==3 & timestep %in% ti,])+
  coord_cartesian(ylim=c(0.12,0.4), xlim=c(-0.2, 0.8))
step <- function(x="cost up to data $t=35$"){
  factor(x, 
         c("pruning at data $t=34$",
           "cost up to data $t=35$",
           "pruning at data $t=35$"))
}
pstep <- function(x){
  step(paste0("pruning at data $t=", x-1, "$"))
}
l <- function(mean, cost, s, cost.type, label){
  data.table(
    mean, cost, cost.type,
    step=step(s),
    cost.type.fac=cfac(cost.type), label)
}
label.dt <- rbind(
  l(0.5,0.2,pstep(35),"compare","$C^{\\geq}_{2,34}$"),
  l(0,0.3,pstep(35),"model","$C_{3,34}$"),
  l(0.25,0.1,pstep(35),
    "minimum","$M_{3,34}=\\min\\{C_{3,34},C^{\\geq}_{2,34}\\}$"),
  l(0.25, 0.35, step(),
    "model","$C_{3,35}=\\ell_{35}+$\n$M_{3,34}$"),
  l(0.25, 0.15, step(),
    "minimum", "$M_{3,34}$"),
  l(0.4, 0.05, step(),
    "add", "$\\ell_{35}$"),
  l(0.2,0.35,pstep(36),"model","$C_{3,35}$"),
  l(0.6,0.26,pstep(36),"compare","$C^{\\geq}_{2,35}$"),
  l(0.25,0.14,pstep(36),
    "minimum","$M_{3,35}=\\min\\{C_{3,35},C^{\\geq}_{2,35}\\}$"))
cost.env <- envelope[timestep==35 & total.segments==3,]
cost.env$step <- step()
cost.env[, cost.type.fac := cfac("minimum")]
compute.cost <- cost.lines[total.segments==3 & timestep==36 & cost.type=="model",]
compute.cost$step <- step()
add.dt <- data.table(
  getLines(gamma.dt[35,]),
  cost.type.fac=cfac("add"),
  step=step())
envelope[, step := pstep(timestep)]
between.intervals[, step := pstep(timestep)]
cost.lines[, step := pstep(timestep)]
cost.between.dots <- rbind(
  compute.cost[, names(cost.env), with=FALSE],
  cost.env)[, {
    .SD[1,]
  }, by=.(
       step, timestep, minimization, cost.type.fac, piece.i
    )][min.mean != min(data.vec),]
fun.colors <- c(
  model="#E41A1C",
  compare="#377EB8",
  "#4DAF4A",
  add="#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF",
  minimum="grey70")
gg.pruning <- ggplot()+
  coord_cartesian(xlim=c(-0.2, 0), ylim=c(0, 0.3))+
  scale_color_manual(
    "cost type",
    values=fun.colors,
    labels=type.code)+
  scale_fill_manual(values=fun.colors)+
  scale_size_manual(
    "cost type",
    values=c(
      model=1,
      add=1,
      compare=1,
      minimum=2.5),
    labels=type.code)+
  theme_bw()+
  theme(
    panel.grid=element_blank(),
    legend.position="bottom",
    legend.box="horizontal",
    panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ step)+
  ## below cost up to data t=35.
  geom_line(aes(mean, cost,
                color=cost.type.fac,
                size=cost.type.fac),
            data=cost.env)+
  geom_line(aes(mean, cost,
                color=cost.type.fac,
                size=cost.type.fac),
            data=add.dt)+
  geom_line(aes(mean, cost, 
                color=cost.type.fac,
                size=cost.type.fac),
            data=compute.cost)+
  geom_point(aes(mean, cost, fill=cost.type.fac),
             shape=21,
             size=0.5,
             data=cost.between.dots)+
  ## above cost up to data t=35.
  geom_line(aes(mean, cost,
                color=cfac("minimum"),
                size=cfac("minimum")),
            data=envelope[total.segments==3 & timestep %in% ti,])+
  geom_line(aes(mean, cost, 
                color=cost.type.fac,
                size=cost.type.fac),
            data=cost.lines[total.segments==3 & timestep %in% ti,])+
  ## geom_point(aes(min.cost.mean, min.cost, color=data.i.fac),
  ##            data=minima[total.segments==3 & timestep %in% ti,])+
  xlab("mean $u_3$ of segment 3")+
  ylab("cost value $C(u_3)$")+
  coord_cartesian(ylim=c(0,0.4), xlim=c(-0.19, 0.75))+
  guides(color=guide_legend(keyheight=2), fill="none")+
  geom_point(aes(mean, cost, fill=cfac(cost.type)),
             shape=21,
             size=0.5,
             data=between.intervals[total.segments==3 & timestep %in% ti,])+
  geom_text(aes(mean, cost, label=label, color=cost.type.fac),
            size=3.5,
            data=label.dt[cost.type!="minimum",])+
  geom_text(aes(mean, cost, label=label),
            size=3.5,
            color="grey30",
            data=label.dt[cost.type=="minimum",])
tikz("figure-2-min-envelope-slides.tex", 4.8, 3)
print(gg.pruning)
dev.off()
tikz("figure-2-min-envelope.tex", 5, 3)
print(gg.pruning)
dev.off()






