library(data.table)

load("cosegData.timings.RData")

disk.timings <- cosegData.timings[algorithm=="on.disk",]
R.timings <- cosegData.timings[algorithm=="in.memory",]
txt <- R.timings[, .SD[which.min(seconds),], by=.(lambda, N)]
lab.df <- data.frame(
  seconds=c(1, 60, 60*2),
  label=c("1 second", "1 minute", "2 minutes"))
ggplot()+
  theme_bw()+
  ## scale_x_continuous("log10(lambda penalty parameter)")+
  ## geom_hline(aes(yintercept=log10(seconds)),
  ##            data=lab.df,
  ##            color="grey")+
  ## geom_text(aes(3, log10(seconds), label=label),
  ##           data=lab.df,
  ##           color="grey50",
  ##           vjust=-0.5)+
  scale_y_log10(breaks=c(range(R.timings$seconds), 1, 10, 60))+
  scale_x_log10()+
  geom_text(aes(lambda, seconds, label=peaks),
             data=txt,
             vjust=1.5,
             hjust=0)+
  geom_point(aes(lambda, seconds),
             shape=21,
             data=R.timings)


R.timings$data <- "observed"
test.dt <- data.table(
  algorithm="in.memory",
  N=c(1e4, 4485563, 115591997),
  label=c(NA, "bedGraph lines in biggest\nchr1 region with no gap", "bases in biggest\nregion with no gap"),
  data="extrapolated")

R.timings[, ratio := as.integer(N/lambda)]
txt[, ratio := as.integer(N/lambda)]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ ratio)+
  scale_y_log10()+
  scale_x_log10("data points to segment")+
  geom_text(aes(N*1.05, memory.kilobytes, label=peaks),
             data=txt,
             vjust=1.5,
             hjust=0)+
  geom_point(aes(N, memory.kilobytes, fill=log10(lambda)),
             shape=21,
             data=R.timings)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_y_log10()+
  scale_x_log10("data points to segment")+
  geom_text(aes(N*1.05, memory.kilobytes, label=peaks),
             data=txt,
             vjust=1.5,
             hjust=0)+
  geom_point(aes(N, memory.kilobytes, fill=log10(lambda)),
             shape=21,
             data=R.timings)

some.times <- R.timings[N>1e5,]
fit.mem <- lm(log(memory.kilobytes) ~ log(N), some.times)
test.dt$pred.memory.kilobytes <- exp(predict(fit.mem, test.dt))
show.dt <- test.dt[-1,]
divide <- 1024*1024
gg.mem <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_y_log10("memory usage (gigabytes)", breaks=as.integer(c(show.dt$pred.memory.kilobytes/divide, 1, 4, 12)))+
  scale_x_log10("N = data points to segment", breaks=c(10^(5:6), show.dt$N), label=scales::comma)+
  geom_line(aes(N, pred.memory.kilobytes/divide), data=test.dt)+
  geom_point(aes(N, pred.memory.kilobytes/divide, fill=data),
             shape=21,
             data=show.dt)+
  geom_text(aes(N/2, pred.memory.kilobytes/divide, label=label),
            hjust=1,
            data=show.dt)+
  geom_point(aes(N, memory.kilobytes/divide, fill=data),
             shape=21,
             data=R.timings)
pdf("figure-cosegData-timings-memory.pdf", h=4.5)
print(gg.mem)
dev.off()

gg.mem.disk <- ggplot()+
  theme_bw()+
  scale_color_manual(values=c(in.memory="black", on.disk="red"))+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_y_log10("memory usage (gigabytes)",
                breaks=c(as.integer(show.dt$pred.memory.kilobytes/divide), 1, 0.1, 0.01, 0.001))+
  scale_x_log10("N = data points to segment", breaks=c(10^(5:6), show.dt$N), label=scales::comma)+
  geom_line(aes(N, pred.memory.kilobytes/divide), data=test.dt)+
  geom_point(aes(N, pred.memory.kilobytes/divide, color=algorithm),
             data=show.dt)+
  geom_text(aes(N/2, pred.memory.kilobytes/divide, label=label),
            hjust=1,
            data=show.dt)+
  geom_point(aes(N, memory.kilobytes/divide, color=algorithm),
             data=cosegData.timings)
pdf("figure-cosegData-timings-memory-disk.pdf", h=4.5)
print(gg.mem.disk)
dev.off()

gg.obs <- ggplot()+
  theme_bw()+
  scale_y_log10(breaks=c(range(R.timings$seconds), 1, 10, 60))+
  scale_x_log10("data points to segment")+
  guides(fill=guide_legend())+
  scale_fill_gradient(low="white", high="black", breaks=6:1)+
  geom_text(aes(ifelse(N==max(N), N/1.05, N*1.05), seconds, label=peaks, color=feasible,
                hjust=ifelse(N==max(N), 1, 0)),
            data=txt,
            size=2.5,
            vjust=1.5)+
  geom_point(aes(N, seconds, fill=log10(lambda)),
             shape=21,
             data=R.timings)
pdf("figure-cosegData-timings-observed.pdf", h=5)
print(gg.obs)
dev.off()

divide <- 1024*1024
gg.disk <- ggplot()+
  theme_bw()+
  scale_y_log10("disk usage (gigabytes)",
                breaks=c(0.1, 1, 10, round(max(cosegData.timings$disk.kilobytes/divide), 1)))+
  scale_color_manual(values=c(in.memory="black", on.disk="red"))+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_x_log10("N = data points to segment", breaks=c(10^(5:6), show.dt$N), label=scales::comma)+
  geom_point(aes(N, disk.kilobytes/divide, color=algorithm),
             data=cosegData.timings)
pdf("figure-cosegData-timings-disk.pdf", h=4.5)
print(gg.disk)
dev.off()

fit <- lm(log(seconds) ~ log(N), R.timings)
test.dt$pred.seconds <- exp(predict(fit, test.dt))
show.dt <- test.dt[-1,]
divide <- 60*60
gg.pred <- ggplot()+
  ##ggtitle("predicted time to segment human chr1")+
  theme_bw()+
  scale_y_log10("hours of computation", breaks=round(show.dt$pred.seconds/divide, 2))+
  scale_x_log10("data points to segment",
                breaks=c(1e5, 1e6, as.integer(show.dt$N)), label=scales::comma)+
  geom_line(aes(N, pred.seconds/divide), data=test.dt)+
  geom_point(aes(N, pred.seconds/divide, fill=data),
             shape=21,
             data=test.dt)+
  geom_text(aes(N/2, pred.seconds/divide, label=label),
            hjust=1,
            data=show.dt)+
  geom_point(aes(N, seconds/divide, fill=data),
             shape=21,
             data=R.timings)
pdf("figure-cosegData-timings-pred.pdf", h=5)
print(gg.pred)
dev.off()

divide <- 60
gg.time.disk <- ggplot()+
  theme_bw()+
  scale_y_log10("minutes of computation", breaks=c(0.1, 1, 10, 30))+
  scale_color_manual(values=c(in.memory="black", on.disk="red"))+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_x_log10("N = data points to segment", breaks=c(10^(5:6), show.dt$N), label=scales::comma)+
  geom_point(aes(N, seconds/divide, color=algorithm),
             data=cosegData.timings)
pdf("figure-cosegData-timings.pdf", h=4.5)
print(gg.time.disk)
dev.off()

