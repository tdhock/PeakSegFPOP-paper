library(data.table)

load("cosegData.timings.RData")

txt <- cosegData.timings[, .SD[which.min(seconds),], by=.(lambda, N)]
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
  scale_y_log10(breaks=c(range(cosegData.timings$seconds), 1, 10, 60))+
  scale_x_log10()+
  geom_text(aes(lambda, seconds, label=peaks),
             data=txt,
             vjust=1.5,
             hjust=0)+
  geom_point(aes(lambda, seconds),
             shape=21,
             data=cosegData.timings)


cosegData.timings$data <- "observed"
test.dt <- data.table(
  N=c(1e4, 13623262, 249230603),
  label=c(NA, "typical bedGraph file", "datum for every base"),
  data="extrapolated")


cosegData.timings[, ratio := as.integer(N/lambda)]
txt[, ratio := as.integer(N/lambda)]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ ratio)+
  scale_y_log10(breaks=c(range(cosegData.timings$seconds), 1, 10, 60))+
  scale_x_log10("data points to segment")+
  geom_text(aes(N*1.05, megabytes, label=peaks),
             data=txt,
             vjust=1.5,
             hjust=0)+
  geom_point(aes(N, megabytes, fill=log10(lambda)),
             shape=21,
             data=cosegData.timings)

some.times <- cosegData.timings[ratio==1 & 1e5 <= N,]
fit.mem <- lm(log(megabytes) ~ log(N), some.times)
test.dt$pred.megabytes <- exp(predict(fit.mem, test.dt))
show.dt <- test.dt[-1,]
divide <- 1024
gg.mem <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  scale_y_log10(breaks=as.integer(c(show.dt$pred.megabytes/divide, 1, 4, 12)))+
  scale_x_log10(breaks=c(10^(5:6), show.dt$N), label=scales::comma)+
  geom_line(aes(N, pred.megabytes/divide), data=test.dt)+
  geom_point(aes(N, pred.megabytes/divide, fill=data),
             shape=21,
             data=show.dt)+
  geom_text(aes(N/1.2, pred.megabytes/divide, label=label),
            hjust=1,
            data=show.dt)+
  geom_point(aes(N, megabytes/divide),
             data=some.times)
pdf("figure-cosegData-timings-memory.pdf", h=5)
print(gg.mem)
dev.off()

gg.obs <- ggplot()+
  theme_bw()+
  scale_y_log10(breaks=c(range(cosegData.timings$seconds), 1, 10, 60))+
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
             data=cosegData.timings)
pdf("figure-cosegData-timings-observed.pdf", h=5)
print(gg.obs)
dev.off()

fit <- lm(log(seconds) ~ log(N), cosegData.timings)
test.dt$pred.seconds <- exp(predict(fit, test.dt))
show.dt <- test.dt[-1,]
divide <- 60*60
gg.pred <- ggplot()+
  ggtitle("predicted time to segment human chr1")+
  theme_bw()+
  scale_y_log10("hours of computation", breaks=round(show.dt$pred.seconds/divide, 1))+
  scale_x_log10("data points to segment", breaks=as.integer(show.dt$N), label=scales::comma)+
  geom_line(aes(N, pred.seconds/divide), data=test.dt)+
  geom_point(aes(N, pred.seconds/divide, fill=data),
             shape=21,
             data=test.dt)+
  geom_text(aes(N/1.2, pred.seconds/divide, label=label),
            hjust=1,
            data=show.dt)+
  geom_point(aes(N, seconds/divide, fill=data),
             shape=21,
             data=cosegData.timings)
pdf("figure-cosegData-timings.pdf", h=5)
print(gg.pred)
dev.off()
