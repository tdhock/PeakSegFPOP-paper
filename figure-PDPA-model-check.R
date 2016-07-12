source("packages.R")

load("PDPA.model.check.RData")

## Models not computable by the cDPA.
PDPA.model.check[, best.cDPA.loss := ifelse({
  is.na(fwd.loss)
}, rev.loss, ifelse(
  is.na(rev.loss), fwd.loss, ifelse(
    fwd.loss<rev.loss, fwd.loss, rev.loss)))]
PDPA.model.check[is.na(fwd.loss),]
PDPA.model.check[is.na(rev.loss),]
PDPA.model.check[is.na(best.cDPA.loss),]

## More likely models found by the PDPA.
PDPA.model.check[, pdpa.better := pdpa.loss-best.cDPA.loss < -1e-8]
abline.dt <- data.table(slope=1, intercept=0)
ggplot()+
  theme_bw()+
  geom_abline(aes(slope=slope, intercept=intercept),
              color="grey",
              data=abline.dt)+
  coord_equal()+
  geom_point(aes(best.cDPA.loss, pdpa.loss, color=pdpa.feasible),
             shape=21,
             data=PDPA.model.check)

## more likely feasible models.
sample.counts <- PDPA.model.check[, list(
  better.models=sum(ifelse(
    is.na(best.cDPA.loss), TRUE, pdpa.loss-best.cDPA.loss < -1e-8 & pdpa.feasible)),
  cDPA.models=sum(!is.na(best.cDPA.loss)),
  PDPA.models=sum(pdpa.feasible)
  ), by=.(chunk, sample.id)]
better.counts <- sample.counts[0<better.models,][order(better.models),]

one.sample <- PDPA.model.check[chunk=="H3K4me3_TDH_other/17" & sample.id=="McGill0016",]
ggplot()+
  geom_line(aes(peaks, pdpa.loss),
            data=one.sample)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=5))+
  geom_point(aes(peaks, pdpa.loss, size=pdpa.feasible),
             data=one.sample)+
  geom_point(aes(peaks, best.cDPA.loss, fill=pdpa.better),
             shape=21,
             data=one.sample)

better.counts[PDPA.models!=10,]
one.sample <- PDPA.model.check[chunk=="H3K36me3_AM_immune/7" & sample.id=="McGill0079",]
ggplot()+
  geom_line(aes(peaks, pdpa.loss),
            data=one.sample)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=5))+
  geom_point(aes(peaks, pdpa.loss, size=pdpa.feasible),
             data=one.sample)+
  geom_point(aes(peaks, best.cDPA.loss, fill=pdpa.better),
             shape=21,
             data=one.sample)

one.sample <- PDPA.model.check[chunk=="H3K4me3_TDH_other/21" & sample.id=="McGill0267",]
ggplot()+
  geom_line(aes(peaks, pdpa.loss),
            data=one.sample)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=5))+
  geom_point(aes(peaks, pdpa.loss, size=pdpa.feasible),
             data=one.sample)+
  geom_point(aes(peaks, best.cDPA.loss, fill=pdpa.better),
             shape=21,
             data=one.sample)

sample.counts[order(cDPA.models),]
one.sample <- PDPA.model.check[chunk=="H3K36me3_AM_immune/19" & sample.id=="McGill0104",]
ggplot()+
  geom_line(aes(peaks, pdpa.loss),
            data=one.sample)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=5))+
  geom_point(aes(peaks, pdpa.loss, size=pdpa.feasible),
             data=one.sample)+
  geom_point(aes(peaks, best.cDPA.loss, fill=pdpa.better),
             shape=21,
             data=one.sample)

sample.counts[PDPA.models < cDPA.models,]
one.sample <- PDPA.model.check[chunk=="H3K36me3_AM_immune/13" & sample.id=="McGill0079",]
ggplot()+
  geom_line(aes(peaks, pdpa.loss),
            data=one.sample)+
  scale_size_manual(values=c("TRUE"=1, "FALSE"=5))+
  geom_point(aes(peaks, pdpa.loss, size=pdpa.feasible),
             data=one.sample)+
  geom_point(aes(peaks, best.cDPA.loss, fill=pdpa.better),
             shape=21,
             data=one.sample)
