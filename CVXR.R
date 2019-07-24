works_with_R(
  "3.6.0",
  data.table="1.12.2",
  PeakSegOptimal="2018.5.25",
  CVXR="0.99.6")
data("H3K4me3_XJ_immune_chunk1", package="PeakSegOptimal")
one <- data.table(H3K4me3_XJ_immune_chunk1)[sample.id==sample.id[1]]
plot(coverage ~ chromStart, one)

weight.vec <- with(one, chromEnd-chromStart)
data.vec <- one$coverage

one[, count := as.integer(coverage)]
fpop <- PeakSegOptimal::PeakSegFPOPchrom(one, 10)

N <- length(data.vec)
means <- CVXR::Variable(N)
states <- CVXR::Int(N)
changes <- CVXR::Int(N-1)
list(
  -1 <= changes,
  changes <= 1,
  0 <= states,
  states <= 1)
## how to code graph constraints? maybe convert infeasible situations
## to Inf objective?
changes * diff(means) <= 0
