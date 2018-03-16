library(PeakSegOptimal)
library(data.table)
input.dt <- data.table(count=as.integer(c(3, 9, 18, 15, 20, 2)), weight=1)
out.vec <- capture.output({
  fit <- PeakSegOptimal::UnconstrainedPDPA(input.dt$count, max.segments=5L)
})
out.str <- gsub("inf", "Inf", paste(out.vec, collapse="\n"))
it.pattern <- paste0(
  "DP changes=(?<changes>[0-9]+) data_i=(?<data_point>[0-9]+)",
  "(?<rest>",
  "(?:\n.+)+",
  ")")
res <- namedCapture::str_match_all_named(out.str, it.pattern)[[1]]
str(res)
cat(res[1, "rest"])
fun.pattern <- paste0(
  "=(?<fun>.*)",
  "(?<table>",
  "(?:\n[^=]+)+",
  ")")
dt <- data.table(namedCapture::str_match_all_named(res[1, "rest"], fun.pattern)[[1]])
dt[, fread(table), by=list(fun)]
fun.dt <- data.table(res)[, {
  data.table(namedCapture::str_match_all_named(rest, fun.pattern)[[1]])[, fread(table), by=list(fun)]
}, by=list(changes, data_point)]
fwrite(fun.dt, "supplementary.unconstrained.csv")
