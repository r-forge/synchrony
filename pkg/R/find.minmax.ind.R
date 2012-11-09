find.minmax.ind <- function (timeseries, mins=FALSE) {
  if (mins)
    ind <- diff(c(Inf, timeseries[,2])) < 0
  else
    ind <- diff(c(-Inf, timeseries[,2])) > 0
  ind <- cumsum(rle(ind)$lengths)
  ind <- ind[seq.int(1, length(ind), 2)]
  ## First and last point cannot be local min/max
  ind = ind[!(ind %in% c(1, NROW(timeseries)))]
  return (ind)
}
