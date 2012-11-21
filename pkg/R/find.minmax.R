find.minmax <- function (timeseries) {
  ## Find local maxima
  max.index <- find.minmax.aux (timeseries)
  ## Find local minima
  min.index <- find.minmax.aux (timeseries, mins=TRUE)
  
  mins=as.data.frame(timeseries[min.index,])
  maxs=as.data.frame(timeseries[max.index,])
  
  if (NCOL(mins)==1)
    mins=t(mins)
  if (NCOL(maxs)==1)
    maxs=t(maxs)
  
  col.names=c("location", "val")
  colnames(mins)=col.names
  rownames(mins)=1:NROW(mins)
  colnames(maxs)=col.names
  rownames(maxs)=1:NROW(maxs)

  return (list(mins=mins, maxs=maxs))
}

find.minmax.aux <- function (timeseries, mins=FALSE) {
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
