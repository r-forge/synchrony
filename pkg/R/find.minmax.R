find.minmax <- function (timeseries) {
  ## Find local maxima
  max.index <- find.minmax.ind (timeseries)
  ## Find local minima
  min.index <- find.minmax.ind (timeseries, mins=TRUE)
  
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
