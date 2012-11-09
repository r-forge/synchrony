peaks <-
function (t1, t2) {
  f1=find.minmax(t1)
  f2=find.minmax(t2)
  
  common.mins=f1$mins$location %in% f2$mins$location
  common.maxs=f1$maxs$location %in% f2$maxs$location
  peaks=(sum(common.mins)+sum(common.maxs))/sum(max(NROW(f1$mins), 
                                                    NROW(f2$mins)) +
                                                max(NROW(f1$maxs), 
                                                      NROW(f2$maxs)))

  locations=sort(c(f1$mins$location[common.mins], 
                      f1$maxs$location[common.maxs]))
  return (list(peaks=peaks, locations=locations))
}
