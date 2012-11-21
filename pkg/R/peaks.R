peaks <- function  (t1, t2, nrands = 0) {
  observed.peaks=peaks.aux (t1, t2)
  
  if (nrands == 0) {
    results = observed.peaks    
  }
  else {
    randomized.peaks=numeric(length=nrands+1)
    prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    for (n in 1:nrands) {
      t1.tmp=cbind(t1[,1], sample(t1[,2]))
      t2.tmp=cbind(t2[,1], sample(t2[,2]))
      randomized.peaks[n]=peaks(t1.tmp, t2.tmp)$peaks
      setTxtProgressBar(prog.bar, n)
    }
    randomized.peaks[n+1]=observed.peaks$peaks
    pval=sum(randomized.peaks >= observed.peaks$peaks)/(nrands+1)
    results=list(pval=pval, rands=randomized.peaks, peaks=observed.peaks$peaks, 
                 locations=observed.peaks$locations)
  }
  
  return (results)
}

peaks.aux <- function (t1, t2) {
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
