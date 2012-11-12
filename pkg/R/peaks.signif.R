peaks.signif <-
function  (t1, t2, nrands=999) {
  observed.peaks=peaks (t1, t2)
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
  return (list(pval=pval, rands=randomized.peaks, peaks=observed.peaks$peaks, 
               locations=observed.peaks$locations))
}
