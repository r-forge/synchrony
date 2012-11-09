phase.sync.signif <- function (t1, t2, nrands = 999, type = 1, nbreaks = 10, 
                                          mins = FALSE) {
  
  p=phase.sync(t1, t2, mins=mins)
  rands=numeric(length=nrands+1)*NA

  if (type == 1) {
    column=3
    breaks=seq(from=0, to=2*pi, length.out=10)
  }
  else {
    column=4
    breaks=seq(from=-pi, to=pi, length.out=10)
  }
  h.obs=hist(p$deltaphase[, column], breaks=breaks, plot=FALSE)
  p.obs=h.obs$counts/sum(h.obs$counts)
  nbins.obs=length(h.obs$counts)
  
  S.obs=-sum(p.obs*log(p.obs), na.rm=TRUE)
  Smax.obs=log(nbins.obs)
  Q.obs=(Smax.obs-S.obs)/Smax.obs
  
  ## Determine transition probabilities
  distr.t1 <- cut(t1[,2], quantile(t1[,2], seq(0, 1, len = (nbreaks+1))), 
                  include.lowest = TRUE, labels=FALSE)
  distr.t2 <- cut(t2[,2], quantile(t2[,2], seq(0, 1, len = (nbreaks+1))), 
                  include.lowest = TRUE, labels=FALSE)  
  trans.t1=matrix(nrow=nbreaks, ncol=nbreaks, 0)
  trans.t2=matrix(nrow=nbreaks, ncol=nbreaks, 0)
  
  for (i in 1:(NROW(t1)-1)) {
    trans.t1[distr.t1[i], distr.t1[i+1]]=trans.t1[distr.t1[i], distr.t1[i+1]]+1      
    trans.t2[distr.t2[i], distr.t2[i+1]]=trans.t2[distr.t2[i], distr.t2[i+1]]+1      
  }
  trans.t1=trans.t1/rowSums(trans.t1)  
  trans.t2=trans.t2/rowSums(trans.t2)  
  
  prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
  for (r in 1:nrands) {
    surr.t1=surrogate.ts(ts=t1, distr.ts=distr.t1, nbreaks=nbreaks)
    surr.t2=surrogate.ts(ts=t2, distr.ts=distr.t2, nbreaks=nbreaks)
    
    p.rand=phase.sync(surr.t1$surr.ts, surr.t2$surr.ts)
    
    rand.h=hist(p.rand$deltaphase[, column], breaks=breaks, plot=FALSE)
    rand.p=rand.h$counts/sum(rand.h$counts, na.rm=TRUE)
    rand.nbins=length(rand.h$counts)  
    rand.S=-sum(rand.p*log(rand.p), na.rm=TRUE)
    rand.Smax=log(rand.nbins)
    rands[r]=(rand.Smax-rand.S)/rand.Smax
    setTxtProgressBar(prog.bar, r)
  }
  rands[r+1]=Q.obs
  pValue = sum (rands >= Q.obs)/(nrands+1)
  
  ## ICDF
  o=sort(rands)
  icdf=data.frame(Q=o, icdf=sapply(o, FUN=function (x) {sum(rands >= x)/(nrands+1)}))
  
  return (list(Q.obs=Q.obs, pval=pValue, rands=rands, phases1=p$phases1,
               phases2=p$phases2, deltaphase=p$deltaphase, icdf=icdf))
}
