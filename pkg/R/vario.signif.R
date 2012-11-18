vario.signif <- function (nrands=999, bins, all.combs, grpdata, data, data2=data, type, 
                          vario, glob.mean=NULL, glob.sd=NULL, glob.N=NULL, is.multivar=FALSE,
                          alternative=c("one.tailed", "two.tailed"), centered, regional.mean,
                          mult.test.corr=c(FALSE, "holm", "hochberg", "sidak", "bonferroni")) {
  
  rands=matrix(NA, nrow=nrands+1, ncol=length(bins))
  prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
  n.cols=NCOL(data)
  
  if (is.null(data2))
    data2=data
  for (i in 1:nrands) {
    s=sample(grpdata)
    
    if (is.multivar==TRUE) {
      vals=data
      rands[i,]=tapply(vals, s, FUN=mean, na.rm=TRUE)
    }
    else {
      for (j in 1:(length(bins))) {        
        tmp=all.combs[s==j,]
        tmp=tmp[complete.cases(tmp),]
        rands[i,j]=vario.func(data[tmp[,1], 3], data2[tmp[,2], 3], glob.mean, 
                              glob.sd, glob.N, is.multivar, type)
      }
    }
    setTxtProgressBar(prog.bar, i)
  }
  rands[nrands+1,]=vario
  rands=rands-rowMeans(rands, na.rm=TRUE)
  crit.val=0
  if (!centered) {
    rands=rands+regional.mean
    crit.val=regional.mean
  }
  
  if (alternative == "two.tailed") {
    pvals=apply(rands, MARGIN=2, 
                function (x) {
                  sum(abs(x - crit.val) >= abs(x[nrands+1] - crit.val))/(nrands+1)
                })
  }
  else {
    pvals=apply(rands, MARGIN=2, function (x) {
      ifelse (x[nrands+1] > crit.val, 
              sum(x >= x[nrands+1])/(nrands+1), 
              sum(x <= x[nrands+1])/(nrands+1))})
  }

  if (mult.test.corr[1] !=FALSE) {
    pvals=p.adjust(pvals, method=mult.test.corr[1])
  }

  colnames(rands)=names(bins)
  names(pvals)=names(bins)
  
  results=list(rands=rands, bins=bins, pvals=pvals, vario=vario, mult.test.corr=mult.test.corr[1])
}