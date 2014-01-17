## Community matrix comm.matrix: n x m matrix with n=time step, m=species
meancorr <- function (comm.matrix, nrands = 0, 
                            alternative=c("two.tailed", "greater", "less"), 
                            method=c("pearson", "kendall", "spearman"), 
                      type=1, quiet = FALSE, ...) {
  comm.matrix=as.matrix(comm.matrix)
  results=list()
  methods=c("pearson", "kendall", "spearman")
  method=match.arg(method, methods)
  
  alternatives=c("two.tailed", "greater", "less")
  alternative=match.arg(tolower(alternative), alternatives)
  
  results$obs=meancorr.aux (comm.matrix, method=method, ...)
  
  if (nrands > 0) {
    nr=NROW(comm.matrix)
    nc=NCOL(comm.matrix)
    if (!quiet)
      prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    results$rands=numeric(length=nrands+1)*NA
    for (i in 1:nrands) {
      if (type==1)
        rand.mat=apply(comm.matrix, 2, sample)
      else {
        lags=sample(1:nr, size=nc, replace=TRUE)
        rand.mat=mlag(comm.matrix, lags)
      }
      results$rands[i]=meancorr.aux(rand.mat, method=method, ...)
      if (!quiet)
        setTxtProgressBar(prog.bar, i)
    }
    results$rands[nrands+1]=results$obs
    
    if (alternative == "two.tailed") {
      pvals=sum(abs(results$rands) >= abs(results$obs))/(nrands+1)
    }
    else {
      if (alternative=="greater")
        pvals=sum(results$rands >= results$obs)/(nrands+1)
      else
        pvals=sum(results$rands <= results$obs)/(nrands+1)
    }
    
    results$pval=pvals
    results$alternative=alternative
  }
  results$method=method
  class(results)="synchrony"
  return (results)
}

meancorr.aux <- function (data, method=method, ...) {
  mean.corr=suppressWarnings(cor(data, method=method, ...))
  mean.corr=mean(mean.corr[lower.tri(mean.corr)], na.rm=TRUE)
  return (mean.corr)
}