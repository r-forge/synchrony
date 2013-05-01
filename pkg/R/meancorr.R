## Community matrix comm.matrix: n x m matrix with n=time step, m=species
meancorr <- function (comm.matrix, nrands = 0, 
                            alternative=c("one.tailed", "two.tailed"), 
                            method=c("pearson", "kendall", "spearman"), ...) {
  comm.matrix=as.matrix(comm.matrix)
  results=list()
  methods=c("pearson", "kendall", "spearman")
  method=match.arg(method, methods)
  
  tails=c("one.tailed", "two.tailed")
  alternative=match.arg(tolower(alternative), tails)
  
  results$obs=meancorr.aux (comm.matrix, method=method, ...)
  
  if (nrands > 0) {
    nr=NROW(comm.matrix)
    nc=NCOL(comm.matrix)      
    prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    results$rands=numeric(length=nrands+1)*NA
    for (i in 1:nrands) {
      rand.mat=apply(comm.matrix, 2, sample)
      results$rands[i]=meancorr.aux(rand.mat, method=method, ...)
      setTxtProgressBar(prog.bar, i)
    }
    results$rands[nrands+1]=results$obs
    
    if (alternative == "two.tailed") {
      pvals=sum(abs(results$rands) >= abs(results$obs)/(nrands+1))
    }
    else {
      pvals=ifelse (results$obs > 0, 
                    sum(results$rands >= results$obs)/(nrands+1), 
                    sum(results$rands <= results$obs)/(nrands+1))
    }
    
    results$pval=sum(results$rands >= results$obs)/(nrands+1)
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