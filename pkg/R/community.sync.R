## Community matrix comm.matrix: n x m matrix with n=time step, m=species
community.sync <- function (comm.matrix, nrands = 0, method=c("pearson", "kendall", "spearman"), 
                            alternative=c("greater", "less"), type=1, quiet=FALSE, ...) {
  alternatives=c("greater", "less")
  alternative=match.arg(tolower(alternative), alternatives)
  
  comm.matrix=as.matrix(comm.matrix)
  results=list()
  results$obs=community.sync.aux (comm.matrix)
  results$meancorr=meancorr(comm.matrix, method=method, ...)$obs
  
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
      results$rands[i]=community.sync.aux(rand.mat)
      if (!quiet)
        setTxtProgressBar(prog.bar, i)
    }
    results$rands[nrands+1]=results$obs
    if (alternative=="greater")
      results$pval=sum(results$rands >= results$obs)/(nrands+1)
    else
      results$pval=sum(results$rands <= results$obs)/(nrands+1)
    results$alternative=alternative
  }
  class(results)="synchrony"
  return (results)
}

community.sync.aux <- function (comm.matrix) {
  species.sd=apply(comm.matrix, MARGIN=2, FUN=sd)
  community.var=var(rowSums(comm.matrix))
  return(community.var/sum(species.sd, na.rm=TRUE)^2)
}
