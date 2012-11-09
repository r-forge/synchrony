## Community matrix comm.matrix: n x m matrix with n=time step, m=species
community.sync.signif <- function (comm.matrix, nrands = 999) {
  comm.matrix=as.matrix(comm.matrix)
  results=list()
  results$obs=community.sync(comm.matrix)
  nr=NROW(comm.matrix)
  nc=NCOL(comm.matrix)
  
  if (nrands==0) {
    results$rands=NA
    results$pval=NA
  }
  else {
    prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    results$rands=numeric(length=nrands+1)*NA
    for (i in 1:nrands) {
      rand.mat=matrix(sample(as.numeric(comm.matrix)), nrow=nr, ncol=nc)
      results$rands[i]=community.sync(rand.mat)
      setTxtProgressBar(prog.bar, i)
    }  
    results$rands[nrands+1]=results$obs
    results$pval=sum(results$rands >= results$obs)/(nrands+1)
  }
  
  return(results)
}