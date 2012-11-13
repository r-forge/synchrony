correlated.matrix <- function (rho = 0, sigma=1, mu=0, ntimes = 200, nspecies = 10) {
  U=matrix(rho, nrow=nspecies, ncol=nspecies)
  diag(U)=1
  # Cholesky decomposition
  A=chol(U)
  community=matrix(rnorm(ntimes*nspecies, sd=1, mean=0), nrow=ntimes, ncol=nspecies) %*%A
  community=scale(community, center=TRUE, scale=TRUE)*sigma+mu
  attr(community, "scaled:center")=NULL
  attr(community, "scaled:scale")=NULL
  
  results=list(rho=rho, sigma=sigma, mu=mu, community=community)
  return (results)  
}
