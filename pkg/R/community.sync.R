## Community matrix comm.matrix: n x m matrix with n=time step, m=species
community.sync <- function (comm.matrix) {
  species.sd=apply(comm.matrix, MARGIN=2, FUN=sd)
  community.var=var(rowSums(comm.matrix))
  sync=community.var/sum(species.sd, na.rm=TRUE)^2
  return (sync)
}