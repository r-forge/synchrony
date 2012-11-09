\name{community.sync}
\alias{community.sync}
\title{
Compute community-wide synchrony
}
\description{
Compute community-wide synchrony. If all species fluctuate in perfect unison, the community-wide
synchrony will be 1. If species undergo uncorrelated fluctuations, the community-wide synchrony 
will be 0.
}
\usage{
community.sync (comm.matrix)
}
\arguments{
 \item{comm.matrix}{
community matrix in wide format where each row contains the abundance at each time step
and each column corresponds to a different species.
}
}
\details{
Loreau and de Mazancourt (2008) show that community-wide synchrony \eqn{\varphi} can be quantified
by computing the temporal variance \eqn{\sigma_{x_T}^2} of the community time series
\eqn{x_T(t)=\sum{x_i(t)}} and the sum of the temporal standard deviation of the time series
across all species \eqn{\left(\sum{\sigma_{x_i}}\right)^2} such that:
\eqn{\varphi=\frac{\sigma_{x_T}^2}{\left(\sum{\sigma_{x_i}}\right)^2}}
}

\value{
Returns community-wide synchrony \eqn{\varphi}
}
\references{
Loreau, M., and C. de Mazancourt. 2008. Species synchrony and its drivers: 
Neutral and nonneutral community dynamics in fluctuating environments. 
\emph{The American Naturalist} 172:E48-E66.
}

\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
## Community matrix for 20 species undergoing random fluctuations 
comm.rand=matrix(runif(100), nrow=5, ncol=20)
community.sync(comm.rand)
## Community matrix for 20 species undergoing synchronized fluctuations 
comm.corr=matrix(rep(comm.rand[,1], 20), nrow=5, ncol=20)
community.sync(comm.corr)
## Compare this metric to Kendall's W
data(bird.traits)
kendall.w(bird.traits)
community.sync(bird.traits)
}
