\name{phase.partnered}
\alias{phase.partnered}
\title{
Phase partnered time series
}
\description{
Create two time series with specific autocorrelation \eqn{\gamma} and cross-correlation
\eqn{\rho}
}
\usage{
phase.partnered (n = 2000, rho = 1, gamma = 1, ts.sd = 0.1, ts.mean = 0)
}
\arguments{
  \item{n}{
number of time steps in time series. Default is \code{2000}.
}
  \item{rho}{
cross-correlation between the two time series (\eqn{-1\le \rho \le 1}).
Default is \code{1}.
}
  \item{gamma}{
autocorrelation of each time series. Gamma (\eqn{\gamma}) describes the relationship between 
frequency \eqn{f} and power \eqn{P}: \eqn{P(f)=1/f^\gamma}


\eqn{-2\le \gamma \le 0}: blue noise and \eqn{0\le \gamma \le 2}: red noise. Default is \code{1}.
}
  \item{ts.sd}{
standard deviation of both time series. Default is \code{0.1}.
}
  \item{ts.mean}{
mean of both time series. Default is \code{0}.
}
}

\value{
Returns an \code{n} x 2 matrix containing the time series
}
\references{
Gouhier, T. C., F. Guichard, and A. Gonzalez. 2010. Synchrony and stability of 
food webs in metacommunities. \emph{The American Naturalist} 175:E16-E34.

Vasseur, D. A. 2007. Environmental colour intensifies the Moran effect when 
population dynamics are spatially heterogeneous. \emph{Oikos} 116:1726-1736.}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
## Positively cross-correlated white noise
pos.corr=phase.partnered(n = 100, rho = 0.7, gamma = 0)
## Negatively cross-correlated white noise
neg.corr=phase.partnered(n = 100, rho = -1, gamma = 0)
par(mfrow=c(2,1))
matplot (pos.corr, t="l", lty=1)
matplot (neg.corr, t="l", lty=1)
}
