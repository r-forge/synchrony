\name{phase.sync.signif}
\alias{phase.sync.signif}
\title{Compute significance of phase synchrony}
\description{
Compute the significance of phase synchrony between two time series via Monte Carlo randomizations
}
\usage{
phase.sync.signif (t1, t2, nrands =  999, type = 1, nbreaks = 10, mins = FALSE)
}
\arguments{
  \item{t1}{
time series 1 in matrix format (\code{n} rows x 2 columns). The first column should contain the time 
steps and the second column should contain the values.
}
  \item{t2}{
time series 2 in matrix format (\code{n} rows x 2 columns). The first column should contain the time 
steps and the second column should contain the values.
}
  \item{nrands}{
number of randomizations to perform (default is 999)  
}
  \item{type}{flag to indicate whether to compute phase difference modulus \eqn{2\pi} 
  between 0 and \eqn{2\pi} (\code{type=1}) or phase difference modulus \eqn{2\pi} between \eqn{-\pi} and \eqn{\pi}. Default is \code{1}.
  } 
  \item{nbreaks}{number of bins to use to group the values in the time series.
  Default is \code{10}.
  } 
    \item{mins}{
use local minima instead of local maxima to compute and the interpolate the phase. Default is
\code{FALSE}.
}
}

\details{
This function determines whether two time series are phase locked by randomly shuffling the time series. Specifically, the distribution of the phase difference of the 
original time series is binned using a histogram and an index of phase locking \eqn{Q} based on Shannon's entropy is computed: \eqn{Q.obs=(S_{max}-S_{obs})/S_{max}}, where
\eqn{S_{obs}=-\sum\limits_{i}^N{p_{i}log(p_{i})}}, \eqn{p_{i}} represents the proportion of phase
differences in each histogram bin, \eqn{N} the total number of bins, and \eqn{S_{max}=log(N)}.
\eqn{Q.obs=0} when there is no phase synchrony and \eqn{Q.obs=1} when there is perfect phase synchrony.
The statistical significance of the phase synchrony \eqn{Q.obs} can then be computed by
comparing the observed \eqn{Q.obs} value to the \eqn{Q} values obtained from \code{nrands} surrogate
time series whose short-term correlations and overall temporal patterns match those of the observed
time series. The p-value is then the proportion of surrogate time series whose \eqn{Q \ge Q.obs}.
}

\value{
Returns a named list containing:

\item{Q.obs }{Observed Q value}
\item{pval }{p-value of observed \code{Q.obs}}
\item{rands }{Q values from all randomizations}
\item{phases1 }{\code{n} x 3 matrix containing the timestep, value, and phase of the 
first time series }
\item{phases2 }{\code{n} x 3 matrix containing the timestep, value, and phase of the 
second time series }
\item{deltaphase }{\code{n} x 4 matrix containing the timestep, raw phase difference, 
phase difference modulus \eqn{2\pi} between 0 and \eqn{2\pi}, phase difference 
modulus \eqn{2\pi} between \eqn{-\pi} and \eqn{\pi}}
\item{icdf }{\code{nrands} x 2 matrix containing the inverse cumulative distribution
of the Q values computed using the surrogate time series}
}

\references{
Cazelles, B., and L. Stone. 2003. Detection of imperfect population synchrony 
in an uncertain world. \emph{Journal of Animal Ecology} 72:953-968.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\seealso{
\code{\link{phase.sync}} and \code{\link{surrogate.ts}}
}

\examples{
t1=cbind(1:100, runif(100))
t2=cbind(1:100, runif(100))
p=phase.sync.signif(t1, t2, nrands=99)
p$pval
hist(p$rands)
plot (p$icdf, t="l")
abline(v=p$Q.obs, col="red", lty=2)
}
