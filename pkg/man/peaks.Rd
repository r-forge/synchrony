\name{peaks}
\alias{peaks}
\title{
Find the proportion of local minima/maxima common to both time series and compute 
its significance via Monte Carlo randomizations
}
\description{
Find the proportion of local minima/maxima common to both time series and compute 
its significance via Monte Carlo randomizations
}
\usage{
peaks (t1, t2, nrands = 0)
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
number of randomizations. Default is \code{0}.
}
}
\value{
Returns a named list containing:

\item{pval }{p-value computed by randomly shuffling both time series \code{nrands} times}
\item{rands }{proportion of local minima/maxima common to both time series for each randomization}
\item{peaks }{proportion of local minima/maxima common to both time series in the observed dataset}
\item{index }{indices of local minima/maxima common to both time series in the observed dataset}
}
\references{
Buonaccorsi, J. P., J. S. Elkinton, S. R. Evans, and A. M. Liebhold. 2001. 
Measuring and testing for spatial synchrony. \emph{Ecology} 82:1668-1679.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
t1=cbind(1:100, runif(100))
t2=cbind(1:100, runif(100))
(p=peaks(t1, t2))

}
