\name{peaks.signif}
\alias{peaks.signif}
\title{
Find proportion of common peaks/troughs and compute its significance
}
\description{
Calculates the proportion of local minima/maxima common to both time series and computes its
significance by randomly shuffling the time series
}
\usage{
peaks.signif (t1, t2, nrands = 999)
}
\arguments{
  \item{t1}{time series 1 in matrix format (\code{n} rows x 2 columns). The first column 
  should contain the time steps and the second column should contain the values.
}
  \item{t2}{time series 2 in matrix format (\code{n} rows x 2 columns). The first column 
  should contain the time steps and the second column should contain the values.
}
  \item{nrands}{
number of randomizations. Default is \code{999}.
}
}

\value{
Returns a named list containing:

\item{pval }{p-value computed by randomly shuffling both time series \code{nrands} times}
\item{rands }{proportion of local minima/maxima common to both time series for each randomization}
\item{peaks }{proportion of local minima/maxima common to both time series in the observed dataset}
\item{locations }{locations/indices of local minima/maxima common to both time series in the observed dataset}
}
\references{
Buonaccorsi, J. P., J. S. Elkinton, S. R. Evans, and A. M. Liebhold. 2001. 
Measuring and testing for spatial synchrony. \emph{Ecology} 82:1668-1679.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\seealso{
\code{\link{peaks}} and \code{\link{find.minmax}}
}

\examples{
t1=cbind(1:100, runif(100))
t2=cbind(1:100, runif(100))
p=peaks.signif(t1, t2)
## Find proportion of time steps where both time series peak together
p$peaks
## Plot (null) distribution of proportion of time steps where both time series
## peak together
hist(p$rand)
## p-value of observed value
p$pval
}
