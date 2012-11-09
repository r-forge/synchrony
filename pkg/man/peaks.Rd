\name{peaks}
\alias{peaks}
\title{
Find common peaks
}
\description{
Find the proportion of local minima/maxima common to both time series
}
\usage{
peaks (t1, t2)
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
}
\value{
Returns a named list containing \code{peaks} and \code{locations}

\item{peaks }{proportion of local minima/maxima common to both time series}
\item{locations }{locations/indices of local minima/maxima common to both time series}
}
\references{
Buonaccorsi, J. P., J. S. Elkinton, S. R. Evans, and A. M. Liebhold. 2001. 
Measuring and testing for spatial synchrony. \emph{Ecology} 82:1668-1679.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\note{
This is the low-level function used by \code{\link{common.peaks.signif}}
}

\seealso{
\code{\link{common.peaks.signif}}
}
\examples{
t1=cbind(1:100, runif(100))
t2=cbind(1:100, runif(100))
p=peaks(t1, t2)
p$peaks
}
