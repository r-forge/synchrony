\name{find.minmax.ind}
\alias{find.minmax.ind}
\title{
Find indices of local minima or maxima of a time series
}
\description{
Find indices of local minima or maxima of a time series
}
\usage{
find.minmax.ind (timeseries, mins = FALSE)
}
\arguments{
 \item{timeseries}{
time series in matrix format (\code{n} rows x 2 columns). The first column should contain the time 
steps and the second column should contain the values.
}
 \item{mins}{
find local minima (\code{mins=TRUE}) or maxima (\code{mins=FALSE}). Default is \code{FALSE}.}
}
\value{
Returns index of local minima or maxima
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\note{
This is the low-level function used by \code{\link{find.minmax}}. The optimized 
version of this function is based on code posted by Ben Bolker on StackOverflow.
}

\seealso{
\code{\link{find.minmax}}
}

\examples{
t1=cbind(1:100, runif(100))
## Find local maxima
maxs=find.minmax.ind(t1)
## Find local minima
mins=find.minmax.ind(t1, mins=TRUE)
}
\keyword{ max }
\keyword{ min }
