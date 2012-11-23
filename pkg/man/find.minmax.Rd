\name{find.minmax}
\alias{find.minmax}
\title{
Find min/max of a time series
}
\description{
Find local minima and maxima of a time series
}
\usage{
find.minmax (timeseries)
}
\arguments{
  \item{timeseries}{
time series in matrix format (\code{n} rows x 2 columns). The first column should contain the time 
steps and the second column should contain the values.
}
}

\value{
Returns a named list containing:

\item{mins }{\code{n} x 3 matrix containing the index, time steps, and the local min values}
\item{maxs }{\code{n} x 3 matrix containing the index, time steps, and the local max values}
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
t1=cbind(1:100, runif(100))
min.max=find.minmax(t1)
min.max$maxs
plot (t1, t="l")
points (min.max$mins, col="blue", bg="blue", pch=19)
points (min.max$maxs, col="red", bg="red", pch=19)
}
\keyword{ max }
\keyword{ min }
