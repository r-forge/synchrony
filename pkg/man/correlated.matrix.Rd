\name{correlated.matrix}
\alias{correlated.matrix}
\title{
correlated.matrix
}
\description{
Create an \code{ntimes} x \code{nspecies} matrix with correlation \eqn{\rho},
standard deviation \eqn{\sigma}, and mean \eqn{\mu}
}
\usage{
correlated.matrix (rho = 0, sigma = 1, mu = 0, ntimes = 200, nspecies = 10)
}
\arguments{
  \item{rho}{
Correlation between the columns of the matrix. Default is \code{0}
}
  \item{sigma}{
Standard deviation of the columns. Default is 1
}
  \item{mu}{
Mean of the columns. Default is 0
}

\item{ntimes}{
Number of rows in the matrix. Default is 200
}
\item{nspecies}{
Number of columns in the matrix. Default is 10
}
}

\value{
Returns a named list containing the following:

\item{rho }{Correlation between the columns}
\item{sigma }{Standard deviation of the columns}
\item{mu }{Mean of the columns}
\item{community }{\code{ntimes} x \code{nspecies} matrix}

}
\references{
Gouhier, T. C., F. Guichard, and A. Gonzalez. 2010. Synchrony and stability of 
food webs in metacommunities. \emph{The American Naturalist} 175:E16-E34.
}

\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
mat=correlated.matrix(rho=0.85, sigma=30, mu=10, nspecies=10)
## Check sd of each column
apply(mat$community, 2, sd)
## Check mean of each column
apply(mat$community, 2, mean)
## Check correlation of matrix
community.sync(mat$community)
}
