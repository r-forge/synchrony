\name{meancorr}
\alias{meancorr}
\title{
Compute mean column-wise correlation and determine its significance via Monte Carlo randomizations
}
\description{
Compute mean column-wise correlation and determine its significance via Monte Carlo randomizations. 
The Monte Carlo randomizations are performed by shuffling the columns of the community
matrix independently.
}
\usage{
meancorr (comm.matrix, nrands = 0, alternative=c("one.tailed", "two.tailed"), 
                            method=c("pearson", "kendall", "spearman"), 
                            quiet = FALSE, ...)
}
\arguments{
 \item{comm.matrix}{
community matrix in wide format where each row contains the abundance at each time step
and each column corresponds to a different species.
}
 \item{nrands}{
number of randomizations to perform (default is 0)
}
\item{alternative}{Conduct a one-tailed or a two-tailed test? Default is \code{one.tailed}.}
\item{method}{Method to compute correlation? Options include \code{pearson}, \code{kendall}, 
and \code{spearman}. Default is \code{pearson}}
\item{quiet}{
Suppress progress bar when set to \code{TRUE}. Default is \code{FALSE}
}
\item{\dots}{
Other parameters to \code{\link{cor}} function.
}
}

\value{
Returns a named list containing:

\item{obs }{the observed mean correlation}
\item{rands }{the mean correlation for each randomization. 
This variable is only returned if \code{nrands > 0}.}
\item{pval }{p-value of observed mean correlation. 
This variable is only returned if \code{nrands > 0}.}
\item{alternative }{One-tailed or two-tailed test?
This variable is only returned if \code{nrands > 0}.}
\item{method }{Method used to compute the mean correlation.}
}

\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
## Community matrix for 20 species undergoing random fluctuations 
comm.rand=matrix(runif(100), nrow=5, ncol=20)
meancorr(comm.rand, nrands=20)$pval
## Community matrix for 20 species undergoing synchronized fluctuations 
comm.corr=matrix(rep(comm.rand[,1], 20), nrow=5, ncol=20)
meancorr(comm.corr, nrands=20)$pval
## On "real" data
data(bird.traits)
meancorr(bird.traits, nrands=20)$pval
}
