\name{vario.signif}
\alias{vario.signif}
\title{
vario.signif
}
\description{
Determine the significance of the empirical variogram via Monte Carlo randomizations
}
\usage{vario.signif (nrands = 999, bins, all.combs, grpdata, data, data2, type, 
                              vario, glob.mean = NULL, glob.sd = NULL, 
                              glob.N = NULL, is.multivar = FALSE,
                              alternative = c("one.tailed", "two.tailed"),
                              centered, regional.mean,
                              mult.test.corr=c(FALSE, "holm", "hochberg", 
                                               "sidak", "bonferroni"))
                              }
                                                 
\arguments{
  \item{nrands}{Number of randomizations}
  \item{bins}{Number of bins or lag distances}
  \item{all.combs}{All combinations of pairs of sites}  
  \item{grpdata}{Grouped sites based on distance bins}
  \item{data}{\code{n} x \code{m} matrix containing y-coordinates (or latitude), 
              x-coordinates (or longitude), and values. The values can either be
              a single column of observations at each site for univariate variograms
              or a matrix of observations at each site for multivariate variograms 
              (e.g., to compute spatial synchrony).
              }
  \item{data2}{\code{n} x \code{m} matrix containing y-coordinates (or latitude), 
              x-coordinates (or longitude), and values. The values can either be
              a single column of observations at each site for univariate variograms
              or a matrix of observations at each site for multivariate variograms 
              (e.g., to compute spatial synchrony).
              }
  \item{type}{Type of variogram to compute. Default is \code{semivar} for semivariance. 
  Other options include \code{cov} for covariance, \code{pearson} for Pearson
correlation, \code{spearman} for Spearman correlation, \code{kendall} for Kendall correlation,
\code{moran} for Moran's I, and \code{geary} for Geary's C}
\item{vario}{Empirical variogram}
\item{glob.mean}{Global mean}
\item{glob.sd}{Global standard deviation}
\item{glob.N}{Global number of points}
\item{is.multivar}{Is the data multivariate? Default is \code{FALSE}}
\item{alternative}{Conduct a one-tailed or a two-tailed test? Default is \code{one.tailed}}
\item{centered}{Is the variogram centered?}
\item{regional.mean}{Regional mean}
\item{mult.test.corr}{Correct for multiple tests? Default is \code{FALSE}. Other options inlcude
\code{holm}, \code{hochberg}, \code{sidak}, and \code{bonferroni}
}
}
\value{
Returns a named list containing the following variables:

\item{rands }{\code{nrands} x \code{nbins} matrix of randomizations}
\item{pvals }{p-values based on one-tailed test}
\item{vario }{Variogram values in each lag/bin}
\item{mult.test.corr }{Corrected for multiple tests?}
\item{alternative }{One-tailed or two-tailed test?}
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\seealso{
\code{\link{vario}, \link{vario.func}}
}
\examples{
## Internal function used by vario
}
\keyword{ variogram }
