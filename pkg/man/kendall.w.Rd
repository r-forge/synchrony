\name{kendall.w}
\alias{kendall.w}
\alias{concordance}
\title{
Kendall's W
}
\description{
Compute Kendall's coefficient of concordance (W)
}
\usage{
kendall.w (data)
}
\arguments{
  \item{data}{
  matrix in wide format where each row represents a different sample and each column 
  represents a different variable.
}
}
\details{
Kendall's W is a non-parametric statistic that ranges from 0 to 1 and measures 
the level of agreement between multiple variables. When the number of observations \eqn{n>10}, 
its significance can be determined by using a \eqn{\chi^2} distribution with \eqn{df=n-1}.
Kendall's W is often used to measure the degree of synchrony among different locations.
}
\value{
Returns a named list containing:

\item{w.uncorrected }{Kendall's W uncorrected for tied ranks}
\item{w.corrected }{Kendall's W corrected for tied ranks}
\item{pval }{p-value of Kendall's W}
\item{spearman.corr }{Spearman's ranked correlation}
}
\references{
Buonaccorsi, J. P., J. S. Elkinton, S. R. Evans, and A. M. Liebhold. 2001. 
Measuring and testing for spatial synchrony. \emph{Ecology} 82:1668-1679.

Gouhier, T. C., and F. Guichard. 2007. Local disturbance cycles and the maintenance 
of spatial heterogeneity across scales in marine metapopulations. \emph{Ecology} 88:647-657.

Gouhier, T. C., F. Guichard, and A. Gonzalez. 2010. Synchrony and stability of 
food webs in metacommunities. \emph{The American Naturalist} 175:E16-E34.

Zar, J. H. 1999. Biostatistical Analysis, Fourth edition. Prentice-Hall, Inc., 
Upper Saddle River, NJ.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
data(bird.traits)
(w=kendall.w(bird.traits))
}
