\name{phase.sync}
\alias{phase.sync}
\title{
Phase synchrony of quasi-periodic time series
}
\description{
Compute the phase synchrony between two quasi-periodic time series by quantifying
their phase difference at each time step
}
\usage{
phase.sync (t1, t2,  nrands =  0, type = 1, nbreaks = 10, mins = FALSE, quiet = FALSE)
}
\arguments{
 \item{t1}{
time series 1 in matrix format (\code{n} rows x 2 columns). 
The first column should contain the time 
steps and the second column should contain the values.
}
  \item{t2}{
time series 2 in matrix format (\code{n} rows x 2 columns). 
The first column should contain the time 
steps and the second column should contain the values.
}
  \item{nrands}{
number of randomizations to perform (default is 0)  
}
  \item{type}{flag to indicate whether to compute phase difference modulus \eqn{2\pi} 
  between 0 and \eqn{2\pi} (\code{type=1}) or phase difference modulus 
  \eqn{2\pi} between \eqn{-\pi} and \eqn{\pi}. Default is \code{1}.
  } 
  \item{nbreaks}{number of bins to use to group the values in the time series.
  Default is \code{10}.
  }
  \item{quiet}{
Suppress progress bar when set to \code{TRUE}. Default is \code{FALSE}
}
  \item{mins}{
use local minima instead of local maxima to compute and the interpolate the phase. Default is
\code{FALSE}.
}
}
\details{
Two time series are phase-locked if the relationship between their phases
remains constant over time. This function computes the phase of successive local maxima
or minima for each time series and then uses linear interpolation to find the phase at
time steps that fall between local maxima/minima. A histogram can be used to determine if the
distribution of the phase difference at each time step is uniform (indicating no phase locking) 
or has a clear peak (indicating phase locking).
}
\value{
Returns a list containing \code{phases1}, \code{phases2}, and \code{deltaphase}

\item{phases1 }{\code{n} x 3 matrix containing the timestep, value, and phase of the 
first time series }
\item{phases2 }{\code{n} x 3 matrix containing the timestep, value, and phase of the 
second time series }
\item{deltaphase }{\code{n} x 4 matrix containing the timestep, raw phase difference, 
phase difference modulus \eqn{2\pi} between 0 and \eqn{2\pi}, phase difference 
modulus \eqn{2\pi} between \eqn{-\pi} and \eqn{\pi}}
}
\references{
Cazelles, B., and L. Stone. 2003. Detection of imperfect population synchrony 
in an uncertain world. \emph{Journal of Animal Ecology} 72:953-968.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\examples{
t1=cbind(1:100, runif(100))
t2=cbind(1:100, runif(100))
## Compute and interpolate phases using successive local minima
sync.mins=phase.sync(t1, t2, mins=TRUE)
## Compute and interpolate phases using successive local maxima
sync.maxs=phase.sync(t1, t2)
## Plot distribution of phase difference
hist(sync.mins$deltaphase$mod_phase_diff_2pi)
}
