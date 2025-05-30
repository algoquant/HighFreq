% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{roll_var}
\alias{roll_var}
\title{Calculate a \emph{matrix} of dispersion (variance) estimates over a rolling
look-back interval attached at the end points of a \emph{time series} or a
\emph{matrix}.}
\usage{
roll_var(
  timeser,
  lookb = 1L,
  startp = 0L,
  endd = 0L,
  step = 1L,
  stub = 0L,
  method = "moment",
  confl = 0.75
)
}
\arguments{
\item{\code{timeser}}{A \emph{time series} or a \emph{matrix} of data.}

\item{\code{lookb}}{The number of end points in the look-back interval
(the default is \code{lookb = 1}).}

\item{\code{startp}}{An \emph{integer} vector of start points (the default
is \code{startp = 0}).}

\item{\code{endd}}{An \emph{integer} vector of end points (the default is
\code{endd = 0}).}

\item{\code{step}}{The number of time periods between the end points (the
default is \code{step = 1}).}

\item{\code{stub}}{An \emph{integer} value equal to the first end point for
calculating the end points (the default is \code{stub = 0}).}

\item{\code{method}}{A \emph{character string} representing the type of the
measure of dispersion (the default is \code{method = "moment"}).}
}
\value{
A \emph{matrix} dispersion (variance) estimates with the same number
  of columns as the input time series \code{timeser}, and the number of rows
  equal to the number of end points.
}
\description{
Calculate a \emph{matrix} of dispersion (variance) estimates over a rolling
look-back interval attached at the end points of a \emph{time series} or a
\emph{matrix}.
}
\details{
The function \code{roll_var()} calculates a \emph{matrix} of dispersion
  (variance) estimates over rolling look-back intervals attached at the end
  points of the \emph{time series} \code{timeser}.
  
  The function \code{roll_var()} performs a loop over the end points, and at
  each end point it subsets the time series \code{timeser} over a look-back
  interval equal to \code{lookb} number of end points.
  
  It passes the subset time series to the function \code{calc_var()}, which
  calculates the dispersion.
  See the function \code{calc_var()} for a description of the dispersion
  methods.
  
  If the arguments \code{endd} and \code{startp} are not given then it
  first calculates a vector of end points separated by \code{step} time
  periods. It calculates the end points along the rows of \code{timeser}
  using the function \code{calc_endpoints()}, with the number of time
  periods between the end points equal to \code{step} time periods.

  For example, the rolling variance at \code{25} day end points, with a
  \code{75} day look-back, can be calculated using the parameters
  \code{step = 25} and \code{lookb = 3}.

  The function \code{roll_var()} with the parameter \code{step = 1}
  performs the same calculation as the function \code{roll_var()} from
  package
  \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
  but it's several times faster because it uses \code{RcppArmadillo} \code{C++}
  code.

  The function \code{roll_var()} is implemented in \code{RcppArmadillo}
  \code{RcppArmadillo} \code{C++} code, which makes it several times faster
  than \code{R} code.
}
\examples{
\dontrun{
# Define time series of returns using package rutils
retp <- na.omit(rutils::etfenv$returns$VTI)
# Calculate the rolling variance at 25 day end points, with a 75 day look-back
varv <- HighFreq::roll_var(retp, lookb=3, step=25)
# Compare the variance estimates over 11-period look-back intervals
all.equal(HighFreq::roll_var(retp, lookb=11)[-(1:10), ], 
  drop(RcppRoll::roll_var(retp, n=11)), check.attributes=FALSE)
# Compare the speed of HighFreq::roll_var() with RcppRoll::roll_var()
library(microbenchmark)
summary(microbenchmark(
  Rcpp=HighFreq::roll_var(retp, lookb=11),
  RcppRoll=RcppRoll::roll_var(retp, n=11),
  times=10))[, c(1, 4, 5)]  # end microbenchmark summary
# Compare the speed of HighFreq::roll_var() with TTR::runMAD()
summary(microbenchmark(
    Rcpp=HighFreq::roll_var(retp, lookb=11, method="quantile"),
    TTR=TTR::runMAD(retp, n = 11),
    times=10))[, c(1, 4, 5)]  # end microbenchmark summary
}  # end dontrun
}
