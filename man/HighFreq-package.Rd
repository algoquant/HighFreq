\name{HighFreq-package}
\alias{HighFreq-package}
\alias{HighFreq}
\docType{package}
\title{
HighFreq R package
}
\description{
This package contains functions for managing high frequency TAQ and OHLC market data.
The functions perform:
\itemize{
  \item reading and writing data, 
  \item managing time zones and alligning time indices, 
  \item chaining and joining time series, 
  \item scrubbing bad data points, 
  \item converting TAQ data to OHLC format, 
  \item aggregating data to lower frequency.
}
This package is inspired by the package 'highfrequency', and follows many of its formats and standards.
For example, the high frequency data is assumed to be stored in separate directories for each symbol, containing *.RData files, each file with one day of data.
The functions are designed for speed and reduced memory usage.  Wherever possible they use vectorized functions from packages 'base', 'xts', 'quantmod', and 'caTools'.
}
\details{
\tabular{ll}{
Package: \tab HighFreq\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-02-08\cr
License: \tab GPL (>= 2)\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Jerzy Pawlowski
Maintainer: Jerzy Pawlowski <algoquant@algoquants.ch>
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
