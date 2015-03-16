\name{HighFreq-package}
\alias{HighFreq-package}
\alias{HighFreq}
\docType{package}
\title{
HighFreq R package
}
\description{
Package \sQuote{\code{HighFreq}} contains functions for managing high frequency \code{TAQ} and \code{OHLC} market data.
The functions perform:
\itemize{
  \item reading and writing data from files, 
  \item managing time zones and alligning time indices, 
  \item chaining and joining time series, 
  \item scrubbing bad data points, 
  \item converting \code{TAQ} data to \code{OHLC} format, 
  \item aggregating data to lower frequency.
}
Package \sQuote{\code{HighFreq}} is inspired by the package \sQuote{\href{http://cran.r-project.org/web/packages/highfrequency/index.html}{highfrequency}}, and follows many of its conventions.  For example, the high frequency data is assumed to be stored in separate directories for each symbol, containing multiple daily \sQuote{\code{*.RData}} files, each file containing one day of \code{TAQ} data.  

Package \sQuote{\code{HighFreq}} also follows many of the conventions of package \sQuote{\href{http://cran.r-project.org/web/packages/quantmod/index.html}{quantmod}}, such column names of \code{OHLC} objects.  

The \sQuote{\code{HighFreq}} functions are designed for speed and reduced memory usage.  Wherever possible they use vectorized functions from packages \code{base}, \sQuote{\href{http://cran.r-project.org/web/packages/xts/index.html}{xts}}, \sQuote{\href{http://cran.r-project.org/web/packages/quantmod/index.html}{quantmod}}, and \sQuote{\href{http://cran.r-project.org/web/packages/caTools/index.html}{caTools}}.  
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
