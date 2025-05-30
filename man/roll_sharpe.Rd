% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HighFreq.R
\name{roll_sharpe}
\alias{roll_sharpe}
\title{Calculate a time series of Sharpe ratios over a rolling look-back interval
for an \emph{OHLC} time series.}
\usage{
roll_sharpe(ohlc, lookb = 11)
}
\arguments{
\item{\code{ohlc}}{An \emph{OHLC} time series of prices in \emph{xts} format.}

\item{\code{lookb}}{The size of the look-back interval, equal to the number of
rows of data used for aggregating the \emph{OHLC} prices.}
}
\value{
An \emph{xts} time series with a single column and the same number of
  rows as the argument \code{ohlc}.
}
\description{
Calculate a time series of Sharpe ratios over a rolling look-back interval
for an \emph{OHLC} time series.
}
\details{
The function \code{roll_sharpe()} calculates the rolling Sharpe 
  ratio defined as the ratio of percentage returns over the look-back
  interval, divided by the average volatility of percentage returns.
}
\examples{
# Calculate rolling Sharpe ratio over SPY
sharper <- roll_sharpe(ohlc=HighFreq::SPY, lookb=11)

}
