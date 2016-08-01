#' @name hf_data
#' @docType data
#' @keywords datasets
#'
#' @title High frequency data sets
#'
#' @description hf_data.RData is a file containing the datasets:
#'
#' \describe{
#'   \item{SPY}{an \code{xts} time series containing 1-minute \code{OHLC} bar
#'   data for the \code{SPY} etf, from 2008-01-02 to 2014-05-19.  SPY contains 625,425 
#'   rows of data, each row contains a single minute bar.}
#'   \item{TLT}{an \code{xts} time series containing 1-minute \code{OHLC} bar
#'   data for the \code{TLT} etf, up to 2014-05-19.}
#'   \item{VXX}{an \code{xts} time series containing 1-minute \code{OHLC} bar
#'   data for the \code{VXX} etf, up to 2014-05-19.}
#' }
#'
#' @format Each \code{xts} time series contains \code{OHLC} data, with each row 
#'   containing a single minute bar:
#' \describe{
#'   \item{Open}{Open price in the bar}
#'   \item{High}{High price in the bar}
#'   \item{Low}{Low price in the bar}
#'   \item{Close}{Close price in the bar}
#'   \item{Volume}{trading volume in the bar}
#' }
#' @source \url{https://wrds-web.wharton.upenn.edu/wrds/}
#'
#' @references Wharton Research Data Service
#' (\href{https://wrds-web.wharton.upenn.edu/wrds/}{WRDS})
#'
#' @usage data(hf_data)  # not required - data is lazy load
#'
#' @examples
#' # data(hf_data)  # not required - data is lazy load
#' head(SPY)
#' \donttest{chart_Series(x=SPY["2009"])}
"SPY"
"TLT"
"VXX"

