#' @name hf_data
#' 
#' @title High frequency data sets
#' 
#' @description "hf_data" is a dataset containing
#' 
#' "sym_bol" is a \code{string} containing the name "SPY".
#' "SPY" is an \code{xts} time series containing \code{OHLC} minute bar data for
#' the SPY etf, from 2008-01-02 to 2014-05-19.
#' 
#' an \code{xts} time series containing 625,425 rows of data, each row contains a single minute bar:
#'
#' @docType data
#' @keywords datasets
#' @usage data(hf_data)
#'
#' @format an \code{xts} time series 625425 rows of data, each row contains a single minute bar:
#' \describe{
#'   \item{Open}{Open price in the bar}
#'   \item{High}{High price in the bar}
#'   \item{Low}{Low price in the bar}
#'   \item{Close}{Close price in the bar}
#'   \item{Volume}{trading volume in the bar}
#' }
#' @source \url{http://www.diamondse.info/}
#' 
#' @references Moore et al. (2013) Genetics 195:1077-1086
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/23979570}{PubMed})
#'
#' @source \href{http://qtlarchive.org/db/q?pg=projdetails&proj=moore_2013b}{QTL Archive}
#'
#' @examples
#' # data(hf_data)  # not needed - data is lazy load
#' head(SPY)
#' \donttest{chart_Series(x=SPY["2009"])}
"hf_data"



