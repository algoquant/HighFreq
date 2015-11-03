#' Identify extreme values in a univariate \code{xts} time series.
#' 
#' Identifies extreme values as those that exceed a multiple of the running volatility.
#' 
#' @param time_series univariate \code{xts} time series.
#' @param vol_window number of data points for estimating running volatility.
#' @param vol_mult volatility multiplier.
#' @return  \code{logical vector}.
#' @details Calculates running volatility as a quantile of values over a sliding
#'   window. Extreme values are those that exceed the product of the volatility 
#'   multiplier times the running volatility. Extreme values are the very tips
#'   of the tails when the distribution of values becomes very fat-tailed. The
#'   volatility multiplier \code{vol_mult} controls the threshold at which
#'   values are identified as extreme. Smaller volatility multiplier values will
#'   cause more values to be identified as extreme.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' # scrub extreme values
#' x_ts <- x_ts[!extreme_values(x_ts, vol_mult=1)]
extreme_values <- function(time_series, vol_window=51, vol_mult=2) {

# calculate volatility as running quantile
  vo_lat <- runquantile(x=abs(as.vector(time_series)), k=vol_window, 
                        probs=0.9, endrule="constant", align="center")
  vo_lat <- xts(vo_lat, order.by=index(time_series))
  colnames(vo_lat) <- "volat"
# carry forward non-zero volatility values
  vo_lat[vo_lat==0] <- NA
  vo_lat[1] <- 1
  vo_lat <- na.locf(vo_lat)
#  vo_lat <- na.locf(vo_lat, fromLast=TRUE)

# extreme value if time_series greater than scaled volatility
  ex_treme <- (abs(time_series) > 2*vol_mult*vo_lat)
  ex_treme[1] <- FALSE
  colnames(ex_treme) <- "suspect"

  cat("date:", format(as.Date(index(first(time_series)))), "\tfound", sum(ex_treme), "extreme values\n")
  ex_treme
}  # end extreme_values




#' Identify isolated price jumps in a univariate \code{xts} time series of prices,
#' based on pairs of large neighboring returns of opposite sign.
#' 
#' @param time_series univariate \code{xts} time series of prices.
#' @inheritParams extreme_values
#' @return  \code{logical vector}.
#' @details Isolated price jumps are single prices that are very different from 
#'   neighboring values.  Price jumps create pairs of large neighboring returns
#'   of opposite sign. The function \code{price_jumps} first calculates simple
#'   returns from prices. Then it calculates the running volatility of returns
#'   as a quantile of returns over a sliding window. Jump prices are identified
#'   as those where neighboring returns both exceed a multiple of the running
#'   volatility, but the sum of those returns doesn't exceed it.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' # scrub jump prices
#' x_ts <- x_ts[!price_jumps(x_ts, vol_mult=1.0)]
price_jumps <- function(time_series, vol_window=51, vol_mult=2) {

# calculate simple returns
  diff_series <- diff(time_series)
  diff_series[1, ] <- 0
  colnames(diff_series) <- "diffs"
  diff_series_fut <- lag(diff_series, -1)
  diff_series_fut[nrow(diff_series_fut)] <- 0
  colnames(diff_series_fut) <- "diff_series_fut"

# calculate vo_lat as running quantile
  vo_lat <- runquantile(x=abs(as.vector(diff_series)), k=vol_window, 
                        probs=0.9, endrule="constant", align="center")
  vo_lat <- xts(vo_lat, order.by=index(diff_series))
  colnames(vo_lat) <- "volat"
# carry forward non-zero vo_lat values
  vo_lat[vo_lat==0] <- NA
  vo_lat[1] <- 1
  vo_lat <- na.locf(vo_lat)
#  vo_lat <- na.locf(vo_lat, fromLast=TRUE)

# suspect value if abs diffs greater than vo_lat, and if abs sum of diffs less than vo_lat
  sus_pect <- (
    (abs(diff_series) > vol_mult*vo_lat) & 
      (abs(diff_series_fut) > vol_mult*vo_lat) & 
      (abs(diff_series+diff_series_fut) < 2*vol_mult*vo_lat)
  )
  sus_pect[1] <- FALSE
  colnames(sus_pect) <- "suspect"
# cat("Parsing", deparse(substitute(taq_data)), "\n")
# cat("Parsing", strsplit(deparse(substitute(taq_data)), split="[.]")[[1]][4], "on date:", format(to_day), "\tfound", sum(sus_pect), "suspect prices\n")
  cat("date:", format(as.Date(index(first(time_series)))), "\tfound", sum(sus_pect), "jump prices\n")
  sus_pect
}  # end price_jumps




#' Scrub a single day of \code{TAQ} data in \code{xts} format, without aggregation.
#' 
#' @inheritParams extreme_values
#' @param taq_data \code{TAQ xts} time series.
#' @param tzone timezone to convert.
#' @return  \code{TAQ xts} time series.
#' @details The function \code{scrub_TAQ} performs the same scrubbing operations
#'   as \code{scrub_agg}, except it doesn't aggregate, and returns the
#'   \code{TAQ} data in \code{xts} format.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"),
#'              to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random TAQ prices
#' x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
#' # create vector of random bid-offer prices
#' bid_offer <- abs(rnorm(length(in_dex)))/10
#' # create TAQ data using cbind
#' taq_data <- cbind(x_ts-bid_offer, x_ts+bid_offer)
#' # add Trade.Price
#' taq_data <- cbind(taq_data, x_ts+rnorm(length(in_dex))/10)
#' # add Volume
#' taq_data <- cbind(taq_data, sample(x=10*(2:18), size=length(in_dex), replace=TRUE))
#' colnames(taq_data) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
#' taq_data <- scrub_TAQ(taq_data)
#' taq_data <- scrub_TAQ(taq_data, vol_window=11, vol_mult=1)
scrub_TAQ <- function(taq_data, vol_window=51, vol_mult=2, tzone="America/New_York") {

# convert timezone of index to New_York
  index(taq_data) <- with_tz(time=index(taq_data), tzone=tzone)
# subset data to NYSE trading hours
  taq_data <- taq_data['T09:30:00/T16:00:00', ]
# return NULL if no data
  if (nrow(taq_data)==0)  return(NULL)
#  to_day <- as.Date(index(first(taq_data)))

# remove duplicate time stamps using 'duplicated'
  taq_data <- taq_data[!duplicated(index(taq_data)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- taq_data[, "Ask.Price"] - taq_data[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- extreme_values(bid_offer, vol_window=vol_window, vol_mult=vol_mult)
# remove suspect values
  taq_data <- taq_data[!sus_pect]
# replace suspect values
# taq_data[sus_pect, "Bid.Price"] <- taq_data[sus_pect, "Trade.Price"]
# taq_data[sus_pect, "Ask.Price"] <- taq_data[sus_pect, "Trade.Price"]

# scrub quotes with suspect price jumps
# calculate mid prices
  mid_prices <- 0.5 * (taq_data[, "Bid.Price"] + taq_data[, "Ask.Price"])
#  mid_prices <- na.omit(mid_prices)
#  colnames(mid_prices) <- "Mid.Price"
# replace NA volumes with zero
  taq_data[is.na(taq_data[, "Volume"]), "Volume"] <- 0
# replace suspect values with NA, and perform 'locf'
  taq_data[price_jumps(mid_prices, vol_window=vol_window, vol_mult=vol_mult), ] <- NA
  na.locf(taq_data)
}  # end scrub_TAQ




#' Scrub a single day of \code{TAQ} data, aggregate it, and convert to \code{OHLC} format.
#' 
#' @inheritParams scrub_TAQ
#' @param period aggregation period.
#' @return  \code{OHLC} time series in \code{xts} format.
#' @details The function \code{scrub_agg} performs:
#' \itemize{
#'   \item index timezone conversion,
#'   \item data subset to trading hours,
#'   \item removal of duplicate time stamps,
#'   \item scrubbing of quotes with suspect bid-offer spreads,
#'   \item scrubbing of quotes with suspect price jumps,
#'   \item cbinding of mid prices with volume data, 
#'   \item aggregation to OHLC using function \code{to.period} from package \code{xts},
#' }
#' Valid 'period' character strings include: "minutes", "3 min", "5 min", "10
#' min", "15 min", "30 min", and "hours". The time index of the output time
#' series is rounded up to the next integer multiple of 'period'.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"), to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random TAQ prices
#' x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
#' # create vector of random bid-offer prices
#' bid_offer <- abs(rnorm(length(in_dex)))/10
#' # create TAQ data using cbind
#' taq_data <- cbind(x_ts-bid_offer, x_ts+bid_offer)
#' # add Trade.Price
#' taq_data <- cbind(taq_data, x_ts+rnorm(length(in_dex))/10)
#' # add Volume
#' taq_data <- cbind(taq_data, sample(x=10*(2:18), size=length(in_dex), replace=TRUE))
#' colnames(taq_data) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
#' # aggregate to ten minutes OHLC data
#' ohlc_data <- scrub_agg(taq_data, period="10 min")
#' chartSeries(ohlc_data, name=sym_bol, theme=chartTheme("white"))
scrub_agg <- function(taq_data, 
                      vol_window=51, 
                      vol_mult=2, 
                      period="minutes", 
                      tzone="America/New_York") {

# convert timezone of index to New_York
  index(taq_data) <- with_tz(time=index(taq_data), tzone=tzone)
# subset data to NYSE trading hours
  taq_data <- taq_data['T09:30:00/T16:00:00', ]
# return NULL if no data
  if (nrow(taq_data)==0)  return(NULL)
#  to_day <- as.Date(index(first(taq_data)))

# remove duplicate time stamps using 'duplicated'
  taq_data <- taq_data[!duplicated(index(taq_data)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- taq_data[, "Ask.Price"] - taq_data[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- extreme_values(bid_offer, vol_window=vol_window, vol_mult=vol_mult)
# remove suspect values
  taq_data <- taq_data[!sus_pect]
# replace suspect values
# taq_data[sus_pect, "Bid.Price"] <- taq_data[sus_pect, "Trade.Price"]
# taq_data[sus_pect, "Ask.Price"] <- taq_data[sus_pect, "Trade.Price"]

# scrub quotes with suspect price jumps
# calculate mid prices
  mid_prices <- 0.5 * (taq_data[, "Bid.Price"] + taq_data[, "Ask.Price"])
#  mid_prices <- na.omit(mid_prices)
  colnames(mid_prices) <- "Mid.Price"
# replace suspect values with NA, and perform 'locf'
  mid_prices[price_jumps(mid_prices, vol_window=vol_window, vol_mult=vol_mult)] <- NA
  mid_prices <- na.locf(mid_prices)
#  mid_prices <- na.locf(mid_prices, fromLast=TRUE)
# cbind mid_prices with volume data, and replace NA volumes with zero
  mid_prices <- cbind(mid_prices, taq_data[index(mid_prices), "Volume"])
  mid_prices[is.na(mid_prices[, "Volume"]), "Volume"] <- 0

# aggregate to OHLC and cumulative volume data
  mid_prices <- switch(period,
                       "minutes"={sec_onds <- 60; to.period(x=mid_prices, period=period)},
                       "3 min"={sec_onds <- 3*60; to.minutes3(x=mid_prices)},
                       "5 min"={sec_onds <- 5*60; to.minutes5(x=mid_prices)},
                       "10 min"={sec_onds <- 10*60; to.minutes10(x=mid_prices)},
                       "15 min"={sec_onds <- 15*60; to.minutes15(x=mid_prices)},
                       "30 min"={sec_onds <- 30*60; to.minutes30(x=mid_prices)},
                       "hours"={sec_onds <- 60*60; to.period(x=mid_prices, period=period)}
                       )  # end switch
# round up times to next period
  index(mid_prices) <- align.time(x=index(mid_prices), n=sec_onds)
  mid_prices
}  # end scrub_agg




#' Calculate returns of either \code{TAQ} or \code{OHLC} data in \code{xts} format.
#' 
#' @param xts_data \code{xts} time series of \code{TAQ} or \code{OHLC} data.
#' @return  \code{xts} time series of returns and volumes.
#' @details The function \code{calc_rets} calculates returns and binds them with
#'   volume data.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"), to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random TAQ prices
#' x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
#' # create vector of random bid-offer prices
#' bid_offer <- abs(rnorm(length(in_dex)))/10
#' # create TAQ data using cbind
#' xts_data <- cbind(x_ts-bid_offer, x_ts+bid_offer)
#' # add Trade.Price
#' xts_data <- cbind(xts_data, x_ts+rnorm(length(in_dex))/10)
#' # add Volume
#' xts_data <- cbind(xts_data, sample(x=10*(2:18), size=length(in_dex), replace=TRUE))
#' colnames(xts_data) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
#' xts_data <- calc_rets(xts_data)
calc_rets <- function(xts_data) {

# return NULL if no data
  if (is.null(xts_data))  return(NULL)
# calculate mid prices
  if(ncol(xts_data)==6)  # TAQ data has 6 columns
    diff_series <- 0.5 * (xts_data[, "Bid.Price"] + xts_data[, "Ask.Price"])
  else
    diff_series <- xts_data[, 4]  # OHLC data
  colnames(diff_series) <- "Mid.Rets"
# calculate returns
  diff_series <- diff(diff_series)/c(1, diff(.index(diff_series)))
  diff_series[1, ] <- 0
# cbind diff_series with volume data
  if(ncol(xts_data)==6)  # TAQ data has 6 columns
    diff_series <- cbind(diff_series, xts_data[, "Volume"])
  else
    diff_series <- cbind(diff_series, xts_data[, 5])
  diff_series
}  # end calc_rets




#' Load, scrub, aggregate, and rbind multiple days of \code{TAQ} data for a
#' single symbol, and save the \code{OHLC} time series to a single \sQuote{\code{*.RData}} file.
#' 
#' @param sym_bol \code{character} string representing symbol or ticker.
#' @param data_dir \code{character} string representing directory containing input \sQuote{\code{*.RData}} files.
#' @param output_dir \code{character} string representing directory containing output \sQuote{\code{*.RData}} files.
#' @inheritParams scrub_agg
#' @return  \code{OHLC} time series in \code{xts} format.
#' @details The function \code{save_scrub_agg} loads multiple days of \code{TAQ} 
#'   data, then scrubs, aggregates, and rbinds them into a \code{OHLC} time 
#'   series, and finally saves it to a single \sQuote{\code{*.RData}} file. The 
#'   \code{OHLC} time series is stored in a variable named 
#'   \sQuote{\code{symbol}}, and then it's saved to a file named 
#'   \sQuote{\code{symbol.RData}} in the \sQuote{\code{output_dir}} directory. 
#'   The \code{TAQ} data files are assumed to be stored in separate directories 
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own 
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}} 
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \code{TAQ} data.
#' @examples
#' \dontrun{
#' save_scrub_agg("SPY")
#' }
save_scrub_agg <- function(sym_bol, 
                      data_dir="E:/mktdata/sec/", 
                      output_dir="E:/output/data/",
                      vol_window=51, 
                      vol_mult=2, 
                      period="minutes", 
                      tzone="America/New_York") {

# create path to directory containing *.RData files
  file_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
  file_list <- list.files(file_dir)
# create paths to *.RData files
  file_names <- file.path(file_dir, file_list)

# load TAQ data one by one, scrub and aggregate it, return list of xts
da_ta <- sapply(file_names, function(file_name) {
  cat("loading", sym_bol, "from file: ", file_name, "\n")
  data_name <- load(file_name)
  scrub_agg(get(data_name),
            vol_window=vol_window, 
            vol_mult=vol_mult, 
            period=period, tzone=tzone)
})  # end sapply

# recursively "rbind" the list into a single xts
  da_ta <- do_call_rbind(da_ta)
# assign column names, i.e. "symbol.High"
  colnames(da_ta) <- lapply(strsplit(colnames(da_ta), split="[.]"), 
                           function(strng) paste(sym_bol, strng[-1], sep="."))

# copy the xts data to a variable with the name 'sym_bol'
  assign(sym_bol, da_ta)

# save the xts data to a file in the output_dir
  save(list=eval(sym_bol), file=file.path(output_dir, paste0(sym_bol, ".RData")))
  invisible(sym_bol)

}  # end save_scrub_agg




#' Load and scrub multiple days of \code{TAQ} data for a single symbol, and save
#' it to multiple \sQuote{\code{*.RData}} files.
#' 
#' @inheritParams save_scrub_agg
#' @return  \code{TAQ} time series in \code{xts} format.
#' @details The function \code{save_TAQ} loads multiple days of \code{TAQ} data,
#'   scrubs it, and saves it to \sQuote{\code{*.RData}} files. It uses the same
#'   file names for output as the input file names.
#'   The \code{TAQ} data files are assumed to be stored in separate directories 
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own 
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}} 
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \code{TAQ} data.
#' @examples
#' \dontrun{
#' save_TAQ("SPY")
#' }
save_TAQ <- function(sym_bol, 
                      data_dir="E:/mktdata/sec/", 
                      output_dir="E:/output/data/",
                      vol_window=51, 
                      vol_mult=2, 
                      tzone="America/New_York") {

# create path to directory containing *.RData files
  data_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
  file_names <- list.files(data_dir)
# create path to directory for writing *.RData files
  output_dir <- file.path(output_dir, sym_bol)

# load TAQ data one-by-one, scrub, and save
  dummy_data <- sapply(file_names, function(file_name) {
    cat("loading", sym_bol, "from file: ", file_name, "\n")
    file_name_in <- file.path(data_dir, file_name)
    data_name <- load(file_name_in)
    file_name_out <- file.path(output_dir, file_name)
# save the xts data to a file in the output_dir
    taq_data <- scrub_TAQ(get(data_name), vol_window=vol_window, vol_mult=vol_mult, tzone=tzone)
    if (!is.null(taq_data)) {
      assign(data_name, taq_data)
      save(list=data_name, file=file_name_out)
      cat("finished saving", sym_bol, "to file: ", file_name, "\n")
    }
    file_name
  })  # end sapply

  invisible(sym_bol)

}  # end save_TAQ




#' Load, scrub, aggregate, and rbind multiple days of \code{TAQ} data for a 
#' single symbol. Calculate returns and save them to a single \sQuote{\code{*.RData}}
#' file.
#' 
#' @inheritParams save_scrub_agg
#' @return  time series of returns and volume in \code{xts} format.
#' @details The function \code{save_rets} loads multiple days of \code{TAQ} 
#'   data, then scrubs, aggregates, and rbinds them into a \code{OHLC} time 
#'   series.  It then calculates returns using function \code{calc_rets}, and
#'   stores them in a variable named \sQuote{\code{symbol.rets}}, and saves them
#'   to a file called \sQuote{\code{symbol.rets.RData}}.
#'   The \code{TAQ} data files are assumed to be stored in separate directories 
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own 
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}} 
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \code{TAQ} data.
#' @examples
#' \dontrun{
#' save_rets("SPY")
#' }
save_rets <- function(sym_bol, 
                      data_dir="E:/mktdata/sec/", 
                      output_dir="E:/output/data/",
                      vol_window=51, 
                      vol_mult=2, 
                      period="minutes", 
                      tzone="America/New_York") {

# create path to directory containing *.RData files
  file_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
  file_list <- list.files(file_dir)
# create paths to *.RData files
  file_names <- file.path(file_dir, file_list)

# load TAQ data into list
  data <- sapply(file_names, function(file_name) {
    cat("loading", sym_bol, "from file: ", file_name, "\n")
    data_name <- load(file_name)
    get(data_name)
  })

# scrub and aggregate the TAQ data
  data <- lapply(data, scrub_agg, vol_window=vol_window, vol_mult=vol_mult, period=period, tzone=tzone)

# calculate returns
  data <- lapply(data, calc_rets)

# recursively "rbind" the list into a single xts
  data <- do_call_rbind(data)
# assign column names, i.e. "symbol.rets"
  colnames(data) <- c(paste(sym_bol, "rets", sep="."), paste(sym_bol, "vol", sep="."))

# copy the xts data to a variable with the name 'sym_bol'
  sym_bol_rets <- paste(sym_bol, "rets", sep=".")
  assign(sym_bol_rets, data)

# save the xts data to a file in the output_dir
  save(list=eval(sym_bol_rets), file=file.path(output_dir, paste0(sym_bol_rets, ".RData")))
  invisible(sym_bol_rets)

}  # end save_rets




#' Load \code{OHLC} time series data for a single symbol, calculate its returns,
#' and save them to a single \sQuote{\code{*.RData}} file, without aggregation.
#' 
#' @inheritParams save_scrub_agg
#' @return  time series of returns and volume in \code{xts} format.
#' @details The function \code{save_rets_ohlc} loads \code{OHLC} time series 
#'   data from a single file.  It then calculates returns using function
#'   \code{calc_rets}, and stores them in a variable named
#'   \sQuote{\code{symbol.rets}}, and saves them to a file called
#'   \sQuote{\code{symbol.rets.RData}}.
#' @examples
#' \dontrun{
#' save_rets_ohlc("SPY")
#' }
save_rets_ohlc <- function(sym_bol, 
                      data_dir="E:/output/data/", 
                      output_dir="E:/output/data/") {

# create path to directory containing sym_bol.RData file
  file_name <- file.path(data_dir, paste0(sym_bol, ".RData"))
# load OHLC data
  cat("loading", sym_bol, "from file: ", file_name, "\n")
  data_name <- load(file_name)

# calculate returns
  da_ta <- calc_rets(get(data_name))

# copy the xts data to a variable with the name 'sym_bol'
  sym_bol_rets <- paste(sym_bol, "rets", sep=".")
  assign(sym_bol_rets, da_ta)

# save the xts data to a file in the output_dir
  cat("saving", sym_bol, "to file: ", paste0(sym_bol_rets, ".RData"), "\n")
  save(list=eval(sym_bol_rets), file=file.path(output_dir, paste0(sym_bol_rets, ".RData")))
  invisible(sym_bol_rets)

}  # end save_rets_ohlc




#' Calculate time series of volatility estimates from a \code{OHLC} time series.
#' 
#' @param log_ohlc \code{OHLC} time series of log prices.
#' @param calc_method \code{character} string representing method for estimating volatility.
#' @return  time series of volatility estimates.
#' @details Calculates volatility estimates from \code{OHLC} values at each
#'   point in time (row).  The methods include Garman-Klass and Rogers-Satchell.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"),
#'               to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random prices
#' x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
#' # add Volume data
#' x_ts <- merge(x_ts, 
#'           volume=sample(x=10*(2:18), 
#'             size=length(in_dex), replace=TRUE))
#' # aggregate to minutes OHLC data
#' oh_lc <- to.period(x=x_ts, period="minutes")
#' # calculate volatility estimates
#' vol_at <- vol_ohlc(oh_lc)
vol_ohlc <- function(log_ohlc, calc_method="rogers.satchell") {
  vol_at <- switch(calc_method,
         "close"={(log_ohlc[, 4]-log_ohlc[, 1])^2},
         "garman.klass"={0.5*(log_ohlc[, 2]-log_ohlc[, 3])^2 - 
                           (2*log(2)-1)*(log_ohlc[, 4]-log_ohlc[, 1])^2},
         "rogers.satchell"={(log_ohlc[, 2]-log_ohlc[, 4])*(log_ohlc[, 2]-log_ohlc[, 1]) + 
                              (log_ohlc[, 3]-log_ohlc[, 4])*(log_ohlc[, 3]-log_ohlc[, 1])},
  )  # end switch
  colnames(vol_at) <- paste0(deparse(substitute(log_ohlc)), ".Volat")
  vol_at
}  # end vol_ohlc




#' Calculate time series of skew estimates from a \code{OHLC} time series.
#' 
#' @param log_ohlc \code{OHLC} time series of log prices.
#' @param calc_method \code{character} string representing method for estimating skew.
#' @return  time series of skew estimates.
#' @details Calculates skew estimates from \code{OHLC} values at each
#'   point in time (row).  The methods include Garman-Klass and Rogers-Satchell.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"),
#'               to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random prices
#' x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
#' # add Volume data
#' x_ts <- merge(x_ts, 
#'           volume=sample(x=10*(2:18), 
#'             size=length(in_dex), replace=TRUE))
#' # aggregate to minutes OHLC data
#' oh_lc <- to.period(x=x_ts, period="minutes")
#' # calculate skew estimates
#' sk_ew <- skew_ohlc(oh_lc)
skew_ohlc <- function(log_ohlc, calc_method="rogers.satchell") {
  sk_ew <- 
    (log_ohlc[, 2]-log_ohlc[, 4])*(log_ohlc[, 2]-log_ohlc[, 1])*(log_ohlc[, 2]-0.5*(log_ohlc[, 4] + log_ohlc[, 1])) + 
    (log_ohlc[, 3]-log_ohlc[, 4])*(log_ohlc[, 3]-log_ohlc[, 1])*(log_ohlc[, 3]-0.5*(log_ohlc[, 4] + log_ohlc[, 1]))
  colnames(sk_ew) <- paste0(deparse(substitute(log_ohlc)), ".Skew")
  sk_ew
}  # end skew_ohlc




#' Calculate the moment of a \code{OHLC} time series.
#' 
#' @param ohlc \code{OHLC} time series of prices.
#' @param mom_fun \code{character} string representing function for estimating
#'   the moment.
#' @param calc_method \code{character} string representing method for estimating
#'   the moment.
#' @param vo_lu \code{logical} should estimate be weighted by trade volume.
#' @return  \code{numeric} value with estimate of the moment.
#' @details Calculates a single number representing the estimate of the moment 
#'   of a \code{OHLC} time series of prices.  By default the estimate is trade
#'   volume weighted.
#' @examples
#' # create time index of one minute intervals over several days
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"),
#'               to=as.POSIXct("2015-02-28 16:00:00"), by="1 min")
#' # create xts of random prices
#' x_ts <- xts(exp(cumsum(0.001*rnorm(length(in_dex)))), order.by=in_dex)
#' # add trade Volume data
#' x_ts <- merge(x_ts, 
#'           volume=sample(x=10*(2:18), 
#'             size=length(in_dex), replace=TRUE))
#' # aggregate to hours OHLC data
#' oh_lc <- to.period(x=x_ts, period="hours")
#' # calculate time series of daily skew estimates
#' sk_ew <- apply.daily(x=oh_lc, FUN=moment_ohlc, mom_fun="skew_ohlc")
moment_ohlc <- function(ohlc, mom_fun="vol_ohlc", calc_method="rogers.satchell", vo_lu=TRUE, ...) {
  log_ohlc <- log(ohlc[, 1:4])
# match "mom_fun" with moment function
  mom_fun <- match.fun(mom_fun)
  mo_ment <- mom_fun(log_ohlc=log_ohlc, calc_method=calc_method)
# weight the estimates by volume
  if (vo_lu) {
    mo_ment <- ohlc[, 5]*mo_ment
    mo_ment <- sum(mo_ment)/sum(ohlc[, 5])
  } else
    mo_ment <- sum(mo_ment)
  mo_ment
}  # end moment_ohlc




#' Calculate the running estimates over time of the moment of a \code{OHLC} time
#' series.
#' 
#' @param ohlc \code{OHLC} time series of prices.
#' @param mom_fun \code{character} string representing function for estimating
#'   the moment.
#' @param calc_method \code{character} string representing method for estimating
#'   the moment.
#' @param n \code{numeric} number of periods for averaging of estimates.
#' @param N \code{numeric} number of periods in a year (to annualize the estimates).
#' @param vo_lu \code{logical} should estimate be weighted by trade volume.
#' @return  \code{numeric} time series of estimates of the moment.
#' @details Calculates a time series of running estimates of the moment of a
#'   \code{OHLC} time series of prices.  By default the estimates are trade 
#'   volume weighted.
#' @examples
#' # create time index of one minute intervals over several days
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"),
#'               to=as.POSIXct("2015-02-28 16:00:00"), by="1 min")
#' # create xts of random prices
#' x_ts <- xts(exp(cumsum(0.001*rnorm(length(in_dex)))), order.by=in_dex)
#' # add trade Volume data
#' x_ts <- merge(x_ts, 
#'           volume=sample(x=10*(2:18), 
#'             size=length(in_dex), replace=TRUE))
#' # aggregate to hours OHLC data
#' oh_lc <- to.period(x=x_ts, period="hours")
#' # calculate time series of running volatility and skew estimates
#' vol_at <- run_moment_ohlc(ohlc=oh_lc)
#' sk_ew <- run_moment_ohlc(ohlc=oh_lc, mom_fun="skew_ohlc")
#' sk_ew <- sk_ew/(vol_at)^(1.5)
#' sk_ew[1, ] <- 0
#' sk_ew <- na.locf(sk_ew)
run_moment_ohlc <- function(ohlc, mom_fun="vol_ohlc", calc="rogers.satchell", n=20, N=260, vo_lu=TRUE, ...) {
  log_ohlc <- log(ohlc[, 1:4])
  # match "mom_fun" with moment function
  mom_fun <- match.fun(mom_fun)
  mo_ment <- mom_fun(log_ohlc=log_ohlc, calc=calc)
  # weight by volume
  if (vo_lu) {
    mo_ment <- ohlc[, 5]*mo_ment
    run_volume <- runSum(ohlc[, 5], n=n)
    mo_ment <- N*runSum(mo_ment, n=n)/(n*run_volume)
  } else
    mo_ment <- N*runSum(mo_ment, n=n)/n
  mo_ment[1:(n-1)] <- 0
  colnames(mo_ment) <- paste(
    strsplit(colnames(ohlc)[1], split="[.]")[[1]][1], 
    "Vol", sep=".")
  mo_ment
}  # end run_moment_ohlc



### Utility Functions

#' Recursively \sQuote{\code{rbind}} a list of objects, such as \code{xts} time series.
#' 
#' Performs the same operation as \code{do.call(rbind, list_var)}, but using 
#' recursion, which is much faster and uses less memory. This is the same 
#' function as \sQuote{\code{\link[qmao]{do.call.rbind}}} from package
#' \sQuote{\href{https://r-forge.r-project.org/R/?group_id=1113}{qmao}}.
#' 
#' @param list_var list of \code{vectors}, \code{matrices}, \code{data
#'   frames}, or \code{time series}.
#' @return  single \code{vector}, \code{matrix}, \code{data frame}, or
#'   \code{time series}.
#' @details Calls lapply in a loop, each time binding neighboring elements and
#'   dividing the length of \code{list_var} by half. The result of performing
#'   \code{do_call_rbind(list_xts)} on a list of \code{xts} time series is
#'   identical to performing \code{do.call(rbind, list_xts)}. But
#'   \code{do.call(rbind, list_xts)} is very slow, and often causes an \sQuote{out of
#'   memory} error.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' # split time series into daily list
#' list_xts <- split(x_ts, "days")
#' # rbind the list back into a time series and compare with the original
#' identical(x_ts, do_call_rbind(list_xts))
do_call_rbind <- function(list_var) {
  while (length(list_var) > 1) {
    # index of odd list elements
    odd_index <- seq(from=1, to=length(list_var), by=2)
    # bind neighboring elements and divide list_var by half
    list_var <- lapply(odd_index, function(in_dex) {
      if (in_dex==length(list_var)) {
        return(list_var[[in_dex]])
      }
      return(rbind(list_var[[in_dex]], list_var[[in_dex+1]]))
    })  # end lapply
  }  # end while
  # list_var has only one element - return it
  list_var[[1]]
}  # end do_call_rbind




#' Calculate the running sum of an \code{xts} time series over a sliding window 
#' (lookback period). 
#' 
#' Performs the same operation as \code{runSum()} from package, 
#' \sQuote{\href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR}}.
#' but using vectorized functions, so it's a little faster.  
#' 
#' @param x_ts an \code{xts} time series containing one or more columns of data.
#' @param win_dow an integer specifying the number of lookback periods. 
#' @return \code{xts} time series with the same dimensions as the input series. 
#' @details For example, if win_dow=3, then the running sum at any point is 
#'   equal to the sum of \code{x_ts} values for that point plus two preceding 
#'   points.
#'   The initial values of run_sum() are equal to cumsum() values, so that
#'   run_sum() doesn't return any NA values.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' foo <- run_sum(x_ts=get("SPY"), win_dow=3)
run_sum <- function(x_ts, win_dow) {
  cum_sum <- cumsum(x_ts)
  out_put <- cum_sum - lag(x=cum_sum, k=win_dow)
  out_put[1:win_dow, ] <- cum_sum[1:win_dow, ]
  out_put
}  # end run_sum




#' Calculate the volume-weighted average price of an \code{OHLC} time series
#' over a sliding window (lookback period).
#' 
#' Performs the same operation as \code{VWAP()} from package, 
#' \sQuote{\href{https://cran.r-project.org/web/packages/TTR/index.html}{VWAP}}.
#' but using vectorized functions, so it's a little faster.  
#' 
#' @param x_ts an \code{xts} time series containing one or more columns of data.
#' @param win_dow an integer specifying the number of lookback periods. 
#' @return \code{xts} time series with a single column containing the VWAP of 
#'   the close prices, with the same number of rows as the input \code{xts} series.
#' @details The volume-weighted average price (VWAP) over a period is defined as
#'   the sum of the prices multiplied by trading volumes, divided by the total
#'   trading volume in that period.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' foo <- v_wap(x_ts=get("SPY"), win_dow=11)
v_wap <- function(x_ts, win_dow) {
  v_wap <- run_sum(x_ts=Cl(x_ts)*Vo(x_ts), win_dow=win_dow)
  vol_ume <- run_sum(x_ts=Vo(x_ts), win_dow=win_dow)
  v_wap <- v_wap/vol_ume
  v_wap[is.na(v_wap)] <- 0
  v_wap
}  # end v_wap




