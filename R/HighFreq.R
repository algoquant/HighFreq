#' Identify extreme values in a univariate \code{xts} time series.
#'
#' Identifies extreme values as those that exceed a multiple of the rolling
#' volatility.
#'
#' @export
#' @param x_ts univariate \code{xts} time series.
#' @param win_dow number of data points for estimating rolling volatility.
#' @param vol_mult volatility multiplier.
#' @return  \code{Boolean} vector with the same number of rows as input time
#'   series.
#' @details Calculates rolling volatility as a quantile of values over a sliding
#'   window. Extreme values are those that exceed the product of the volatility
#'   multiplier times the rolling volatility. Extreme values are the very tips
#'   of the tails when the distribution of values becomes very fat-tailed. The
#'   volatility multiplier \code{vol_mult} controls the threshold at which
#'   values are identified as extreme. Smaller volatility multiplier values will
#'   cause more values to be identified as extreme.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' # scrub extreme values
#' x_ts <- x_ts[!extreme_values(x_ts, vol_mult=1)]

extreme_values <- function(x_ts, win_dow=51, vol_mult=2) {
# calculate volatility as rolling quantile
  vo_lat <- caTools::runquantile(x=abs(as.vector(x_ts)), k=win_dow,
                        probs=0.9, endrule="constant", align="center")
  vo_lat <- xts(vo_lat, order.by=index(x_ts))
  colnames(vo_lat) <- "volat"
# carry forward non-zero volatility values
  vo_lat[vo_lat==0] <- NA
  vo_lat[1] <- 1
  vo_lat <- na.locf(vo_lat)
#  vo_lat <- na.locf(vo_lat, fromLast=TRUE)

# extreme value if x_ts greater than scaled volatility
  ex_treme <- (abs(x_ts) > 2*vol_mult*vo_lat)
  ex_treme[1] <- FALSE
  colnames(ex_treme) <- "suspect"

  cat("date:", format(as.Date(index(first(x_ts)))), "\tfound", sum(ex_treme), "extreme values\n")
  ex_treme
}  # end extreme_values




#' Identify isolated price jumps in a univariate \code{xts} time series of
#' prices, based on pairs of large neighboring returns of opposite sign.
#'
#' @export
#' @param x_ts univariate \code{xts} time series of prices.
#' @inheritParams extreme_values
#' @return  \code{Boolean} vector with the same number of rows as input time
#'   series.
#' @details Isolated price jumps are single prices that are very different from
#'   neighboring values.  Price jumps create pairs of large neighboring returns
#'   of opposite sign. The function \code{price_jumps} first calculates simple
#'   returns from prices. Then it calculates the rolling volatility of returns
#'   as a quantile of returns over a sliding window. Jump prices are identified
#'   as those where neighboring returns both exceed a multiple of the rolling
#'   volatility, but the sum of those returns doesn't exceed it.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' # scrub jump prices
#' x_ts <- x_ts[!price_jumps(x_ts, vol_mult=1.0)]

price_jumps <- function(x_ts, win_dow=51, vol_mult=2) {
# calculate simple returns
  diff_series <- rutils:::diff.xts(x_ts)
  diff_series[1, ] <- 0
  colnames(diff_series) <- "diffs"
  diff_series_fut <- lag(diff_series, -1)
  diff_series_fut[NROW(diff_series_fut)] <- 0
  colnames(diff_series_fut) <- "diff_series_fut"

# calculate vo_lat as rolling quantile
  vo_lat <- caTools::runquantile(x=abs(as.vector(diff_series)), k=win_dow,
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
  cat("date:", format(as.Date(index(first(x_ts)))), "\tfound", sum(sus_pect), "jump prices\n")
  sus_pect
}  # end price_jumps




#' Scrub a single day of \code{TAQ} data in \code{xts} format, without
#' aggregation.
#'
#' @export
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
#' taq_data <- scrub_TAQ(taq_data, win_dow=11, vol_mult=1)

scrub_TAQ <- function(taq_data, win_dow=51, vol_mult=2, tzone="America/New_York") {
# convert timezone of index to New_York
  index(taq_data) <- with_tz(time=index(taq_data), tzone=tzone)
# subset data to NYSE trading hours
  taq_data <- taq_data['T09:30:00/T16:00:00', ]
# return NULL if no data
  if (NROW(taq_data)==0)  return(NULL)
#  to_day <- as.Date(index(first(taq_data)))

# remove duplicate time stamps using 'duplicated'
  taq_data <- taq_data[!duplicated(index(taq_data)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- taq_data[, "Ask.Price"] - taq_data[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- extreme_values(bid_offer, win_dow=win_dow, vol_mult=vol_mult)
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
  taq_data[price_jumps(mid_prices, win_dow=win_dow, vol_mult=vol_mult), ] <- NA
  na.locf(taq_data)
}  # end scrub_TAQ




#' Scrub a single day of \code{TAQ} data, aggregate it, and convert to
#' \code{OHLC} format.
#'
#' @export
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
                      win_dow=51,
                      vol_mult=2,
                      period="minutes",
                      tzone="America/New_York") {
# convert timezone of index to New_York
  index(taq_data) <- with_tz(time=index(taq_data), tzone=tzone)
# subset data to NYSE trading hours
  taq_data <- taq_data['T09:30:00/T16:00:00', ]
# return NULL if no data
  if (NROW(taq_data)==0)  return(NULL)
#  to_day <- as.Date(index(first(taq_data)))

# remove duplicate time stamps using 'duplicated'
  taq_data <- taq_data[!duplicated(index(taq_data)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- taq_data[, "Ask.Price"] - taq_data[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- extreme_values(bid_offer, win_dow=win_dow, vol_mult=vol_mult)
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
  mid_prices[price_jumps(mid_prices, win_dow=win_dow, vol_mult=vol_mult)] <- NA
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
#' @export
#' @param x_ts \code{xts} time series of \code{TAQ} or \code{OHLC} data.
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
#' x_ts <- cbind(x_ts-bid_offer, x_ts+bid_offer)
#' # add Trade.Price
#' x_ts <- cbind(x_ts, x_ts+rnorm(length(in_dex))/10)
#' # add Volume
#' x_ts <- cbind(x_ts, sample(x=10*(2:18), size=length(in_dex), replace=TRUE))
#' colnames(x_ts) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
#' x_ts <- calc_rets(x_ts)

calc_rets <- function(x_ts) {
# return NULL if no data
  if (is.null(x_ts))  return(NULL)
# calculate mid prices
  if(NCOL(x_ts)==6)  # TAQ data has 6 columns
    diff_series <- 0.5 * (x_ts[, "Bid.Price"] + x_ts[, "Ask.Price"])
  else
    diff_series <- x_ts[, 4]  # OHLC data
  colnames(diff_series) <- "Mid.Rets"
# calculate returns
  diff_series <- rutils:::diff.xts(diff_series)/c(1, diff(.index(diff_series)))
  diff_series[1, ] <- 0
# cbind diff_series with volume data
  if(NCOL(x_ts)==6)  # TAQ data has 6 columns
    diff_series <- cbind(diff_series, x_ts[, "Volume"])
  else
    diff_series <- cbind(diff_series, x_ts[, 5])
  diff_series
}  # end calc_rets




#' Load, scrub, aggregate, and rbind multiple days of \code{TAQ} data for a
#' single symbol, and save the \code{OHLC} time series to a single
#' \sQuote{\code{*.RData}} file.
#'
#' @export
#' @param sym_bol \code{character} string representing symbol or ticker.
#' @param data_dir \code{character} string representing directory containing
#'   input \sQuote{\code{*.RData}} files.
#' @param output_dir \code{character} string representing directory containing
#'   output \sQuote{\code{*.RData}} files.
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
#'   \sQuote{\code{*.RData}} files, each file containing one day of \code{TAQ}
#'   data.
#' @examples
#' \dontrun{
#' save_scrub_agg("SPY")
#' }

save_scrub_agg <- function(sym_bol,
                      data_dir="E:/mktdata/sec/",
                      output_dir="E:/output/data/",
                      win_dow=51,
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
da_ta <- lapply(file_names, function(file_name) {
  cat("loading", sym_bol, "from file: ", file_name, "\n")
  data_name <- load(file_name)
  scrub_agg(get(data_name),
            win_dow=win_dow,
            vol_mult=vol_mult,
            period=period, tzone=tzone)
})  # end sapply

# recursively "rbind" the list into a single xts
  da_ta <- caTools::do_call_rbind(da_ta)
# assign column names, i.e. "symbol.High"
  colnames(da_ta) <- sapply(strsplit(colnames(da_ta), split="[.]"),
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
#' @export
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
                      win_dow=51,
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
    taq_data <- scrub_TAQ(get(data_name), win_dow=win_dow, vol_mult=vol_mult, tzone=tzone)
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
#' @export
#' @inheritParams save_scrub_agg
#' @return  Time series of returns and volume in \code{xts} format.
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
                      win_dow=51,
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
  taq_data <- sapply(file_names, function(file_name) {
    cat("loading", sym_bol, "from file: ", file_name, "\n")
    data_name <- load(file_name)
    get(data_name)
  })

# scrub and aggregate the TAQ data
  ohlc_data <- lapply(taq_data, scrub_agg, win_dow=win_dow, vol_mult=vol_mult, period=period, tzone=tzone)

# calculate returns
  ohlc_data <- lapply(ohlc_data, calc_rets)

# recursively "rbind" the list into a single xts
  ohlc_data <- caTools::do_call_rbind(ohlc_data)
# assign column names, i.e. "symbol.rets"
  colnames(ohlc_data) <- c(paste(sym_bol, "rets", sep="."), paste(sym_bol, "vol", sep="."))

# copy the xts data to a variable with the name 'sym_bol'
  sym_bol_rets <- paste(sym_bol, "rets", sep=".")
  assign(sym_bol_rets, ohlc_data)

# save the xts data to a file in the output_dir
  save(list=eval(sym_bol_rets), file=file.path(output_dir, paste0(sym_bol_rets, ".RData")))
  invisible(sym_bol_rets)

}  # end save_rets




#' Load \code{OHLC} time series data for a single symbol, calculate its returns,
#' and save them to a single \sQuote{\code{*.RData}} file, without aggregation.
#'
#' @export
#' @inheritParams save_scrub_agg
#' @return  Time series of returns and volume in \code{xts} format.
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




#' Adjusts an \code{OHLC} time series to make open prices equal to the close
#' prices from the previous period.
#'
#' @export
#' @param x_ts an \code{xts} time series containing one or more columns of data.
#' @return \code{xts} \code{OHLC} time series with open prices equal to the
#'   close prices from the previous period.
#' @details Adds or subtracts a price adjustment to make all open prices equal
#'   to the close prices from the previous period.  The adjustment preserves the
#'   price differences within each bar of \code{OHLC} prices, and so preserves
#'   open to close returns, variance estimates, etc.
#' @examples
#' # define end points at 10-minute intervals (SPY is minutely bars)
#' end_points <- rutils::end_points(SPY["2009"], inter_val=10)
#' # aggregate over 10-minute end_points:
#' open_close(x_ts=SPY["2009"], end_points=end_points)
#' # aggregate over days:
#' open_close(x_ts=SPY["2009"], period="days")
#' # equivalent to:
#' to.period(x=SPY["2009"], period="days", name=rutils::na_me(SPY))

open_close <- function(x_ts) {
  op_en <- Op(x_ts)
  clo_se <- lag.xts(Cl(x_ts), k=-1)
  which(!(op_en==clo_se))
}  # end open_close




#' Calculate a time series of variance estimates for an \code{OHLC} time series.
#'
#' Calculates variance estimates for each bar of \code{OHLC} prices at each
#' point in time (row), using the squared differences of \code{OHLC} prices at
#' each point in time.
#' 
#' @export
#' @param ohlc \code{OHLC} time series of prices.
#' @param calc_method \code{character} string representing method for estimating
#'   variance.  The methods include: 
#' \itemize{
#'   \item "close" close to close, 
#'   \item "garman.klass" Garman-Klass, 
#'   \item "garman.klass_yz" Garman-Klass with account for close-to-open jumps, 
#'   \item "rogers.satchell" Rogers-Satchell, 
#'   \item "yang.zhang" Yang-Zhang, 
#' }
#' @return  A single-column \code{xts} time series of variance estimates, with
#'   the same number of rows as the input time series.
#' @details Performs a similar operation as function \code{volatility()} from
#'   package \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR},
#'   but without calculating a running sum using \code{runSum()}.  It's also a 
#'   little faster because it performs less data validation.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"),
#'               to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random prices
#' x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
#' # aggregate to minutes OHLC data
#' oh_lc <- to.period(x=x_ts, period="minutes")
#' # calculate variance estimates
#' vari_ance <- vari_ance(oh_lc)
#' # calculate variance estimates for SPY
#' vari_ance <- vari_ance(SPY, calc_method="yang.zhang")

vari_ance <- function(ohlc, calc_method="garman.klass_yz") {
  sym_bol <- deparse(substitute(ohlc))
  ohlc <- log(ohlc[, 1:4])
  vari_ance <- switch(calc_method,
         "close"={(ohlc[, 4]-ohlc[, 1])^2},
         "garman.klass"={0.5*(ohlc[, 2]-ohlc[, 3])^2 -
                           (2*log(2)-1)*(ohlc[, 4]-ohlc[, 1])^2},
         "rogers.satchell"={(ohlc[, 2]-ohlc[, 4])*(ohlc[, 2]-ohlc[, 1]) +
                              (ohlc[, 3]-ohlc[, 4])*(ohlc[, 3]-ohlc[, 1])},
         "garman.klass_yz"={(ohlc[, 1]-rutils:::lag.xts(ohlc[, 4]))^2 +
                            0.5*(ohlc[, 2]-ohlc[, 3])^2 -
                            (2*log(2)-1)*(ohlc[, 4]-ohlc[, 1])^2},
         "yang.zhang"={c_o <- ohlc[, 1]-rutils:::lag.xts(ohlc[, 4]); 
                       o_c <- ohlc[, 1]-ohlc[, 4]; 
                       (c_o-sum(c_o)/NROW(c_o))^2 + 
                       0.67*(o_c-sum(o_c)/NROW(o_c))^2 + 
                       0.33*((ohlc[, 2]-ohlc[, 4])*(ohlc[, 2]-ohlc[, 1]) + 
                               (ohlc[, 3]-ohlc[, 4])*(ohlc[, 3]-ohlc[, 1]))}
  )  # end switch
  colnames(vari_ance) <- paste0(sym_bol, ".Variance")
  vari_ance
}  # end vari_ance




#' Calculate time series of skew estimates from a \code{OHLC} time series.
#'
#' @export
#' @param ohlc \code{OHLC} time series of prices.
#' @param calc_method \code{character} string representing method for estimating
#'   skew.
#' @return  Time series of skew estimates.
#' @details Calculates skew estimates from \code{OHLC} prices at each
#'   point in time (row).  The methods include Garman-Klass and Rogers-Satchell.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"),
#'               to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random prices
#' x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
#' # aggregate to minutes OHLC data
#' oh_lc <- to.period(x=x_ts, period="minutes")
#' # calculate skew estimates
#' sk_ew <- skew_ohlc(oh_lc)

skew_ohlc <- function(ohlc, calc_method="rogers.satchell") {
  sym_bol <- deparse(substitute(ohlc))
  ohlc <- log(ohlc[, 1:4])
  sk_ew <-
    (ohlc[, 2]-ohlc[, 4])*(ohlc[, 2]-ohlc[, 1])*(ohlc[, 2]-0.5*(ohlc[, 4] + ohlc[, 1])) +
    (ohlc[, 3]-ohlc[, 4])*(ohlc[, 3]-ohlc[, 1])*(ohlc[, 3]-0.5*(ohlc[, 4] + ohlc[, 1]))
  colnames(sk_ew) <- paste0(sym_bol, ".Skew")
  sk_ew
}  # end skew_ohlc




#' Calculate time series of Hurst exponent estimates for a \code{OHLC} time
#' series.
#'
#' @export
#' @param ohlc \code{OHLC} time series of prices.
#' @param calc_method \code{character} string representing method for estimating
#'   Hurst exponent.
#' @return  Time series of Hurst exponent estimates.
#' @details Calculates Hurst exponent estimates from \code{OHLC} prices at each
#'   point in time (bar).  The methods include Garman-Klass and Rogers-Satchell.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"),
#'               to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random prices
#' x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
#' # aggregate to minutes OHLC data
#' oh_lc <- to.period(x=x_ts, period="minutes")
#' # calculate Hurst exponent
#' hurst_exp <- hurst_ohlc(oh_lc)

hurst_ohlc <- function(ohlc, calc_method="rogers.satchell") {
  hurst_exp <- switch(calc_method,
                   "close"={(ohlc[, 4]-ohlc[, 1])^2},
                   "garman.klass"={0.5*(ohlc[, 2]-ohlc[, 3])^2 -
                       (2*log(2)-1)*(ohlc[, 4]-ohlc[, 1])^2},
                   "rogers.satchell"={(ohlc[, 2]-ohlc[, 3])/abs(ohlc[, 4]-ohlc[, 1])}
  )  # end switch
  hurst_exp <- ifelse(ohlc[, 4]==ohlc[, 1], 0, log(hurst_exp))
  colnames(hurst_exp) <- paste0(deparse(substitute(ohlc)), ".Hurst")
  hurst_exp
}  # end hurst_ohlc




#' Calculate the aggregation (weighted sum) of a statistical estimator over a
#' \code{OHLC} time series.
#'
#' @export
#' @param ohlc \code{OHLC} time series of prices and trading volumes.
#' @param esti_mator \code{character} string representing function for
#'   estimating the moment.
#' @param weight_ed \code{Boolean} should estimate be weighted by trade volume 
#'   (default is \code{TRUE})?
#' @param ... additional parameters to esti_mator function.
#' @return  Single \code{numeric} value equal to the weighted sum of an
#'   estimator over the time series.
#' @details Calculates a single number representing the weighted sum of an
#'   estimator over the \code{OHLC} time series of prices.  By default the sum
#'   is trade volume weighted.
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
#' sk_ew <- apply.daily(x=oh_lc, FUN=agg_regate, esti_mator="skew_ohlc")

agg_regate <- function(ohlc, esti_mator="vari_ance", weight_ed=TRUE, ...) {
# match "esti_mator" with moment function
  esti_mator <- match.fun(esti_mator)
  agg_regations <- esti_mator(ohlc=ohlc, ...)
# weight the estimates by volume
  if (weight_ed) {
    agg_regations <- ohlc[, 5]*agg_regations
    agg_regations <- sum(agg_regations)/sum(ohlc[, 5])
  } else
    agg_regations <- sum(agg_regations)
  agg_regations
}  # end agg_regate




#' Calculate the rolling aggregations of a statistical estimator over a
#' \code{OHLC} time series.
#'
#' @export
#' @param ohlc \code{OHLC} time series of prices and trading volumes.
#' @param esti_mator \code{character} string representing function for
#'   estimating the moment.
#' @param n \code{numeric} number of periods for averaging of estimates.
#' @param N \code{numeric} number of periods in a year (to annualize the
#'   estimates).
#' @param weight_ed \code{Boolean} should estimate be weighted by trade volume
#'   (default \code{TRUE})?
#' @param ... additional parameters to esti_mator function.
#' @return  \code{numeric} time series of rolling aggregations, with the same
#'   number of rows as the input \code{xts} series.
#' @details Calculates a time series of rolling aggregations of an estimator
#'   over a \code{OHLC} time series of prices or returns, etc.  By default the
#'   estimates are trade volume weighted.
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
#' # calculate time series of rolling variance and skew estimates
#' vari_ance <- roll_agg(ohlc=oh_lc)
#' sk_ew <- roll_agg(ohlc=oh_lc, esti_mator="skew_ohlc")
#' sk_ew <- sk_ew/(vari_ance)^(1.5)
#' sk_ew[1, ] <- 0
#' sk_ew <- na.locf(sk_ew)

roll_agg <- function(ohlc, esti_mator="vari_ance",
                          n=20, N=260, weight_ed=TRUE, ...) {
# match "esti_mator" with moment function
  esti_mator <- match.fun(esti_mator)
  agg_regations <- esti_mator(ohlc=ohlc, ...)
# weight by volume
  if (weight_ed) {
    agg_regations <- ohlc[, 5]*agg_regations
    roll_volume <- roll_sum(ohlc[, 5], win_dow=n)
    agg_regations <- N*roll_sum(agg_regations, win_dow=n)/roll_volume
  } else
    agg_regations <- N*roll_sum(agg_regations, win_dow=n)/n
  agg_regations[1:(n-1)] <- 0
  colnames(agg_regations) <- paste(rutils::na_me(ohlc), "Vol", sep=".")
  agg_regations
}  # end roll_agg



### Utility Functions

#' Calculate the rolling sum of an \code{xts} time series over a sliding window
#' (lookback period).
#'
#' Performs the same operation as function \code{runSum()} from package
#' \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR},
#' but using vectorized functions, so it's a little faster.
#'
#' @export
#' @param x_ts an \code{xts} time series containing one or more columns of data.
#' @param win_dow an integer specifying the number of lookback periods.
#' @return \code{xts} time series with the same dimensions as the input series.
#' @details For example, if win_dow=3, then the rolling sum at any point is
#'   equal to the sum of \code{x_ts} values for that point plus two preceding
#'   points.
#'   The initial values of roll_sum() are equal to cumsum() values, so that
#'   roll_sum() doesn't return any NA values.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' roll_sum(x_ts=get("SPY"), win_dow=3)

roll_sum <- function(x_ts, win_dow) {
  cum_sum <- cumsum(x_ts)
  out_put <- cum_sum - lag(x=cum_sum, k=win_dow)
  out_put[1:win_dow, ] <- cum_sum[1:win_dow, ]
  out_put
}  # end roll_sum




#' Perform daily, weekly, monthly, and yearly seasonality aggregations over a
#' univariate \code{xts} time series.
#'
#' @export
#' @param x_ts univariate \code{xts} time series.
#' @return \code{xts} time series with aggregations over the seasonality
#'   interval.
#' @details An example of a daily seasonality aggregation is the average price
#'   of a stock between 9:30AM and 10:00AM every day, over many days.
#' @examples
#' season_ality(x_ts=get("SPY"))

season_ality <- function(x_ts) {
  in_dex <- format(index(x_ts), "%H:%M")
  agg_regation <- tapply(X=x_ts, INDEX=in_dex, FUN=mean)
  agg_regation <- structure(as.vector(agg_regation), names=names(agg_regation))
  agg_regation <- xts(x=agg_regation,
      order.by=as.POSIXct(paste(Sys.Date(), unique(in_dex))))
  colnames(agg_regation) <- colnames(x_ts)
  agg_regation
}  # end season_ality




#' Calculate the volume-weighted average price of an \code{OHLC} time series
#' over a sliding window (lookback period).
#'
#' Performs the same operation as function \code{VWAP()} from package
#' \href{https://cran.r-project.org/web/packages/TTR/index.html}{VWAP},
#' but using vectorized functions, so it's a little faster.
#'
#' @export
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
#' v_wap(x_ts=get("SPY"), win_dow=11)

v_wap <- function(x_ts, win_dow) {
  v_wap <- roll_sum(x_ts=Cl(x_ts)*Vo(x_ts), win_dow=win_dow)
  vol_ume <- roll_sum(x_ts=Vo(x_ts), win_dow=win_dow)
  v_wap <- v_wap/vol_ume
  v_wap[is.na(v_wap)] <- 0
  v_wap
}  # end v_wap




#' Aggregates an \code{OHLC} time series to lower periodicity.
#'
#' Given an \code{OHLC} time series at high periodicity (say seconds), 
#' calculates the \code{OHLC} prices at lower periodicity (say minutes).
#'
#' @export
#' @param x_ts an \code{xts} time series containing one or more columns of data.
#' @param period aggregation interval ("seconds", "minutes", "hours", "days",
#'   "weeks", "months", "quarters", and "years").
#' @param k number of periods to aggregate over (for example if period="minutes"
#'   and k=2, then aggregate over two minute intervals.)
#' @param end_points an integer vector of end points.
#' @return \code{xts} \code{OHLC} time series of lower periodicity defined by
#'   end_points.
#' @details #' Performs a similar aggregation as function \code{to.period()}
#'   from package 
#'   \href{https://cran.r-project.org/web/packages/xts/index.html}{xts}, but has
#'   the flexibility to aggregate to a user-specified vector of end points. The
#'   function \code{to_period()} simply calls the compiled function 
#'   \code{toPeriod()} (from package 
#'   \href{https://cran.r-project.org/web/packages/xts/index.html}{xts}), to 
#'   perform the actual aggregations.  If \code{end_points} are passed in 
#'   explicitly, then the \code{period} argument is ignored.
#' @examples
#' # define end points at 10-minute intervals (SPY is minutely bars)
#' end_points <- rutils::end_points(SPY["2009"], inter_val=10)
#' # aggregate over 10-minute end_points:
#' to_period(x_ts=SPY["2009"], end_points=end_points)
#' # aggregate over days:
#' to_period(x_ts=SPY["2009"], period="days")
#' # equivalent to:
#' to.period(x=SPY["2009"], period="days", name=rutils::na_me(SPY))

to_period <- function(x_ts,
                      period="minutes", k=1,
                      end_points=xts::endpoints(x_ts, period, k)) {
  .Call("toPeriod", x_ts, as.integer(end_points), TRUE, NCOL(x_ts),
        FALSE, FALSE, colnames(x_ts), PACKAGE="xts")
}  # end to_period



