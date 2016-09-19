#' Calculate a random \code{TAQ} time series of prices and trading volumes, in
#' \code{xts} format.
#'
#' Calculate a \code{TAQ} time series of prices and trading volumes, using
#' random log-normal prices and a time index.
#'
#' @export
#' @param in_dex time index for the \code{TAQ} time series.
#' @param bid_offer the bid-offer spread expressed as a fraction of the prices.
#'   The default value is equal to 0.001 (10bps).
#' @return An \code{xts} time series, with time index equal to the input
#'   \code{in_dex} time index, and with four columns containing the bid, ask,
#'   and trade prices, and the trade volume.
#' @details The function \code{random_taq()} calculates an \code{xts} time
#'   series with four columns containing random log-normal prices: the bid, ask,
#'   and trade prices, and the trade volume.
#'   If \code{in_dex} isn't supplied as an argument, then by default it's
#'   equal to the secondly index over the two previous calendar days.
#' @examples
#' # create secondly TAQ time series of random prices
#' ta_q <- HighFreq::random_taq()
#' # create random TAQ time series from SPY index
#' ta_q <- HighFreq::random_taq(in_dex=SPY["2012-02-13/2012-02-15"])

random_taq <- function(
  in_dex=seq(from=as.POSIXct(paste(Sys.Date()-3, "09:30:00")),
             to=as.POSIXct(paste(Sys.Date()-1, "16:00:00")), by="1 sec"),
  bid_offer=0.001, ...) {
  # create synthetic xts of random log-normal prices
  ta_q <- xts(exp(cumsum(0.001*rnorm(length(in_dex)))), order.by=in_dex)
  # create vector of random bid-offer spreads
  bid_offer <- bid_offer*(1 + runif(length(in_dex)))/2
  # create TAQ data from bid and offer prices
  ta_q <- merge(ta_q*(1-bid_offer), ta_q*(1+bid_offer))
  # add traded price to TAQ data
  r_unif <- runif(length(in_dex))
  ta_q <- merge(ta_q, r_unif*ta_q[, 1] + (1-r_unif)*ta_q[, 2])
  # add trade volume column
  ta_q <- merge(ta_q, sample(x=10*(2:18), size=length(in_dex), replace=TRUE))
  colnames(ta_q) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
  ta_q
}  # end random_taq




#' Calculate a random \code{OHLC} time series of prices and trading volumes, in
#' \code{xts} format.
#'
#' Calculate a random \code{OHLC} time series of prices and trading volumes,
#' either by generating random log-normal prices, or by randomly sampling from
#' an input time series.
#'
#' @export
#' @param oh_lc \code{OHLC} time series of prices and trading volumes, in
#'   \code{xts} format.
#' @param re_duce \code{Boolean} should \code{oh_lc} time series be transformed
#'   to reduced form? (default is \code{TRUE})
#' @return An \code{xts} time series with the same dimensions and the same time
#'   index as the input \code{oh_lc} time series.
#' @details If the input \code{oh_lc} time series is \code{NULL} (the default),
#'   then a synthetic minutely \code{OHLC} time series of random log-normal
#'   prices is calculated, over the two previous calendar days.
#'   If the input \code{oh_lc} time series is not \code{NULL}, then the rows of
#'   \code{oh_lc} are randomly sampled, to produce a random time series.
#'   If \code{re_duce} is \code{TRUE} (the default), then the \code{oh_lc} time
#'   series is first transformed to reduced form, then randomly sampled, and
#'   finally converted to standard form.
#'   Note: randomly sampling from an intraday time series over multiple days
#'   will cause the overnight price jumps to be re-arranged into intraday price
#'   jumps.  This will cause moment estimates to become inflated compared to the
#'   original time series.
#' @examples
#' # create minutely synthetic OHLC time series of random prices
#' oh_lc <- HighFreq::random_ohlc()
#' # create random time series from SPY by randomly sampling it
#' oh_lc <- HighFreq::random_ohlc(oh_lc=SPY["2012-02-13/2012-02-15"])

random_ohlc <- function(oh_lc=NULL, re_duce=TRUE, ...) {
  if (is.null(oh_lc)) {
    # create time index of one second intervals over several days
    in_dex <- seq(from=as.POSIXct(paste(Sys.Date()-3, "09:30:00")),
                  to=as.POSIXct(paste(Sys.Date()-1, "16:00:00")), by="1 sec")
    # create synthetic xts of random log-normal prices
    x_ts <- xts(exp(cumsum(0.01/sqrt(86400)*rnorm(length(in_dex)))), order.by=in_dex)
    # add trade volume column
    x_ts <- merge(x_ts, volume=sample(x=10*(2:18), size=length(in_dex), replace=TRUE))
    # aggregate to minutes OHLC data
    to.period(x=x_ts, period="minutes")
  } else {
    oh_lc <- log(oh_lc)  # transform to normal
    if (re_duce)  # calculate reduced form of oh_lc
      oh_lc <- rutils::diff_ohlc(oh_lc)
    # randomly sample from the rows of oh_lc
    oh_lc <- xts(coredata(oh_lc)[c(1, sample(x=2:NROW(oh_lc), replace=TRUE)), ], order.by=index(oh_lc))
    # return standard form of randomized oh_lc
    exp(rutils::diff_ohlc(oh_lc, re_duce=FALSE))
  }
}  # end random_ohlc




#' Calculate single period returns from either \code{TAQ} or \code{OHLC} prices.
#'
#' @export
#' @param x_ts \code{xts} time series of either \code{TAQ} or \code{OHLC} data.
#' @param col_umn the column number to extract from the \code{OHLC} data.
#'   (default is \code{4}, or the \code{close} prices column)
#' @return A single-column \code{xts} time series of returns.
#' @details Calculates single period returns for either \code{TAQ} or
#'   \code{OHLC} data, as the ratio of the differenced prices divided by the
#'   time index differences.
#'   Identifies the \code{x_ts} time series as \code{TAQ} data when it has six
#'   columns, otherwise assumes it's \code{OHLC} data.
#'   By default, for \code{OHLC} data, it differences the \code{close} prices,
#'   but can also difference other prices depending on the value of
#'   \code{col_umn}.
#' @examples
#' # calculate close to close returns
#' re_turns <- HighFreq::run_returns(x_ts=SPY)
#' # calculate open to open returns
#' re_turns <- HighFreq::run_returns(x_ts=SPY, col_umn=1)

run_returns <- function(x_ts, col_umn=4) {
  # return NULL if no data
  if (is.null(x_ts))  return(NULL)
  # calculate mid prices
  if(NCOL(x_ts)==6)  # TAQ data has 6 columns
    re_turns <- 0.5 * (x_ts[, "Bid.Price"] + x_ts[, "Ask.Price"])
  else
    re_turns <- x_ts[, col_umn]  # OHLC data
  # calculate returns
  re_turns <- 86400*rutils::diff_xts(re_turns)/c(1, diff(.index(re_turns)))
  re_turns[1, ] <- 0
  colnames(re_turns) <- paste0(rutils::na_me(x_ts), ".returns")
  re_turns
}  # end run_returns




#' Identify extreme values in a single-column \code{xts} time series.
#'
#' Identifies extreme values as those that exceed a multiple of the rolling
#' volatility.
#'
#' @export
#' @param x_ts single-column \code{xts} time series.
#' @param win_dow number of data points for estimating rolling volatility.
#' @param vol_mult volatility multiplier.
#' @return A \code{Boolean} vector with the same number of rows as input time
#'   series.
#' @details Calculates the rolling volatility as a quantile of values over a
#'   rolling window. Extreme values are those that exceed the product of the
#'   volatility multiplier times the rolling volatility. Extreme values are the
#'   very tips of the tails when the distribution of values becomes very
#'   fat-tailed. The volatility multiplier \code{vol_mult} controls the
#'   threshold at which values are identified as extreme. Smaller volatility
#'   multiplier values will cause more values to be identified as extreme.
#' @examples
#' # create local copy of SPY TAQ data
#' ta_q <- SPY_TAQ
#' # scrub quotes with suspect bid-offer spreads
#' bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#' sus_pect <- extreme_values(bid_offer, win_dow=win_dow, vol_mult=vol_mult)
#' # remove suspect values
#' ta_q <- ta_q[!sus_pect]

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




#' Identify isolated price jumps in a single-column \code{xts} time series of
#' prices, based on pairs of large neighboring returns of opposite sign.
#'
#' @export
#' @param x_ts single-column \code{xts} time series of prices.
#' @inheritParams extreme_values
#' @return A \code{Boolean} vector with the same number of rows as input time
#'   series.
#' @details Isolated price jumps are single prices that are very different from
#'   neighboring values.  Price jumps create pairs of large neighboring returns
#'   of opposite sign. The function \code{price_jumps()} first calculates simple
#'   returns from prices. Then it calculates the rolling volatility of returns
#'   as a quantile of returns over a rolling window. Jump prices are identified
#'   as those where neighboring returns both exceed a multiple of the rolling
#'   volatility, but the sum of those returns doesn't exceed it.
#' @examples
#' # create local copy of SPY TAQ data
#' ta_q <- SPY_TAQ
#' # calculate mid prices
#' mid_prices <- 0.5 * (ta_q[, "Bid.Price"] + ta_q[, "Ask.Price"])
#' # replace whole rows containing suspect price jumps with NA, and perform locf()
#' ta_q[price_jumps(mid_prices, win_dow=31, vol_mult=1.0), ] <- NA
#' ta_q <- na.locf(ta_q)

price_jumps <- function(x_ts, win_dow=51, vol_mult=2) {
# calculate simple returns
  diff_series <- rutils::diff_xts(x_ts)
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
# cat("Parsing", deparse(substitute(ta_q)), "\n")
# cat("Parsing", strsplit(deparse(substitute(ta_q)), split="[.]")[[1]][4], "on date:", format(to_day), "\tfound", sum(sus_pect), "suspect prices\n")
  cat("date:", format(as.Date(index(first(x_ts)))), "\tfound", sum(sus_pect), "jump prices\n")
  sus_pect
}  # end price_jumps




#' Scrub a single day of \code{TAQ} data in \code{xts} format, without
#' aggregation.
#'
#' @export
#' @inheritParams extreme_values
#' @param ta_q \code{TAQ} time series in \code{xts} format.
#' @param tzone timezone to convert.
#' @return A \code{TAQ} time series in \code{xts} format.
#' @details The function \code{scrub_taq()} performs the same scrubbing
#'   operations as \code{scrub_agg}, except it doesn't aggregate, and returns
#'   the \code{TAQ} data in \code{xts} format.
#' @examples
# scrub a single day of TAQ data without aggregating it
#' ta_q <- HighFreq::scrub_taq(ta_q=SPY_TAQ, win_dow=11, vol_mult=1)
#' # create random TAQ prices and scrub them
#' ta_q <- HighFreq::random_taq()
#' ta_q <- HighFreq::scrub_taq(ta_q=ta_q)
#' ta_q <- HighFreq::scrub_taq(ta_q=ta_q, win_dow=11, vol_mult=1)

scrub_taq <- function(ta_q, win_dow=51, vol_mult=2, tzone="America/New_York") {
# convert timezone of index to New_York
  index(ta_q) <- lubridate::with_tz(time=index(ta_q), tzone=tzone)
# subset data to NYSE trading hours
  ta_q <- ta_q['T09:30:00/T16:00:00', ]
# return NULL if no data
  if (NROW(ta_q)==0)  return(NULL)
#  to_day <- as.Date(index(first(ta_q)))

# remove duplicate time stamps using duplicated()
  ta_q <- ta_q[!duplicated(index(ta_q)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- extreme_values(bid_offer, win_dow=win_dow, vol_mult=vol_mult)
# remove suspect values
  ta_q <- ta_q[!sus_pect]
# replace suspect values
# ta_q[sus_pect, "Bid.Price"] <- ta_q[sus_pect, "Trade.Price"]
# ta_q[sus_pect, "Ask.Price"] <- ta_q[sus_pect, "Trade.Price"]

# scrub quotes with suspect price jumps
# calculate mid prices
  mid_prices <- 0.5 * (ta_q[, "Bid.Price"] + ta_q[, "Ask.Price"])
#  mid_prices <- na.omit(mid_prices)
#  colnames(mid_prices) <- "Mid.Price"
# replace NA volumes with zero
  ta_q[is.na(ta_q[, "Volume"]), "Volume"] <- 0
# replace whole rows containing suspect price jumps with NA, and perform locf()
  ta_q[price_jumps(mid_prices, win_dow=win_dow, vol_mult=vol_mult), ] <- NA
  na.locf(ta_q)
}  # end scrub_taq




#' Scrub a single day of \code{TAQ} data, aggregate it, and convert to
#' \code{OHLC} format.
#'
#' @export
#' @inheritParams scrub_taq
#' @param period aggregation period.
#' @return A \code{OHLC} time series in \code{xts} format.
#' @details The function \code{scrub_agg()} performs:
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
#' # create random TAQ prices
#' ta_q <- HighFreq::random_taq()
#' # aggregate to ten minutes OHLC data
#' oh_lc <- HighFreq::scrub_agg(ta_q, period="10 min")
#' chart_Series(oh_lc, name="random prices")
#' # scrub and aggregate a single day of SPY TAQ data to OHLC
#' oh_lc <- HighFreq::scrub_agg(ta_q=SPY_TAQ)
#' chart_Series(oh_lc, name=sym_bol)

scrub_agg <- function(ta_q, win_dow=51, vol_mult=2,
                      period="minutes", tzone="America/New_York") {
# convert timezone of index to New_York
  index(ta_q) <- lubridate::with_tz(time=index(ta_q), tzone=tzone)
# subset data to NYSE trading hours
  ta_q <- ta_q['T09:30:00/T16:00:00', ]
# return NULL if no data
  if (NROW(ta_q)==0)  return(NULL)
#  to_day <- as.Date(index(first(ta_q)))

# remove duplicate time stamps using duplicated()
  ta_q <- ta_q[!duplicated(index(ta_q)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- extreme_values(bid_offer, win_dow=win_dow, vol_mult=vol_mult)
# remove suspect values
  ta_q <- ta_q[!sus_pect]
# replace suspect values
# ta_q[sus_pect, "Bid.Price"] <- ta_q[sus_pect, "Trade.Price"]
# ta_q[sus_pect, "Ask.Price"] <- ta_q[sus_pect, "Trade.Price"]

# scrub quotes with suspect price jumps
# calculate mid prices
  mid_prices <- 0.5 * (ta_q[, "Bid.Price"] + ta_q[, "Ask.Price"])
#  mid_prices <- na.omit(mid_prices)
  colnames(mid_prices) <- "Mid.Price"
# replace whole rows containing suspect price jumps with NA, and perform locf()
  mid_prices[price_jumps(mid_prices, win_dow=win_dow, vol_mult=vol_mult)] <- NA
  mid_prices <- na.locf(mid_prices)
#  mid_prices <- na.locf(mid_prices, fromLast=TRUE)
# cbind mid_prices with volume data, and replace NA volumes with zero
  mid_prices <- cbind(mid_prices, ta_q[index(mid_prices), "Volume"])
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
#' @return An \code{OHLC} time series in \code{xts} format.
#' @details The function \code{save_scrub_agg()} loads multiple days of
#'   \code{TAQ} data, then scrubs, aggregates, and rbinds them into a
#'   \code{OHLC} time series, and finally saves it to a single
#'   \sQuote{\code{*.RData}} file. The \code{OHLC} time series is stored in a
#'   variable named \sQuote{\code{symbol}}, and then it's saved to a file named
#'   \sQuote{\code{symbol.RData}} in the \sQuote{\code{output_dir}} directory.
#'   The \code{TAQ} data files are assumed to be stored in separate directories
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}}
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \code{TAQ}
#'   data.
#' @examples
#' \dontrun{
#' # set data directories
#' data_dir <- "C:/Develop/data/hfreq/src/"
#' output_dir <- "C:/Develop/data/hfreq/scrub/"
#' sym_bol <- "SPY"
#' # aggregate SPY TAQ data to 15-min OHLC bar data, and save the data to a file
#' save_scrub_agg(sym_bol=sym_bol, data_dir=data_dir, output_dir=output_dir, period="15 min")
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
  sym_bol <- load(file_name)
  scrub_agg(get(sym_bol),
            win_dow=win_dow,
            vol_mult=vol_mult,
            period=period, tzone=tzone)
})  # end sapply

# recursively "rbind" the list into a single xts
  da_ta <- rutils::do_call_rbind(da_ta)
# assign column names, i.e. "symbol.High"
  colnames(da_ta) <- sapply(strsplit(colnames(da_ta), split="[.]"),
                           function(strng) paste(sym_bol, strng[-1], sep="."))

# copy the xts data to a variable with the name 'sym_bol'
  assign(sym_bol, da_ta)

# save the xts data to a file in the output_dir
  save(list=sym_bol, file=file.path(output_dir, paste0(sym_bol, ".RData")))
  invisible(sym_bol)

}  # end save_scrub_agg




#' Load and scrub multiple days of \code{TAQ} data for a single symbol, and save
#' it to multiple \sQuote{\code{*.RData}} files.
#'
#' @export
#' @inheritParams save_scrub_agg
#' @return A \code{TAQ} time series in \code{xts} format.
#' @details The function \code{save_taq()} loads multiple days of \code{TAQ}
#'   data, scrubs it, and saves the scrubbed TAQ data to individual
#'   \sQuote{\code{*.RData}} files. It uses the same file names for output as
#'   the input file names. The \code{TAQ} data files are assumed to be stored in
#'   separate directories for each \sQuote{\code{symbol}}. Each
#'   \sQuote{\code{symbol}} has its own directory (named \sQuote{\code{symbol}})
#'   in the \sQuote{\code{data_dir}} directory.
#'   Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \code{TAQ}
#'   data.
#' @examples
#' \dontrun{
#' save_taq("SPY")
#' }

save_taq <- function(sym_bol,
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
    sym_bol <- load(file_name_in)
    file_name_out <- file.path(output_dir, file_name)
# save the xts data to a file in the output_dir
    ta_q <- scrub_taq(get(sym_bol), win_dow=win_dow, vol_mult=vol_mult, tzone=tzone)
    if (!is.null(ta_q)) {
      assign(sym_bol, ta_q)
      save(list=sym_bol, file=file_name_out)
      cat("finished saving", sym_bol, "to file: ", file_name, "\n")
    }
    file_name
  })  # end sapply

  invisible(sym_bol)

}  # end save_taq




#' Load, scrub, aggregate, and rbind multiple days of \code{TAQ} data for a
#' single symbol. Calculate returns and save them to a single \sQuote{\code{*.RData}}
#' file.
#'
#' @export
#' @inheritParams save_scrub_agg
#' @return A time series of returns and volume in \code{xts} format.
#' @details The function \code{save_rets} loads multiple days of \code{TAQ}
#'   data, then scrubs, aggregates, and rbinds them into a \code{OHLC} time
#'   series.  It then calculates returns using function \code{run_returns}, and
#'   stores them in a variable named \sQuote{\code{symbol.rets}}, and saves them
#'   to a file called \sQuote{\code{symbol.rets.RData}}.
#'   The \code{TAQ} data files are assumed to be stored in separate directories
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}}
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \code{TAQ}
#'   data.
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
  ta_q <- sapply(file_names, function(file_name) {
    cat("loading", sym_bol, "from file: ", file_name, "\n")
    sym_bol <- load(file_name)
    get(sym_bol)
  })

# scrub and aggregate the TAQ data
  oh_lc <- lapply(ta_q, scrub_agg,
                      win_dow=win_dow,
                      vol_mult=vol_mult,
                      period=period,
                      tzone=tzone)

# calculate returns
  oh_lc <- lapply(oh_lc, run_returns)

# recursively "rbind" the list into a single xts
  oh_lc <- rutils::do_call_rbind(oh_lc)
# assign column names, i.e. "symbol.rets"
  colnames(oh_lc) <-
    c(paste(sym_bol, "rets", sep="."), paste(sym_bol, "vol", sep="."))

# copy the xts data to a variable with the name 'sym_bol'
  sym_bol_rets <- paste(sym_bol, "rets", sep=".")
  assign(sym_bol_rets, oh_lc)

# save the xts data to a file in the output_dir
  save(list=eval(sym_bol_rets),
       file=file.path(output_dir, paste0(sym_bol_rets, ".RData")))
  invisible(sym_bol_rets)

}  # end save_rets




#' Load \code{OHLC} time series data for a single symbol, calculate its returns,
#' and save them to a single \sQuote{\code{*.RData}} file, without aggregation.
#'
#' @export
#' @inheritParams save_scrub_agg
#' @return A time series of returns and volume in \code{xts} format.
#' @details The function \code{save_rets_ohlc()} loads \code{OHLC} time series
#'   data from a single file.  It then calculates returns using function
#'   \code{run_returns}, and stores them in a variable named
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
  sym_bol <- load(file_name)

# calculate returns
  da_ta <- run_returns(get(sym_bol))

# copy the xts data to a variable with the name 'sym_bol'
  sym_bol_rets <- paste(sym_bol, "rets", sep=".")
  assign(sym_bol_rets, da_ta)

# save the xts data to a file in the output_dir
  cat("saving", sym_bol, "to file: ", paste0(sym_bol_rets, ".RData"), "\n")
  save(list=eval(sym_bol_rets), file=file.path(output_dir, paste0(sym_bol_rets, ".RData")))
  invisible(sym_bol_rets)

}  # end save_rets_ohlc




#' Calculate a time series of variance estimates for an \code{OHLC} time series,
#' assuming zero drift.
#'
#' Calculates the variance estimates for each bar of \code{OHLC} prices at each
#' point in time (row), using the squared differences of \code{OHLC} prices at
#' each point in time.
#'
#' @export
#' @param oh_lc an \code{OHLC} time series of prices in \code{xts} format.
#' @param calc_method \code{character} string representing method for estimating
#'   variance.  The methods include:
#'   \itemize{
#'     \item "close" close to close,
#'     \item "garman_klass" Garman-Klass,
#'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
#'     \item "rogers_satchell" Rogers-Satchell,
#'     \item "yang_zhang" Yang-Zhang,
#'    }
#'    (default is \code{"yang_zhang"})
#' @return An \code{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{run_variance()} calculates a time series of
#'   variance estimates from \code{OHLC} prices, one for each bar of \code{OHLC}
#'   data.
#'
#'   The user can choose from several different variance estimation methods. The
#'   methods \code{"close"}, \code{"garman_klass_yz"}, and \code{"yang_zhang"}
#'   do account for close-to-open price jumps, while the methods
#'   \code{"garman_klass"} and \code{"rogers_satchell"} do not account for
#'   close-to-open price jumps. The default method is \code{"yang_zhang"}, which
#'   theoretically has the lowest standard error among unbiased estimators. All
#'   the methods are implemented assuming zero drift, for two reasons. First,
#'   the drift in daily or intraday data is insignificant compared to the
#'   volatility. Second, the purpose of the function \code{run_variance()} is to
#'   produce technical indicators, rather than statistical estimates.
#'
#'   The variance is scaled to the scale of the time index of the \code{OHLC} 
#'   time series.  For example, if the time index is in seconds, then the 
#'   variance is expressed in units equal to the variance per second, if the 
#'   time index is in days, then the variance is equal to the variance per day.
#'   The function \code{run_variance()} performs a similar operation to the 
#'   function \code{volatility()} from package 
#'   \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR}, but it 
#'   assumes zero drift, and doesn't calculate a running sum using 
#'   \code{runSum()}.  It's also a little faster because it performs less data 
#'   validation.
#' @examples
#' # create minutely OHLC time series of random prices
#' oh_lc <- HighFreq::random_ohlc()
#' # calculate variance estimates for oh_lc
#' var_running <- HighFreq::run_variance(oh_lc)
#' # calculate variance estimates for SPY
#' var_running <- HighFreq::run_variance(SPY, calc_method="yang_zhang")
#' # calculate SPY variance without overnight jumps
#' var_running <- HighFreq::run_variance(SPY, calc_method="rogers_satchell")

run_variance <- function(oh_lc, calc_method="yang_zhang") {
  sym_bol <- rutils::na_me(oh_lc)
  oh_lc <- log(oh_lc[, 1:4])
  vari_ance <- switch(calc_method,
         "close"={rutils::diff_xts(oh_lc[, 4])^2},
         "garman_klass"={0.5*(oh_lc[, 2]-oh_lc[, 3])^2 -
                         (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^2},
         "rogers_satchell"={(oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                            (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1])},
         "garman_klass_yz"={(oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]))^2 +
                            0.5*(oh_lc[, 2]-oh_lc[, 3])^2 -
                            (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^2},
         "yang_zhang"={(oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]))^2 +
                       0.67*(oh_lc[, 1]-oh_lc[, 4])^2 +
                       0.33*((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                               (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]))}
  )  # end switch
  vari_ance <- vari_ance/c(1, diff(.index(oh_lc)))^2
  vari_ance[1, ] <- 0
  vari_ance <- na.locf(vari_ance)
  colnames(vari_ance) <- paste0(sym_bol, ".Variance")
  vari_ance
}  # end run_variance




#' Calculate time series of skew estimates from a \code{OHLC} time series,
#' assuming zero drift.
#'
#' @export
#' @param oh_lc an \code{OHLC} time series of prices in \code{xts} format.
#' @param calc_method \code{character} string representing method for estimating
#'   skew.
#' @return A time series of skew estimates.
#' @details The function \code{run_skew()} calculates a time series of skew
#'   estimates from \code{OHLC} prices, one for each bar of \code{OHLC} data.
#'   The skew estimates are scaled to the time scale of the index of the
#'   \code{OHLC} time series.  For example, if the time index is in seconds,
#'   then the estimates are equal to the skew per second, if the time index is
#'   in days, then the estimates are equal to the skew per day.
#'   Currently only the \code{"close"} skew estimation method is correct
#'   (assuming zero drift), while the \code{"rogers_satchell"} method produces a
#'   skew-like indicator, proportional to the skew. The default method is
#'   \code{"rogers_satchell"}.
#' @examples
#' # calculate time series of skew estimates for SPY
#' sk_ew <- HighFreq::run_skew(SPY)

run_skew <- function(oh_lc, calc_method="rogers_satchell") {
  sym_bol <- rutils::na_me(oh_lc)
  oh_lc <- log(oh_lc[, 1:4])
  sk_ew <- switch(calc_method,
                  "close"={rutils::diff_xts(oh_lc[, 4])^3},
                  "garman_klass"={0.5*(oh_lc[, 2]-oh_lc[, 3])^3 -
                      (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^3},
                  "rogers_satchell"={
                    (oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1])*(oh_lc[, 2]-0.5*(oh_lc[, 4] + oh_lc[, 1])) +
                      (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1])*(oh_lc[, 3]-0.5*(oh_lc[, 4] + oh_lc[, 1]))},
                  "garman_klass_yz"={(oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]))^3 +
                      0.5*(oh_lc[, 2]-oh_lc[, 3])^3 -
                      (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^3},
                  "yang_zhang"={c_o <- oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]);
                  o_c <- oh_lc[, 1]-oh_lc[, 4];
                  (c_o-sum(c_o)/NROW(c_o))^3 +
                    0.67*(o_c-sum(o_c)/NROW(o_c))^3 +
                    0.33*((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                            (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]))}
  )  # end switch
  sk_ew <- sk_ew/c(1, diff(.index(oh_lc)))^3
  sk_ew[1, ] <- 0
  sk_ew <- na.locf(sk_ew)
  colnames(sk_ew) <- paste0(sym_bol, ".Skew")
  sk_ew
}  # end run_skew




#' Calculate time series of Sharpe-like statistics for each bar of a \code{OHLC}
#' time series.
#'
#' @export
#' @param oh_lc an \code{OHLC} time series of prices in \code{xts} format.
#' @param calc_method \code{character} string representing method for estimating
#'   the Sharpe-like exponent.
#' @return An \code{xts} time series with the same number of rows as the
#'   argument \code{oh_lc}.
#' @details The function \code{run_sharpe()} calculates Sharpe-like statistics
#'   for each bar of a \code{OHLC} time series.
#'   The Sharpe-like statistic is defined as the ratio of the difference between
#'   \code{Close} minus \code{Open} prices divided by the difference between
#'   \code{High} minus \code{Low} prices.
#'   This statistic may also be interpreted as something like a Hurst exponent
#'   for a single bar of data.
#'   The motivation for the Sharpe-like statistic is the notion that if prices
#'   are trending in the same direction inside a given time bar of data, then
#'   this statistic is close to either 1 or -1.
#' @examples
#' # calculate time series of running Sharpe ratios for SPY
#' sharpe_running <- run_sharpe(SPY)

run_sharpe <- function(oh_lc, calc_method="close") {
  sharpe_ratio <- switch(calc_method,
                   "close"={(oh_lc[, 4]-oh_lc[, 1])/(oh_lc[, 2]-oh_lc[, 3])},
                   "method2"={(oh_lc[, 4]-oh_lc[, 1])/(oh_lc[, 2]-oh_lc[, 3])}
  )  # end switch
  sharpe_ratio <- ifelse(oh_lc[, 2]==oh_lc[, 3], 0, sharpe_ratio)
  colnames(sharpe_ratio) <- paste0(rutils::na_me(oh_lc), ".Sharpe")
  sharpe_ratio
}  # end run_sharpe




#' Calculate the aggregation (weighted average) of a statistical estimator over
#' a \code{OHLC} time series.
#'
#' @export
#' @param oh_lc \code{OHLC} time series of prices and trading volumes, in
#'   \code{xts} format.
#' @param mo_ment \code{character} string representing function for
#'   estimating the moment.
#' @param weight_ed \code{Boolean} should estimate be weighted by the trading
#'   volume? (default is \code{TRUE})
#' @param ... additional parameters to the mo_ment function.
#' @return A single \code{numeric} value equal to the volume weighted average of
#'   an estimator over the time series.
#' @details The function \code{agg_regate()} calculates a single number
#'   representing the volume weighted average of an estimator over the
#'   \code{OHLC} time series of prices.  By default the sum is trade volume
#'   weighted.
#' @examples
#' # calculate weighted average variance for SPY (single number)
#' vari_ance <- agg_regate(oh_lc=SPY, mo_ment="run_variance")
#' # calculate time series of daily skew estimates for SPY
#' skew_daily <- apply.daily(x=SPY, FUN=agg_regate, mo_ment="run_skew")

agg_regate <- function(oh_lc, mo_ment="run_variance", weight_ed=TRUE, ...) {
# match "mo_ment" with moment function
  mo_ment <- match.fun(mo_ment)
  agg_regations <- mo_ment(oh_lc, ...)
# weight the estimates by volume
  if (weight_ed) {
    agg_regations <- oh_lc[, 5]*agg_regations
    agg_regations <- sum(agg_regations)/sum(oh_lc[, 5])
  } else
    agg_regations <- sum(agg_regations)/length(agg_regations)
  agg_regations
}  # end agg_regate




#' Calculate the volume-weighted average price of an \code{OHLC} time series
#' over a rolling window (lookback period).
#'
#' Performs the same operation as function \code{VWAP()} from package
#' \href{https://cran.r-project.org/web/packages/TTR/index.html}{VWAP},
#' but using vectorized functions, so it's a little faster.
#'
#' @export
#' @param oh_lc an \code{OHLC} time series of prices in \code{xts} format.
#' @param x_ts single-column \code{xts} time series.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for calculating the average price.
#' @return An \code{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{roll_vwap()} calculates the volume-weighted
#'   average closing price, defined as the sum of the prices multiplied by
#'   trading volumes in the lookback window, divided by the sum of trading
#'   volumes in the window. If the argument \code{x_ts} is passed in explicitly,
#'   then its volume-weighted average value over time is calculated.
#' @examples
#' # calculate and plot rolling volume-weighted average closing prices (VWAP)
#' prices_rolling <- roll_vwap(oh_lc=SPY["2013-11"], win_dow=11)
#' chart_Series(SPY["2013-11-12"], name="SPY prices")
#' add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)
#' legend("top", legend=c("SPY prices", "VWAP prices"),
#' bg="white", lty=c(1, 1), lwd=c(2, 2),
#' col=c("black", "red"), bty="n")
#' # calculate running returns
#' returns_running <- run_returns(x_ts=SPY)
#' # calculate the rolling volume-weighted average returns
#' roll_vwap(oh_lc=SPY, x_ts=returns_running, win_dow=11)

roll_vwap <- function(oh_lc, x_ts=oh_lc[, 4], win_dow) {
  roll_vwap <- rutils::roll_sum(x_ts=x_ts*oh_lc[, 5], win_dow=win_dow)
  volume_rolling <- rutils::roll_sum(x_ts=oh_lc[, 5], win_dow=win_dow)
  roll_vwap <- roll_vwap/volume_rolling
  roll_vwap[is.na(roll_vwap)] <- 0
  colnames(roll_vwap) <- paste0(rutils::na_me(oh_lc), ".VWAP")
  roll_vwap
}  # end roll_vwap




#' Calculate a vector of statistics over an \code{OHLC} time series, and
#' calculate a rolling mean over the statistics.
#'
#' @export
#' @param oh_lc \code{OHLC} time series of prices and trading volumes, in
#'   \code{xts} format.
#' @param mo_ment \code{character} string representing a function for
#'   estimating statistics of a single bar of \code{OHLC} data, such as
#'   volatility, skew, and higher moments.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for calculating the rolling mean.
#' @param weight_ed \code{Boolean} should statistic be weighted by trade volume?
#'   (default \code{TRUE})
#' @param ... additional parameters to the mo_ment function.
#' @return An \code{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{roll_moment()} calculates a vector of statistics
#'   over an \code{OHLC} time series, such as volatility, skew, and higher
#'   moments.  The statistics could also be any other aggregation of a single
#'   bar of \code{OHLC} data, for example the \code{High} price minus the
#'   \code{Low} price squared.  The length of the vector of statistics is equal
#'   to the number of rows of the argument \code{oh_lc}. Then it calculates a
#'   trade volume weighted rolling mean over the vector of statistics over and
#'   calculate statistics.
#' @examples
#' # calculate time series of rolling variance and skew estimates
#' var_rolling <- roll_moment(oh_lc=SPY, win_dow=21)
#' skew_rolling <- roll_moment(oh_lc=SPY, mo_ment="run_skew", win_dow=21)
#' skew_rolling <- skew_rolling/(var_rolling)^(1.5)
#' skew_rolling[1, ] <- 0
#' skew_rolling <- na.locf(skew_rolling)

roll_moment <- function(oh_lc, mo_ment="run_variance", win_dow=11, weight_ed=TRUE, ...) {
# match "mo_ment" with moment function
  mo_ment <- match.fun(mo_ment)
  agg_regations <- mo_ment(oh_lc, ...)
# weight by volume
  if (weight_ed) {
    agg_regations <- oh_lc[, 5]*agg_regations
    volume_rolling <- rutils::roll_sum(oh_lc[, 5], win_dow=win_dow)
    agg_regations <- rutils::roll_sum(agg_regations, win_dow=win_dow)/volume_rolling
    agg_regations[is.na(agg_regations)] <- 0
  } else
    agg_regations <- rutils::roll_sum(agg_regations, win_dow=win_dow)/win_dow
  colnames(agg_regations) <- paste(rutils::na_me(oh_lc), "Vol", sep=".")
  agg_regations
}  # end roll_moment




#' Calculate the rolling Sharpe ratio over a rolling lookback window for an
#' \code{OHLC} time series.
#'
#' @export
#' @param oh_lc an \code{OHLC} time series of prices in \code{xts} format.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for aggregating the \code{OHLC} prices.
#' @return An \code{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{roll_sharpe()} calculates the rolling Sharpe
#'   ratio as the ratio of absolute returns over the lookback window (not
#'   percentage returns), divided by the average volatility of returns.
#' @examples
#' # calculate rolling Sharpe ratio over SPY
#' sharpe_rolling <- roll_sharpe(oh_lc=SPY, win_dow=10)

roll_sharpe <- function(oh_lc, win_dow=11) {
  var_ohlc_agg <- (oh_lc[, 4] - rutils::lag_xts(oh_lc[, 1], k=(win_dow-1)))/win_dow
  var_ohlc <- sqrt(rutils::roll_sum(run_variance(oh_lc), win_dow=win_dow)/win_dow)
  sharpe_rolling <- ifelse(var_ohlc==0,
                       1.0,
                       var_ohlc_agg/var_ohlc)
  colnames(sharpe_rolling) <- paste0(rutils::na_me(oh_lc), ".Sharpe")
  na.locf(sharpe_rolling)
}  # end roll_sharpe




#' Calculate the rolling Hurst exponent over a rolling lookback window or the
#' end points of an \code{OHLC} time series.
#'
#' @export
#' @param oh_lc an \code{OHLC} time series of prices in \code{xts} format.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for aggregating the \code{OHLC} prices.
#' @param off_set the number of bars of data in the first, stub window.
#' @param roll_end_points \code{Boolean} should the Hurst exponent be calculated
#'   using aggregations over the end points, or by rolling over a lookback
#'   window? (default is \code{FALSE})
#' @return An \code{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{roll_hurst()} calculates the rolling Hurst
#'   exponent in two different ways, depending on the argument
#'   \code{roll_end_points}.
#'
#'   If \code{roll_end_points} is \code{FALSE} (the default), then the rolling
#'   Hurst exponent is calculated as the logarithm of the ratios of two rolling
#'   price range estimates. The Hurst exponent is defined as the logarithm of
#'   the ratio of the range of aggregated prices, divided by the average range
#'   of prices in each bar. The aggregated prices are calculated over
#'   overlapping windows, and the Hurst exponent values are calculated at each
#'   point in time.
#'
#'   If \code{roll_end_points} is \code{TRUE}, then the rolling Hurst exponent is
#'   calculated as the logarithm of the ratios of two rolling variance
#'   estimates. The Hurst exponent is defined as the logarithm of the ratio of
#'   the variance of aggregated returns, divided by the variance of simple
#'   returns. The aggregated returns are calculated over non-overlapping windows
#'   spanned by the end points, using the function \code{to_period()}. The Hurst
#'   exponent values are calculated only at the end points. The non-overlapping
#'   aggregation windows can be shifted by using the argument \code{off_set},
#'   which produces a slightly different series of rolling hurst exponent
#'   values.
#'
#' @examples
#' # calculate rolling Hurst over SPY
#' hurst_rolling <- roll_hurst(oh_lc=SPY, win_dow=10)
#' # calculate Hurst over end points of SPY
#' hurst_rolling <- roll_hurst(oh_lc=SPY, win_dow=10, off_set=0, roll_end_points=TRUE)
#' # calculate a series of rolling hurst values using argument off_set
#' hurst_rolling <- lapply(0:9, roll_hurst, oh_lc=SPY, win_dow=10, roll_end_points=TRUE)
#' hurst_rolling <- rutils::do_call_rbind(hurst_rolling)
#' # remove daily warmup periods
#' hurst_rolling <- hurst_rolling["T09:41:00/T16:00:00"]
#' chart_Series(x=hurst_rolling["2012-02-13"],
#'   name=paste(colnames(hurst_rolling), "10-minute aggregations"))

roll_hurst <- function(oh_lc, win_dow=11, off_set=0, roll_end_points=FALSE) {
  if(roll_end_points) {  # roll over end points
    # aggregate oh_lc to a lower periodicity specified by end_points
    ohlc_agg <- rutils::to_period(oh_lc=oh_lc,
                                  end_points=rutils::end_points(oh_lc, inter_val=win_dow, off_set=off_set))
    var_ohlc_agg <- run_variance(ohlc_agg)
    in_dex <- index(ohlc_agg)
    var_ohlc <- rutils::roll_sum(run_variance(oh_lc), win_dow=win_dow)[in_dex]/win_dow
  }
  else {  # roll over overlapping windows
    max_hi <- TTR::runMax(x=oh_lc[, 2], n=win_dow)
    min_lo <- -TTR::runMax(x=-oh_lc[, 3], n=win_dow)
    var_ohlc_agg <- max_hi - min_lo
    var_ohlc_agg[1:(win_dow-1), ] <- 0
    var_ohlc <- rutils::roll_sum((oh_lc[, 2] - oh_lc[, 3]), win_dow=win_dow)/win_dow
  }  # end if
  hurst_rolling <- ifelse((var_ohlc==0) | (var_ohlc_agg==0),
                       1.0,
                       log(var_ohlc_agg/var_ohlc)/log(win_dow))
  colnames(hurst_rolling) <- paste0(rutils::na_me(oh_lc), ".Hurst")
  na.locf(hurst_rolling)
}  # end roll_hurst




#' Apply an aggregation function over a rolling lookback window and the end
#' points of an \code{OHLC} time series.
#'
#' @export
#' @param oh_lc \code{OHLC} time series of prices and trading volumes, in
#'   \code{xts} format.
#' @param agg_fun \code{character} string representing an aggregation function
#'   to be applied over a rolling lookback window.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for applying the aggregation function.
#' @param by_columns \code{Boolean} should the function \code{agg_fun()} be
#'   applied column-wise (individually), or should it be applied to all the
#'   columns combined? (default is \code{FALSE})
#' @param end_points an integer vector of end points.
#' @param ... additional parameters to the agg_fun function.
#' @return An \code{xts} time series with the same number of rows as the
#'   argument \code{oh_lc}.
#' @details The function \code{roll_apply()} applies an aggregation function
#'   over a rolling lookback window and the end points of an \code{OHLC} time
#'   series.
#'
#'   Performs similar operations to the functions \code{rollapply()} and
#'   \code{period.apply()} from package
#'   \href{https://cran.r-project.org/web/packages/xts/index.html}{xts}, and
#'   also the function \code{apply.rolling()} from package
#'   \href{https://cran.r-project.org/web/packages/PerformanceAnalytics/index.html}{PerformanceAnalytics}.
#'   (The function \code{rollapply()} isn't exported from the package
#'   \href{https://cran.r-project.org/web/packages/xts/index.html}{xts}.)
#'
#'   But the function \code{roll_apply()} is faster because it performs less
#'   type-checking and other overhead. Unlike the other functions,
#'   \code{roll_apply()} doesn't produce any leading \code{NA} values.
#'
#'   The function \code{roll_apply()} can be called in two different ways,
#'   depending on the argument \code{end_points}.
#'   If the argument \code{end_points} isn't explicitly passed to
#'   \code{roll_apply()}, then the default value is used, and
#'   \code{roll_apply()} performs aggregations over overlapping windows at each
#'   point in time.
#'   If the argument \code{end_points} is explicitly passed to
#'   \code{roll_apply()}, then \code{roll_apply()} performs aggregations over
#'   overlapping windows spanned by the end_points.
#'
#'   The aggregation function \code{agg_fun} can return either a single value or
#'   a vector of values. If the aggregation function \code{agg_fun} returns a
#'   single value, then \code{roll_apply()} returns an \code{xts} time series
#'   with a single column. If the aggregation function \code{agg_fun} returns a
#'   vector of values, then \code{roll_apply()} returns an \code{xts} time
#'   series with multiple columns equal to the length of the vector returned by
#'   the aggregation function \code{agg_fun}.
#'
#' @examples
#' # extract a single day of SPY data
#' x_ts <- SPY["2012-02-13"]
#' win_dow <- 11
#' # calculate the rolling sums of the columns of x_ts
#' agg_regations <- roll_apply(x_ts, agg_fun=sum, win_dow=win_dow, by_columns=TRUE)
#' # apply a vector-valued aggregation function over a rolling window
#' agg_function <- function(x_ts)  c(max(x_ts[, 2]), min(x_ts[, 3]))
#' agg_regations <- roll_apply(x_ts, agg_fun=agg_function, win_dow=win_dow)
#' # define end points at 11-minute intervals (SPY is minutely bars)
#' end_points <- rutils::end_points(x_ts, inter_val=win_dow)
#' # calculate the rolling sums of the columns of x_ts over end_points
#' agg_regations <- roll_apply(x_ts, agg_fun=sum, win_dow=2, end_points=end_points, by_columns=TRUE)
#' # apply a vector-valued aggregation function over the end_points of x_ts
#' agg_regations <- roll_apply(x_ts, agg_fun=agg_function, win_dow=2, end_points=end_points)

roll_apply <- function(oh_lc, agg_fun="run_variance", win_dow=11,
                       end_points=(0:NROW(oh_lc)), by_columns=FALSE, ...) {
  # match "agg_fun" with some aggregation function
  agg_fun <- match.fun(agg_fun)
  len_gth <- length(end_points)
  # define start_points as lag of end_points
  start_points <-  end_points[c(rep_len(1, win_dow-1), 1:(len_gth-win_dow+1))] +
    (NROW(oh_lc) > (len_gth+1))
  # perform aggregations over length of end_points
  agg_regations <- if(by_columns)
    sapply(oh_lc, function(col_umn)
      sapply(2:len_gth, function(in_dex)
        agg_fun(.subset_xts(col_umn,
                            start_points[in_dex]:end_points[in_dex]), ...)
      ))  # end sapply
  else {  # not by_columns
    agg_regations <- sapply(2:len_gth, function(in_dex)
      agg_fun(.subset_xts(oh_lc,
                          start_points[in_dex]:end_points[in_dex]), ...)
    )  # end sapply
    # coerce agg_regations into matrix and transpose it
    if (is.vector(agg_regations))
      agg_regations <- t(agg_regations)
    agg_regations <- t(agg_regations)
  }  # end if
  # coerce agg_regations into xts series
  xts(agg_regations, order.by=index(oh_lc[end_points]))
}  # end roll_apply




#' Perform seasonality aggregations over a single-column \code{xts} time series.
#'
#' @export
#' @param x_ts single-column \code{xts} time series.
#' @param in_dex vector of \code{character} strings representing points in time,
#'   of the same length as the argument \code{x_ts}.
#' @return An \code{xts} time series with mean aggregations over the seasonality
#'   interval.
#' @details The function \code{season_ality()} calculates the mean of values
#'   observed at the same points in time specified by the argument
#'   \code{in_dex}. An example of a daily seasonality aggregation is the average
#'   price of a stock between 9:30AM and 10:00AM every day, over many days. The
#'   argument \code{in_dex} is passed into function \code{tapply()}, and must be
#'   the same length as the argument \code{x_ts}.
#' @examples
#' # calculate running variance of each minutely OHLC bar of data
#' x_ts <- run_variance(SPY)
#' # remove overnight variance spikes at "09:31"
#' in_dex <- format(index(x_ts), "%H:%M")
#' x_ts <- x_ts[!in_dex=="09:31", ]
#' # calculate daily seasonality of variance
#' var_seasonal <- season_ality(x_ts=x_ts)
#' chart_Series(x=var_seasonal, name=paste(colnames(var_seasonal),
#'   "daily seasonality of variance"))

season_ality <- function(x_ts, in_dex=format(index(x_ts), "%H:%M")) {
# aggregate the mean
  agg_regation <- tapply(X=x_ts, INDEX=in_dex, FUN=mean)
# coerce from array to named vector
  agg_regation <- structure(as.vector(agg_regation), names=names(agg_regation))
# coerce to xts
  agg_regation <- xts(x=agg_regation,
      order.by=as.POSIXct(paste(Sys.Date(), names(agg_regation))))
  colnames(agg_regation) <- colnames(x_ts)
  agg_regation
}  # end season_ality


