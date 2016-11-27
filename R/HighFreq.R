#' Calculate a random \emph{TAQ} time series of prices and trading volumes, in
#' \emph{xts} format.
#'
#' Calculate a \emph{TAQ} time series of random prices following geometric
#' Brownian motion, combined with random trading volumes.
#'
#' @export
#' @param vol_at volatility per period of the \code{in_dex} time index (default
#'   is \code{6.5e-05} per second, or about \code{0.01=1.0\%} per day).
#' @param dri_ft drift per period of the \code{in_dex} time index (default is
#'   0.0).
#' @param in_dex time index for the \emph{TAQ} time series.
#' @param bid_offer the bid-offer spread expressed as a fraction of the prices
#'   (default is 0.001=10bps).
#' @return An \emph{xts} time series, with time index equal to the input
#'   \code{in_dex} time index, and with four columns containing the bid, ask,
#'   and trade prices, and the trade volume.
#' @details The function \code{random_taq()} calculates an \emph{xts} time 
#'   series with four columns containing random prices following geometric 
#'   Brownian motion: the bid, ask, and trade prices, combined with random trade
#'   volume data.
#'   If \code{in_dex} isn't supplied as an argument, then by default it's
#'   equal to the secondly index over the two previous calendar days.
#' @examples
#' # create secondly TAQ time series of random prices
#' ta_q <- HighFreq::random_taq()
#' # create random TAQ time series from SPY index
#' ta_q <- HighFreq::random_taq(in_dex=index(SPY["2012-02-13/2012-02-15"]))

random_taq <- function(vol_at=6.5e-5, dri_ft=0.0, 
  in_dex=seq(from=as.POSIXct(paste(Sys.Date()-3, "09:30:00")),
             to=as.POSIXct(paste(Sys.Date()-1, "16:00:00")), by="1 sec"),
  bid_offer=0.001, ...) {
  len_gth <- NROW(in_dex)
  # create xts of random prices following geometric Brownian motion
  ta_q <- xts(exp(cumsum(vol_at*rnorm(len_gth) + dri_ft - vol_at^2/2)), 
              order.by=in_dex)
  # create vector of random bid-offer spreads
  bid_offer <- bid_offer*(1 + runif(len_gth))/2
  # create TAQ data from bid and offer prices
  ta_q <- merge(ta_q*(1-bid_offer), ta_q*(1+bid_offer))
  # add traded price to TAQ data
  r_unif <- runif(len_gth)
  ta_q <- merge(ta_q, r_unif*ta_q[, 1] + (1-r_unif)*ta_q[, 2])
  # add trade volume column
  ta_q <- merge(ta_q, sample(x=10*(2:18), size=len_gth, replace=TRUE))
  colnames(ta_q) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
  ta_q
}  # end random_taq




#' Calculate a random \emph{OHLC} time series of prices and trading volumes, in
#' \emph{xts} format.
#'
#' Calculate a random \emph{OHLC} time series either by simulating random prices
#' following geometric Brownian motion, or by randomly sampling from an input 
#' time series.
#'
#' @export
#' @param oh_lc \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format (default is \emph{NULL}).
#' @param vol_at volatility per period of the \code{in_dex} time index (default
#'   is \code{6.5e-05} per second, or about \code{0.01=1.0\%} per day).
#' @param dri_ft drift per period of the \code{in_dex} time index (default is
#'   0.0).
#' @param in_dex time index for the \emph{OHLC} time series.
#' @param re_duce \emph{Boolean} argument: should \code{oh_lc} time series be
#'   transformed to reduced form? (default is \code{TRUE})
#' @return An \emph{xts} time series with the same dimensions and the same time
#'   index as the input \code{oh_lc} time series.
#' @details If the input \code{oh_lc} time series is \emph{NULL} (the default), 
#'   then the function \code{random_ohlc()} simulates a minutely \emph{OHLC} 
#'   time series of random prices following geometric Brownian motion, over the
#'   two previous calendar days.
#'   
#'   If the input \code{oh_lc} time series is not \emph{NULL}, then the rows of
#'   \code{oh_lc} are randomly sampled, to produce a random time series.
#'   
#'   If \code{re_duce} is \code{TRUE} (the default), then the \code{oh_lc} time
#'   series is first transformed to reduced form, then randomly sampled, and
#'   finally converted to standard form.
#'   
#'   Note: randomly sampling from an intraday time series over multiple days
#'   will cause the overnight price jumps to be re-arranged into intraday price
#'   jumps.  This will cause moment estimates to become inflated compared to the
#'   original time series.
#' @examples
#' # create minutely synthetic OHLC time series of random prices
#' oh_lc <- HighFreq::random_ohlc()
#' # create random time series from SPY by randomly sampling it
#' oh_lc <- HighFreq::random_ohlc(oh_lc=SPY["2012-02-13/2012-02-15"])

random_ohlc <- function(oh_lc=NULL, re_duce=TRUE, vol_at=6.5e-5, dri_ft=0.0, 
    in_dex=seq(from=as.POSIXct(paste(Sys.Date()-3, "09:30:00")), 
      to=as.POSIXct(paste(Sys.Date()-1, "16:00:00")), by="1 sec"), ...) {
  if (is.null(oh_lc)) {
    len_gth <- NROW(in_dex)
    # create xts of random prices following geometric Brownian motion
    x_ts <- xts(exp(cumsum(vol_at*rnorm(len_gth) + dri_ft - vol_at^2/2)), order.by=in_dex)
    # add trade volume column
    x_ts <- merge(x_ts, volume=sample(x=10*(2:18), size=len_gth, replace=TRUE))
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




#' Remove overnight close-to-open price jumps from an \emph{OHLC} time series, 
#' by adding adjustment terms to its prices.
#'
#' @export
#' @param oh_lc \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#' @return An \emph{OHLC} time series with the same dimensions and the same time
#'   index as the input \code{oh_lc} time series.
#' @details The function \code{remove_jumps()} removes the overnight 
#'   close-to-open price jumps from an \emph{OHLC} time series, by adjusting its
#'   prices so that the first \emph{Open} price of the day is equal to the last 
#'   \emph{Close} price of the previous day.
#'   
#'   The function \code{remove_jumps()} adds adjustment terms to all the 
#'   \emph{OHLC} prices, so that intra-day returns and volatilities are not
#'   affected.
#' 
#'   The function \code{remove_jumps()} identifies overnight periods as those 
#'   that are greater than 60 seconds. This assumes that intra-day periods
#'   between neighboring bars of data are 60 seconds or less.
#'   
#'   The time index of the \code{oh_lc} time series is assumed to be in 
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'   
#' @examples
#' # remove overnight close-to-open price jumps from SPY data
#' oh_lc <- remove_jumps(SPY)

remove_jumps <- function(oh_lc) {
  # find time index of the periods greater than 60 seconds
  peri_ods <- rutils::diff_it(as.vector(xts::.index(oh_lc)))
  which_periods <- which(peri_ods > 60)
  # calculate cumulative sum of overnight price jumps
  jump_s <- numeric(NROW(oh_lc))
  jump_s[which_periods] <- as.vector(oh_lc[which_periods, 1]) - as.vector(oh_lc[which_periods-1, 4])
  jump_s <- cumsum(jump_s)
  # subtract overnight price jumps from OHLC
  oh_lc[, 1:4] <- coredata(oh_lc[, 1:4]) - jump_s
  oh_lc
}  # end remove_jumps




#' Calculate single period percentage returns from either \emph{TAQ} or
#' \emph{OHLC} prices.
#'
#' @export
#' @param x_ts \emph{xts} time series of either \emph{TAQ} or \emph{OHLC} data.
#' @param lag integer equal to the number of time periods of lag. (default is 1)
#' @param col_umn the column number to extract from the \emph{OHLC} data.
#'   (default is \code{4}, or the \emph{Close} prices column)
#' @param sca_le \emph{Boolean} argument: should the returns be divided by the 
#'   number of seconds in each period? (default is \code{TRUE})
#' @return A single-column \emph{xts} time series of returns.
#' @details The function \code{run_returns()} calculates the percentage returns 
#'   for either \emph{TAQ} or \emph{OHLC} data, defined as the difference of log
#'   prices.  Multi-period returns can be calculated by setting the \code{lag} 
#'   parameter to values greater than \code{1} (the default).
#'   
#'   If \code{sca_le} is \code{TRUE} (the default), then the returns are divided
#'   by the differences of the time index (which scales the returns to units of
#'   returns per second.)
#'   
#'   The time index of the \code{x_ts} time series is assumed to be in 
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'   
#'   If \code{sca_le} is \code{TRUE} (the default), then the returns are
#'   expressed in the scale of the time index of the \code{x_ts} time series. 
#'   For example, if the time index is in seconds, then the returns are given in
#'   units of returns per second.  If the time index is in days, then the
#'   returns are equal to the returns per day.
#'   
#'   The function \code{run_returns()} identifies the \code{x_ts} time series as
#'   \emph{TAQ} data when it has six columns, otherwise assumes it's \emph{OHLC}
#'   data. By default, for \emph{OHLC} data, it differences the \emph{Close}
#'   prices, but can also difference other prices depending on the value of 
#'   \code{col_umn}.
#' @examples
#' # calculate secondly returns from TAQ data
#' re_turns <- HighFreq::run_returns(x_ts=SPY_TAQ)
#' # calculate close to close returns
#' re_turns <- HighFreq::run_returns(x_ts=SPY)
#' # calculate open to open returns
#' re_turns <- HighFreq::run_returns(x_ts=SPY, col_umn=1)

run_returns <- function(x_ts, lag=1, col_umn=4, sca_le=TRUE) {
  # return NULL if no data
  if (is.null(x_ts))  return(NULL)
  # calculate mid prices
  if(NCOL(x_ts)==6)  # TAQ data has 6 columns
    re_turns <- 0.5 * (x_ts[, "Bid.Price"] + x_ts[, "Ask.Price"])
  else
    re_turns <- x_ts[, col_umn]  # OHLC data
  # calculate returns
  re_turns <- rutils::diff_xts(log(re_turns), lag=lag)
  if (sca_le)
    re_turns <- re_turns / c(rep(1, lag), diff(xts::.index(re_turns), lag=lag))
  re_turns[1:lag, ] <- 0
  colnames(re_turns) <- paste0(rutils::na_me(x_ts), ".returns")
  re_turns
}  # end run_returns




#' Calculate a \emph{Boolean} vector that identifies extreme tail values in a
#' single-column \emph{xts} time series or vector, over a rolling window.
#'
#' @export
#' @param x_ts A single-column \emph{xts} time series, or a \emph{numeric} or
#'   \emph{Boolean} vector.
#' @param win_dow number of data points for estimating rolling quantile.
#' @param vol_mult quantile multiplier.
#' @return A \emph{Boolean} vector with the same number of rows as the input
#'   time series or vector.
#' @details The function \code{which_extreme()} calculates a \emph{Boolean} 
#'   vector, with \code{TRUE} for values that belong to the extreme tails
#'   of the distribution of values.
#'   
#'   The function \code{which_extreme()} applies a version of the Hampel median 
#'   filter to identify extreme values, but instead of using the median absolute
#'   deviation (MAD), it uses the \code{0.9} quantile values calculated over a
#'   rolling window.
#'   
#'   Extreme values are defined as those that exceed the product of the 
#'   multiplier times the rolling quantile. Extreme values belong to the fat 
#'   tails of the recent (trailing) distribution of values, so they are present 
#'   only when the trailing distribution of values has fat tails.  If the
#'   trailing distribution of values is closer to normal (without fat tails),
#'   then there are no extreme values.
#'   
#'   The quantile multiplier \code{vol_mult} controls the threshold at which
#'   values are identified as extreme. Smaller quantile multiplier values will
#'   cause more values to be identified as extreme.
#'   
#' @examples
#' # create local copy of SPY TAQ data
#' ta_q <- SPY_TAQ
#' # scrub quotes with suspect bid-offer spreads
#' bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#' sus_pect <- which_extreme(bid_offer, win_dow=win_dow, vol_mult=vol_mult)
#' # remove suspect values
#' ta_q <- ta_q[!sus_pect]

which_extreme <- function(x_ts, win_dow=51, vol_mult=2) {
# calculate volatility as rolling quantile
  quan_tile <- caTools::runquantile(x=abs(as.vector(x_ts)), k=win_dow,
                        probs=0.9, endrule="constant", align="center")
#  quan_tile <- xts(quan_tile, order.by=index(x_ts))
#  colnames(quan_tile) <- "volat"
# carry forward non-zero volatility values
  quan_tile[quan_tile==0] <- NA
  quan_tile[1] <- 1
  quan_tile <- rutils::na_locf(quan_tile)
#  quan_tile <- rutils::na_locf(quan_tile, fromLast=TRUE)

# extreme value if x_ts greater than scaled volatility
  ex_treme <- (abs(x_ts) > 2*vol_mult*quan_tile)
  ex_treme[1] <- FALSE
#  colnames(ex_treme) <- "suspect"

  cat("date:", format(as.Date(index(first(x_ts)))), "\tfound", sum(ex_treme), "extreme values\n")
  ex_treme
}  # end which_extreme




#' Calculate a \emph{Boolean} vector that identifies isolated jumps (spikes) in
#' a single-column \emph{xts} time series or vector, over a rolling window.
#'
#' @export
#' @inheritParams which_extreme
#' @return A \emph{Boolean} vector with the same number of rows as the input
#'   time series or vector.
#' @details The function \code{which_jumps()} calculates a \emph{Boolean}
#'   vector, with \code{TRUE} for values that are isolated jumps (spikes).
#'   
#'   The function \code{which_jumps()} applies a version of the Hampel median 
#'   filter to identify jumps, but instead of using the median absolute 
#'   deviation (MAD), it uses the \code{0.9} quantile of returns calculated over
#'   a rolling window.
#'   This is in contrast to function \code{which_extreme()}, which applies a
#'   Hampel filter to the values themselves, instead of the returns.
#'   Returns are defined as simple differences between neighboring values.
#'   
#'   Jumps (or spikes), are defined as isolated values that are very different
#'   from the neighboring values, either before or after.  Jumps create
#'   pairs of large neighboring returns of opposite sign.
#'   
#'   Jumps (spikes) must satisfy two conditions:
#'   \enumerate{
#'     \item Neighboring returns both exceed a multiple of the rolling quantile,
#'     \item The sum of neighboring returns doesn't exceed that multiple.
#'   }
#'   
#'   The quantile multiplier \code{vol_mult} controls the threshold at which
#'   values are identified as jumps. Smaller quantile multiplier values will
#'   cause more values to be identified as jumps.
#'   
#' @examples
#' # create local copy of SPY TAQ data
#' ta_q <- SPY_TAQ
#' # calculate mid prices
#' mid_prices <- 0.5 * (ta_q[, "Bid.Price"] + ta_q[, "Ask.Price"])
#' # replace whole rows containing suspect price jumps with NA, and perform locf()
#' ta_q[which_jumps(mid_prices, win_dow=31, vol_mult=1.0), ] <- NA
#' ta_q <- zoo::na.locf(ta_q)

which_jumps <- function(x_ts, win_dow=51, vol_mult=2) {
# calculate simple returns
  re_turns <- rutils::diff_it(as.vector(x_ts))
#  re_turns[1] <- 0
#  colnames(re_turns) <- "diffs"
  rets_advanced <- rutils::lag_it(re_turns, -1)
#  rets_advanced[NROW(rets_advanced)] <- 0
#  colnames(rets_advanced) <- "rets_advanced"

# calculate volatility as the rolling quantile of returns
  quan_tile <- caTools::runquantile(x=abs(re_turns), k=win_dow,
                        probs=0.9, endrule="constant", align="center")
#  quan_tile <- xts(quan_tile, order.by=index(re_turns))
#  colnames(quan_tile) <- "volat"
# carry forward non-zero quan_tile values
  quan_tile[quan_tile==0] <- NA
  quan_tile[1] <- 1
  quan_tile <- rutils::na_locf(quan_tile)
#  quan_tile <- rutils::na_locf(quan_tile, fromLast=TRUE)

# value is suspect if abs re_turns greater than quan_tile, 
# and if abs sum of re_turns less than quan_tile
  sus_pect <- ((abs(re_turns) > vol_mult*quan_tile) &
      (abs(rets_advanced) > vol_mult*quan_tile) &
      (abs(re_turns+rets_advanced) < 2*vol_mult*quan_tile))
  sus_pect[1] <- FALSE
#  colnames(sus_pect) <- "suspect"
# cat("Parsing", deparse(substitute(ta_q)), "\n")
# cat("Parsing", strsplit(deparse(substitute(ta_q)), split="[.]")[[1]][4], "on date:", format(to_day), "\tfound", sum(sus_pect), "suspect prices\n")
  cat("date:", format(as.Date(index(first(x_ts)))), "\tfound", sum(sus_pect), "jump prices\n")
  sus_pect
}  # end which_jumps




#' Scrub a single day of \emph{TAQ} data in \emph{xts} format, without
#' aggregation.
#'
#' @export
#' @inheritParams which_extreme
#' @param ta_q \emph{TAQ} time series in \emph{xts} format.
#' @param tzone timezone to convert.
#' @return A \emph{TAQ} time series in \emph{xts} format.
#' @details The function \code{scrub_taq()} performs the same scrubbing
#'   operations as \code{scrub_agg}, except it doesn't aggregate, and returns
#'   the \emph{TAQ} data in \emph{xts} format.
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
  ta_q <- ta_q["T09:30:00/T16:00:00", ]
# return NULL if no data
  if (NROW(ta_q)==0)  return(NULL)
#  to_day <- as.Date(index(first(ta_q)))

# remove duplicate time stamps using duplicated()
  ta_q <- ta_q[!duplicated(index(ta_q)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- which_extreme(bid_offer, win_dow=win_dow, vol_mult=vol_mult)
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
  ta_q[which_jumps(mid_prices, win_dow=win_dow, vol_mult=vol_mult), ] <- NA
  rutils::na_locf(ta_q)
}  # end scrub_taq




#' Scrub a single day of \emph{TAQ} data, aggregate it, and convert to
#' \emph{OHLC} format.
#'
#' @export
#' @inheritParams scrub_taq
#' @param period aggregation period.
#' @return A \emph{OHLC} time series in \emph{xts} format.
#' @details The function \code{scrub_agg()} performs:
#' \itemize{
#'   \item index timezone conversion,
#'   \item data subset to trading hours,
#'   \item removal of duplicate time stamps,
#'   \item scrubbing of quotes with suspect bid-offer spreads,
#'   \item scrubbing of quotes with suspect price jumps,
#'   \item cbinding of mid prices with volume data,
#'   \item aggregation to OHLC using function \code{to.period()} from package \emph{xts},
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
  ta_q <- ta_q["T09:30:00/T16:00:00", ]
# return NULL if no data
  if (NROW(ta_q)==0)  return(NULL)
#  to_day <- as.Date(index(first(ta_q)))

# remove duplicate time stamps using duplicated()
  ta_q <- ta_q[!duplicated(index(ta_q)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- which_extreme(bid_offer, win_dow=win_dow, vol_mult=vol_mult)
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
  mid_prices[which_jumps(mid_prices, win_dow=win_dow, vol_mult=vol_mult)] <- NA
  mid_prices <- rutils::na_locf(mid_prices)
#  mid_prices <- rutils::na_locf(mid_prices, fromLast=TRUE)
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




#' Load, scrub, aggregate, and rbind multiple days of \emph{TAQ} data for a
#' single symbol, and save the \emph{OHLC} time series to a single
#' \sQuote{\code{*.RData}} file.
#'
#' @export
#' @param sym_bol \emph{character} string representing symbol or ticker.
#' @param data_dir \emph{character} string representing directory containing
#'   input \sQuote{\code{*.RData}} files.
#' @param output_dir \emph{character} string representing directory containing
#'   output \sQuote{\code{*.RData}} files.
#' @inheritParams scrub_agg
#' @return An \emph{OHLC} time series in \emph{xts} format.
#' @details The function \code{save_scrub_agg()} loads multiple days of
#'   \emph{TAQ} data, then scrubs, aggregates, and rbinds them into a
#'   \emph{OHLC} time series, and finally saves it to a single
#'   \sQuote{\code{*.RData}} file. The \emph{OHLC} time series is stored in a
#'   variable named \sQuote{\code{symbol}}, and then it's saved to a file named
#'   \sQuote{\code{symbol.RData}} in the \sQuote{\code{output_dir}} directory.
#'   The \emph{TAQ} data files are assumed to be stored in separate directories
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}}
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
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




#' Load and scrub multiple days of \emph{TAQ} data for a single symbol, and save
#' it to multiple \sQuote{\code{*.RData}} files.
#'
#' @export
#' @inheritParams save_scrub_agg
#' @return A \emph{TAQ} time series in \emph{xts} format.
#' @details The function \code{save_taq()} loads multiple days of \emph{TAQ}
#'   data, scrubs it, and saves the scrubbed TAQ data to individual
#'   \sQuote{\code{*.RData}} files. It uses the same file names for output as
#'   the input file names. The \emph{TAQ} data files are assumed to be stored in
#'   separate directories for each \sQuote{\code{symbol}}. Each
#'   \sQuote{\code{symbol}} has its own directory (named \sQuote{\code{symbol}})
#'   in the \sQuote{\code{data_dir}} directory.
#'   Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
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




#' Load, scrub, aggregate, and rbind multiple days of \emph{TAQ} data for a 
#' single symbol. Calculate returns and save them to a single
#' \sQuote{\code{*.RData}} file.
#'
#' @export
#' @inheritParams save_scrub_agg
#' @return A time series of returns and volume in \emph{xts} format.
#' @details The function \code{save_rets} loads multiple days of \emph{TAQ}
#'   data, then scrubs, aggregates, and rbinds them into a \emph{OHLC} time
#'   series.  It then calculates returns using function \code{run_returns()}, and
#'   stores them in a variable named \sQuote{\code{symbol.rets}}, and saves them
#'   to a file called \sQuote{\code{symbol.rets.RData}}.
#'   The \emph{TAQ} data files are assumed to be stored in separate directories
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}}
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
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




#' Load \emph{OHLC} time series data for a single symbol, calculate its returns,
#' and save them to a single \sQuote{\code{*.RData}} file, without aggregation.
#'
#' @export
#' @inheritParams save_scrub_agg
#' @return A time series of returns and volume in \emph{xts} format.
#' @details The function \code{save_rets_ohlc()} loads \emph{OHLC} time series
#'   data from a single file.  It then calculates returns using function
#'   \code{run_returns()}, and stores them in a variable named
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




#' Calculate a time series of point estimates of variance for an \emph{OHLC}
#' time series, using different range estimators for variance.
#'
#' Calculates the point variance estimates from individual bars of \emph{OHLC} 
#' prices (rows of data), using the squared differences of \emph{OHLC} prices at
#' each point in time, without averaging them over time.
#'
#' @export
#' @param oh_lc an \emph{OHLC} time series of prices in \emph{xts} format.
#' @param calc_method \emph{character} string representing the method for
#'   estimating variance.  The methods include:
#'   \itemize{
#'     \item "close" close to close,
#'     \item "garman_klass" Garman-Klass,
#'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
#'     \item "rogers_satchell" Rogers-Satchell,
#'     \item "yang_zhang" Yang-Zhang,
#'    }
#'    (default is \code{"yang_zhang"})
#' @param sca_le \emph{Boolean} argument: should the returns be divided by the 
#'   number of seconds in each period? (default is \code{TRUE})
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{run_variance()} calculates a time series of point
#'   variance estimates of percentage returns, from \emph{OHLC} prices, without
#'   averaging them over time. For example, the method \code{"close"} simply
#'   calculates the squares of the differences of the log \emph{Close} prices.
#'   
#'   The other methods calculate the squares of other possible differences of 
#'   the log \emph{OHLC} prices.  This way the point variance estimates only
#'   depend on the price differences within individual bars of data (and
#'   possibly from the neighboring bars.)
#'   All the methods are implemented assuming zero drift, since the calculations
#'   are performed only for a single bar of data, at a single point in time.
#'   
#'   The user can choose from several different variance estimation methods. The
#'   methods \code{"close"}, \code{"garman_klass_yz"}, and \code{"yang_zhang"}
#'   do account for close-to-open price jumps, while the methods
#'   \code{"garman_klass"} and \code{"rogers_satchell"} do not account for
#'   close-to-open price jumps. The default method is \code{"yang_zhang"}, which
#'   theoretically has the lowest standard error among unbiased estimators. 
#'
#'   The point variance estimates can be passed into function \code{roll_vwap()}
#'   to perform averaging, to calculate rolling variance estimates.  This is
#'   appropriate only for the methods \code{"garman_klass"} and
#'   \code{"rogers_satchell"}, since they don't require subtracting the rolling
#'   mean from the point variance estimates.
#'   
#'   The point variance estimates can also be considered to be technical
#'   indicators, and can be used as inputs into trading models.
#'   
#'   If \code{sca_le} is \code{TRUE} (the default), then the variance is divided
#'   by the squared differences of the time index (which scales the variance to 
#'   units of variance per second squared.) This is useful for example, when
#'   calculating intra-day variance from minutely bar data, because dividing
#'   returns by the number of seconds decreases the effect of overnight price 
#'   jumps.
#'   
#'   If \code{sca_le} is \code{TRUE} (the default), then the variance is 
#'   expressed in the scale of the time index of the \emph{OHLC} time series. 
#'   For example, if the time index is in seconds, then the variance is given in
#'   units of variance per second squared.  If the time index is in days, then
#'   the variance is equal to the variance per day squared.
#'   
#'   The time index of the \code{oh_lc} time series is assumed to be in 
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'   
#'   The function \code{run_variance()} performs similar calculations to the
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

run_variance <- function(oh_lc, calc_method="yang_zhang", sca_le=TRUE) {
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
         "yang_zhang"={co_eff <- 0.34/(1.34 + (win_dow + 1)/(win_dow - 1))
                            (oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]))^2 +
                            co_eff*(oh_lc[, 1]-oh_lc[, 4])^2 +
                            (1-co_eff)*((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                               (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]))}
  )  # end switch
  if (sca_le)
    vari_ance <- vari_ance/c(1, diff(xts::.index(oh_lc)))^2
  vari_ance[1, ] <- 0
  vari_ance <- rutils::na_locf(vari_ance)
  colnames(vari_ance) <- paste0(sym_bol, ".Variance")
  vari_ance
}  # end run_variance




#' Calculate time series of skew estimates from a \emph{OHLC} time series,
#' assuming zero drift.
#'
#' @export
#' @param oh_lc an \emph{OHLC} time series of prices in \emph{xts} format.
#' @param calc_method \emph{character} string representing method for estimating
#'   skew.
#' @return A time series of skew estimates.
#' @details The function \code{run_skew()} calculates a time series of skew
#'   estimates from \emph{OHLC} prices, one for each bar of \emph{OHLC} data.
#'   The skew estimates are expressed in the time scale of the index of the
#'   \emph{OHLC} time series.  
#'   For example, if the time index is in seconds, then the skew is given in 
#'   units of skew per second.  If the time index is in days, then the skew is 
#'   equal to the skew per day.
#'   
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
  sk_ew <- sk_ew/c(1, diff(xts::.index(oh_lc)))^3
  sk_ew[1, ] <- 0
  sk_ew <- rutils::na_locf(sk_ew)
  colnames(sk_ew) <- paste0(sym_bol, ".Skew")
  sk_ew
}  # end run_skew




#' Calculate time series of Sharpe-like statistics for each bar of a \emph{OHLC}
#' time series.
#'
#' @export
#' @param oh_lc an \emph{OHLC} time series of prices in \emph{xts} format.
#' @param calc_method \emph{character} string representing method for estimating
#'   the Sharpe-like exponent.
#' @return An \emph{xts} time series with the same number of rows as the
#'   argument \code{oh_lc}.
#' @details The function \code{run_sharpe()} calculates Sharpe-like statistics
#'   for each bar of a \emph{OHLC} time series.
#'   The Sharpe-like statistic is defined as the ratio of the difference between
#'   \emph{Close} minus \emph{Open} prices divided by the difference between
#'   \emph{High} minus \emph{Low} prices.
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
#' a \emph{OHLC} time series.
#'
#' @export
#' @param oh_lc \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#' @param mo_ment \emph{character} string representing function for
#'   estimating the moment.
#' @param weight_ed \emph{Boolean} argument: should estimate be weighted by
#'   the trading volume? (default is \code{TRUE})
#' @param ... additional parameters to the mo_ment function.
#' @return A single \emph{numeric} value equal to the volume weighted average of
#'   an estimator over the time series.
#' @details The function \code{agg_regate()} calculates a single number
#'   representing the volume weighted average of an estimator over the
#'   \emph{OHLC} time series of prices.  By default the sum is trade volume
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
    agg_regations <- sum(agg_regations)/NROW(agg_regations)
  agg_regations
}  # end agg_regate




#' Calculate the volume-weighted average price of an \emph{OHLC} time series
#' over a rolling window (lookback period).
#'
#' Performs the same operation as function \code{VWAP()} from package
#' \href{https://cran.r-project.org/web/packages/TTR/index.html}{VWAP},
#' but using vectorized functions, so it's a little faster.
#'
#' @export
#' @param oh_lc an \emph{OHLC} time series of prices in \emph{xts} format.
#' @param x_ts single-column \emph{xts} time series.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for calculating the average price.
#' @return An \emph{xts} time series with a single column and the same number of
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




#' Calculate a vector of statistics over an \emph{OHLC} time series, and
#' calculate a rolling mean over the statistics.
#'
#' @export
#' @param oh_lc \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#' @param mo_ment \emph{character} string representing a function for
#'   estimating statistics of a single bar of \emph{OHLC} data, such as
#'   volatility, skew, and higher moments.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for calculating the rolling mean.
#' @param weight_ed \emph{Boolean} argument: should statistic be weighted by
#'   trade volume? (default \code{TRUE})
#' @param ... additional parameters to the mo_ment function.
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{roll_moment()} calculates a vector of statistics
#'   over an \emph{OHLC} time series, such as volatility, skew, and higher
#'   moments.  The statistics could also be any other aggregation of a single
#'   bar of \emph{OHLC} data, for example the \emph{High} price minus the
#'   \emph{Low} price squared.  The length of the vector of statistics is equal
#'   to the number of rows of the argument \code{oh_lc}. Then it calculates a
#'   trade volume weighted rolling mean over the vector of statistics over and
#'   calculate statistics.
#' @examples
#' # calculate time series of rolling variance and skew estimates
#' var_rolling <- roll_moment(oh_lc=SPY, win_dow=21)
#' skew_rolling <- roll_moment(oh_lc=SPY, mo_ment="run_skew", win_dow=21)
#' skew_rolling <- skew_rolling/(var_rolling)^(1.5)
#' skew_rolling[1, ] <- 0
#' skew_rolling <- rutils::na_locf(skew_rolling)

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




#' Calculate a time series of variance estimates over a rolling lookback window 
#' for an \emph{OHLC} time series of prices, using different range estimators
#' for variance.
#'
#' @export
#' @param oh_lc an \emph{OHLC} time series of prices in \emph{xts} format.
#' @param calc_method \emph{character} string representing method for estimating
#'   variance.  The methods include:
#'   \itemize{
#'     \item "close" close to close,
#'     \item "garman_klass" Garman-Klass,
#'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
#'     \item "rogers_satchell" Rogers-Satchell,
#'     \item "yang_zhang" Yang-Zhang,
#'    }
#'    (default is \code{"yang_zhang"})
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for calculating the variance. 
#' @param sca_le \emph{Boolean} argument: should the returns be divided by the 
#'   number of seconds in each period? (default is \code{TRUE})
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{roll_variance()} calculates a time series of 
#'   variance estimates of percentage returns, from \emph{OHLC} prices, using
#'   several different variance estimation methods based on the range of
#'   \emph{OHLC} prices.
#'   
#'   If \code{sca_le} is \code{TRUE} (the default), then the variance is divided
#'   by the squared differences of the time index (which scales the variance to 
#'   units of variance per second squared.) This is useful for example, when
#'   calculating intra-day variance from minutely bar data, because dividing
#'   returns by the number of seconds decreases the effect of overnight price 
#'   jumps.
#'   
#'   If \code{sca_le} is \code{TRUE} (the default), then the variance is 
#'   expressed in the scale of the time index of the \emph{OHLC} time series. 
#'   For example, if the time index is in seconds, then the variance is given in
#'   units of variance per second squared.  If the time index is in days, then
#'   the variance is equal to the variance per day squared.
#'   
#'   The time index of the \code{oh_lc} time series is assumed to be in 
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'   
#'   The methods \code{"close"}, \code{"garman_klass_yz"}, and
#'   \code{"yang_zhang"} do account for close-to-open price jumps, while the
#'   methods \code{"garman_klass"} and \code{"rogers_satchell"} do not account
#'   for close-to-open price jumps. 
#'   
#'   The default method is \code{"yang_zhang"}, which theoretically has the
#'   lowest standard error among unbiased estimators.
#'
#'   The function \code{roll_variance()} performs the same calculations as the 
#'   function \code{volatility()} from package 
#'   \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR}, but 
#'   it's a little faster because it uses function RcppRoll::roll_sd(), and it 
#'   performs less data validation.
#' @examples
#' # create minutely OHLC time series of random prices
#' oh_lc <- HighFreq::random_ohlc()
#' # calculate variance estimates for oh_lc over a 21 period window
#' var_rolling <- HighFreq::roll_variance(oh_lc, win_dow=21)
#' # calculate variance estimates for SPY
#' var_rolling <- HighFreq::roll_variance(SPY, calc_method="yang_zhang")
#' # calculate SPY variance without accounting for overnight jumps
#' var_rolling <- HighFreq::roll_variance(SPY, calc_method="rogers_satchell")

roll_variance <- function(oh_lc, win_dow=11, calc_method="yang_zhang", sca_le=TRUE) {
  sym_bol <- rutils::na_me(oh_lc)
  oh_lc <- log(oh_lc[, 1:4])
  vari_ance <- switch(calc_method,
                      "close"={xts(c(rep(0, win_dow-1), 
                            RcppRoll::roll_var(rutils::diff_xts(oh_lc[, 4]), n=win_dow, align="left")),
                              order.by=index(oh_lc))},
                      "garman_klass"={rutils::roll_sum(
                            0.5*(oh_lc[, 2]-oh_lc[, 3])^2 -
                            (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^2, win_dow=win_dow) / win_dow},
                      "rogers_satchell"={rutils::roll_sum((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                            (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]), win_dow=win_dow) / win_dow},
                      "garman_klass_yz"={rutils::roll_sum(
                            (oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]))^2 + 0.5*(oh_lc[, 2]-oh_lc[, 3])^2 -
                            (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^2, win_dow=win_dow) / win_dow},
                      "yang_zhang"={co_eff <- 0.34/(1.34 + (win_dow + 1)/(win_dow - 1))
                            c(rep(0, win_dow-1), 
                            RcppRoll::roll_var(oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]), n=win_dow, align="left") +
                            co_eff*RcppRoll::roll_var(oh_lc[, 1]-oh_lc[, 4], n=win_dow, align="left")) +
                            (1-co_eff)*rutils::roll_sum((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                              (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]), win_dow=win_dow) / win_dow}
  )  # end switch
  if (sca_le)
    vari_ance <- vari_ance/c(1, diff(xts::.index(oh_lc)))^2
  vari_ance[1, ] <- 0
  vari_ance <- rutils::na_locf(vari_ance)
  colnames(vari_ance) <- paste0(sym_bol, ".Variance")
  vari_ance
}  # end roll_variance




#' Calculate the variance of an \emph{OHLC} time series, using different range
#' estimators for variance.
#'
#' @export
#' @param oh_lc an \emph{OHLC} time series of prices in \emph{xts} format.
#' @param calc_method \emph{character} string representing method for estimating
#'   variance.  The methods include:
#'   \itemize{
#'     \item "close" close to close,
#'     \item "garman_klass" Garman-Klass,
#'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
#'     \item "rogers_satchell" Rogers-Satchell,
#'     \item "yang_zhang" Yang-Zhang,
#'    }
#'    (default is \code{"yang_zhang"})
#' @param sca_le \emph{Boolean} argument: should the returns be divided by the 
#'   number of seconds in each period? (default is \code{TRUE})
#' @return A \emph{numeric} value equal to the variance.
#' @details The function \code{calc_variance()} calculates the variance estimate
#'   from \emph{OHLC} prices, using several different variance estimation
#'   methods based on the range of \emph{OHLC} prices.
#'   
#'   The methods \code{"close"}, \code{"garman_klass_yz"}, and
#'   \code{"yang_zhang"} do account for close-to-open price jumps, while the
#'   methods \code{"garman_klass"} and \code{"rogers_satchell"} do not account
#'   for close-to-open price jumps. 
#'   
#'   The default method is \code{"yang_zhang"}, which theoretically has the
#'   lowest standard error among unbiased estimators.
#'
#'   If \code{sca_le} is \code{TRUE} (the default), then the variance is divided
#'   by the squared differences of the time index (which scales the variance to 
#'   units of variance per second squared.) This is useful for example, when 
#'   calculating variance from minutely bar data, because dividing returns by
#'   the number of seconds decreases the effect of overnight price jumps.
#'   
#'   If \code{sca_le} is \code{TRUE} (the default), then the variance is 
#'   expressed in the scale of the time index of the \emph{OHLC} time series. 
#'   For example, if the time index is in seconds, then the variance is given in
#'   units of variance per second squared.  If the time index is in days, then
#'   the variance is equal to the variance per day squared.
#'   
#'   The function \code{calc_variance()} performs the same calculations as the 
#'   function \code{run_variance()} and then calculates the average of the spot
#'   variance estimates.
#' @examples
#' # create minutely OHLC time series of random prices
#' oh_lc <- HighFreq::random_ohlc()
#' # calculate variance of oh_lc
#' vari_ance <- HighFreq::calc_variance(oh_lc)
#' # calculate variance of SPY
#' vari_ance <- HighFreq::calc_variance(SPY, calc_method="yang_zhang")
#' # calculate variance of SPY without accounting for overnight jumps
#' vari_ance <- HighFreq::calc_variance(SPY, calc_method="rogers_satchell")

calc_variance <- function(oh_lc, calc_method="yang_zhang", sca_le=TRUE) {
  if (sca_le)
    in_dex <- c(1, diff(xts::.index(oh_lc)))
  else
    in_dex <- rep(1, NROW(oh_lc))
  oh_lc <- log(oh_lc[, 1:4])
  switch(calc_method,
         "close"={var(rutils::diff_it(as.vector(oh_lc[, 4]))/in_dex)},
         "garman_klass"={sum((0.5*(oh_lc[, 2]-oh_lc[, 3])^2 -
                          (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^2)/in_dex^2) / NROW(oh_lc)},
         "rogers_satchell"={sum(((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                            (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]))/in_dex^2) / NROW(oh_lc)},
         "garman_klass_yz"={sum(((oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]))^2 + 0.5*(oh_lc[, 2]-oh_lc[, 3])^2 -
                          (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^2)/in_dex^2) / NROW(oh_lc)},
         "yang_zhang"={co_eff <- 0.34/(1.34 + (NROW(oh_lc) + 1)/(NROW(oh_lc) - 1))
                        drop(var((oh_lc[, 1]-rutils::lag_xts(oh_lc[, 4]))/in_dex) +
                          co_eff*var((oh_lc[, 1]-oh_lc[, 4])/in_dex) +
                          (1-co_eff)*sum(((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                                          (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]))/in_dex^2) / NROW(oh_lc))}
  )  # end switch
}  # end calc_variance




#' Calculate a time series of Sharpe ratios over a rolling lookback window for
#' an \emph{OHLC} time series.
#'
#' @export
#' @param oh_lc an \emph{OHLC} time series of prices in \emph{xts} format.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for aggregating the \emph{OHLC} prices.
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{roll_sharpe()} calculates the rolling Sharpe 
#'   ratio defined as the ratio of percentage returns over the lookback window,
#'   divided by the average volatility of percentage returns.
#' @examples
#' # calculate rolling Sharpe ratio over SPY
#' sharpe_rolling <- roll_sharpe(oh_lc=SPY, win_dow=11)

roll_sharpe <- function(oh_lc, win_dow=11) {
  re_turns <- run_returns(oh_lc, lag=win_dow, sca_le=FALSE)
  var_rolling <- sqrt(roll_variance(oh_lc, win_dow=win_dow, sca_le=FALSE))
  sharpe_rolling <- ifelse(var_rolling==0,
                           1.0,
                           re_turns/var_rolling)
  colnames(sharpe_rolling) <- paste0(rutils::na_me(oh_lc), ".Sharpe")
  rutils::na_locf(sharpe_rolling)
}  # end roll_sharpe




#' Calculate a time series of \emph{Hurst} exponents over a rolling lookback
#' window.
#'
#' @export
#' @param oh_lc an \emph{OHLC} time series of prices in \emph{xts} format.
#' @param win_dow the size of the lookback window, equal to the number of bars
#'   of data used for aggregating the \emph{OHLC} prices.
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#' @details The function \code{roll_hurst()} calculates a time series of
#'   \emph{Hurst} exponents from \emph{OHLC} prices, over a rolling lookback
#'   window.
#'   
#'   The \emph{Hurst} exponent is defined as the logarithm of the ratio of the
#'   price range, divided by the standard deviation of returns, and divided by
#'   the logarithm of the window length.
#'   
#'   The function \code{roll_hurst()} doesn't use the same definition as the 
#'   rescaled range definition of the \emph{Hurst} exponent.
#'   First, because the price range is calculated using \emph{High} and 
#'   \emph{Low} prices, which produces bigger range values, and higher
#'   \emph{Hurst} exponent estimates.
#'   Second, because the \emph{Hurst} exponent is estimated using a single
#'   aggregation window, instead of multiple windows in the rescaled range
#'   definition.
#'   
#'   The rationale for using a different definition of the \emph{Hurst} exponent
#'   is that it's designed to be a technical indicator for use as input into
#'   trading models, rather than an estimator for statistical analysis.
#'   
#' @examples
#' # calculate rolling Hurst for SPY in March 2009
#' hurst_rolling <- roll_hurst(oh_lc=SPY["2009-03"], win_dow=11)
#' chart_Series(hurst_rolling["2009-03-10/2009-03-12"], name="SPY hurst_rolling")

roll_hurst <- function(oh_lc, win_dow=11) {
  ran_ge <- c(rep(0, win_dow-1), (RcppRoll::roll_max(x=log(oh_lc[, 2]), n=win_dow) + 
               RcppRoll::roll_max(x=-log(oh_lc[, 3]), n=win_dow)))
  var_rolling <- sqrt(roll_variance(oh_lc, win_dow=win_dow, sca_le=FALSE))
  hurst_rolling <- ifelse((var_rolling==0) | (ran_ge==0),
                          0.5,
                          log(ran_ge/var_rolling)/log(win_dow))
  colnames(hurst_rolling) <- paste0(rutils::na_me(oh_lc), ".Hurst")
  rutils::na_locf(hurst_rolling)
}  # end roll_hurst




#' Apply an aggregation function over a rolling lookback window and the end
#' points of an \emph{OHLC} time series.
#'
#' @export
#' @param x_ts \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#' @param agg_fun \emph{character} string representing an aggregation function
#'   to be applied over a rolling lookback window.
#' @param win_dow the size of the lookback window, equal to the number of bars 
#'   of data used for applying the aggregation function (including the current
#'   bar).
#' @param by_columns \emph{Boolean} argument: should the function
#'   \code{agg_fun()} be applied column-wise (individually), or should it be
#'   applied to all the columns combined? (default is \code{FALSE})
#' @param end_points an integer vector of end points.
#' @param ... additional parameters to the agg_fun function.
#' @return An \emph{xts} time series with the same number of rows as the
#'   argument \code{x_ts}.
#' @details The function \code{roll_apply()} applies an aggregation function
#'   over a rolling lookback window and the end points of an \emph{OHLC} time
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
#'   \code{roll_apply()} doesn't produce any leading \emph{NA} values.
#'
#'   The function \code{roll_apply()} can be called in two different ways,
#'   depending on the argument \code{end_points}.
#'   If the argument \code{end_points} isn't explicitly passed to
#'   \code{roll_apply()}, then the default value is used, and
#'   \code{roll_apply()} performs aggregations over overlapping windows at each
#'   point in time.
#'   If the argument \code{end_points} is explicitly passed to 
#'   \code{roll_apply()}, then \code{roll_apply()} performs aggregations over 
#'   windows spanned by the end_points.  If win_dow=2 then the aggregations are
#'   performed over non-overlapping windows, otherwise they are performed over
#'   overlapping windows.
#'
#'   The aggregation function \code{agg_fun()} can return either a single value
#'   or a vector of values. If the aggregation function \code{agg_fun()} returns
#'   a single value, then \code{roll_apply()} returns an \emph{xts} time series 
#'   with a single column. If the aggregation function \code{agg_fun()} returns
#'   a vector of values, then \code{roll_apply()} returns an \emph{xts} time 
#'   series with multiple columns equal to the length of the vector returned by 
#'   the aggregation function \code{agg_fun()}.
#'
#' @examples
#' # extract a single day of SPY data
#' oh_lc <- SPY["2012-02-13"]
#' win_dow <- 11
#' # calculate the rolling sums of oh_lc columns over a rolling window
#' agg_regations <- roll_apply(oh_lc, agg_fun=sum, win_dow=win_dow, by_columns=TRUE)
#' # apply a vector-valued aggregation function over a rolling window
#' agg_function <- function(oh_lc)  c(max(oh_lc[, 2]), min(oh_lc[, 3]))
#' agg_regations <- roll_apply(oh_lc, agg_fun=agg_function, win_dow=win_dow)
#' # define end points at 11-minute intervals (SPY is minutely bars)
#' end_points <- rutils::end_points(oh_lc, inter_val=win_dow)
#' # calculate the sums of oh_lc columns over end_points using non-overlapping windows
#' agg_regations <- roll_apply(oh_lc, agg_fun=sum, win_dow=2, 
#'                             end_points=end_points, by_columns=TRUE)
#' # apply a vector-valued aggregation function over the end_points of oh_lc
#' # using overlapping windows
#' agg_regations <- roll_apply(oh_lc, agg_fun=agg_function, 
#'                             win_dow=5, end_points=end_points)

roll_apply <- function(x_ts, agg_fun="run_variance", win_dow=11,
                       end_points=(0:NROW(x_ts)), by_columns=FALSE, ...) {
  # match "agg_fun" with some aggregation function
  agg_fun <- match.fun(agg_fun)
  len_gth <- NROW(end_points)
  # define start_points as lag of end_points
  start_points <-  end_points[c(rep_len(1, win_dow-1), 1:(len_gth-win_dow+1))] +
    (NROW(x_ts) > (len_gth+1))
  # perform aggregations over length of end_points
  agg_regations <- if(by_columns)
    sapply(x_ts, function(col_umn)
      sapply(2:len_gth, function(in_dex)
        agg_fun(.subset_xts(col_umn,
                            start_points[in_dex]:end_points[in_dex]), ...)
      ))  # end sapply
  else {  # not by_columns
    agg_regations <- sapply(2:len_gth, function(in_dex)
      agg_fun(.subset_xts(x_ts,
                          start_points[in_dex]:end_points[in_dex]), ...)
    )  # end sapply
    # coerce agg_regations into matrix and transpose it
    if (is.vector(agg_regations))
      agg_regations <- t(agg_regations)
    agg_regations <- t(agg_regations)
  }  # end if
  # coerce agg_regations into xts series
  xts(agg_regations, order.by=index(x_ts[end_points]))
}  # end roll_apply




#' Perform seasonality aggregations over a single-column \emph{xts} time series.
#'
#' @export
#' @param x_ts single-column \emph{xts} time series.
#' @param in_dex vector of \emph{character} strings representing points in time,
#'   of the same length as the argument \code{x_ts}.
#' @return An \emph{xts} time series with mean aggregations over the seasonality
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

season_ality <- function(x_ts, in_dex=format(zoo::index(x_ts), "%H:%M")) {
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


