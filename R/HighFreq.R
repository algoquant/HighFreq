################################################
###
###  scrubbing and aggregating HFREQ data
###
################################################


###########
# functions


### recursively "rbind" a list of xts time series - same as do.call.rbind
# call lapply in a loop to divide list_var by half, binding neighboring elements
#' Recursively "rbind" a list of objects.
#' 
#' Performs the same operation as \code{do.call(rbind, list_var)}, but using recursion, which is much faster and uses less memory.
#' This is the same function as '\code{do.call.rbind}' from package 'qmao'.
#' 
#' @param list_var A list of \code{vectors}, \code{matrices}, \code{data frames}, or \code{time series}.
#' @return  a single \code{vector}, \code{matrix}, \code{data frame}, or \code{time series}.
#' @details Performing a \code{do.call(rbind, x_ts)} on a list of 'xts' time series is very slow, and often causes an "out of memory" error.
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


### identify extreme values in a univariate xts time series
#' Identify extreme values in a univariate xts time series.
#' 
#' Calculates running volatility defined as a quantile over a sliding window.
#' Identifies extreme values as those that exceed the scaled running volatility.
#' 
#' @param time_series univariate \code{xts} time series.
#' @param vol_window number of data points for estimating running volatility.
#' @param thresh_old level of prices above which values are considered extreme.
#' @return  \code{logical vector}.
#' @details Extreme values are those that exceed a fixed multiple of the running volatility.
#' @details Extreme values are the very tips of the tails when the distribution of values becomes very fat-tailed.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' # scrub extreme values
#' x_ts <- x_ts[!extreme_values(x_ts, thresh_old=1)]
extreme_values <- function(time_series, vol_window=51, thresh_old=2) {

# calculate volatility as running quantile
  vo_lat <- runquantile(x=abs(as.vector(time_series)), k=vol_window, probs=0.9, endrule="constant", align="center")
  vo_lat <- xts(vo_lat, order.by=index(time_series))
  colnames(vo_lat) <- "volat"
# carry forward non-zero volatility values
  vo_lat[vo_lat==0] <- NA
  vo_lat[1] <- 1
  vo_lat <- na.locf(vo_lat)
#  vo_lat <- na.locf(vo_lat, fromLast=TRUE)

# extreme value if time_series greater than scaled volatility
  ex_treme <- (abs(time_series) > 2*thresh_old*vo_lat)
  ex_treme[1] <- FALSE
  colnames(ex_treme) <- "suspect"

  cat("date:", format(as.Date(index(first(time_series)))), "\tfound", sum(ex_treme), "extreme values\n")
  ex_treme
}  # end extreme_values



### identify jump values in a univariate xts price time series
#' Identify jump values in a univariate xts time series.
#' 
#' Calculates the running volatility of simple returns, defined as a quantile over a sliding window.
#' Identifies jump values as those returns that exceed the scaled running volatility, but their sum with the lagged returns doesn't exceed it.
#' 
#' @param time_series univariate \code{xts} time series.
#' @param vol_window number of data points for estimating running volatility.
#' @param thresh_old level of change above which values are considered extreme.
#' @return  \code{logical vector}.
#' @details Jump values are those where neighboring values both exceed the scaled running volatility (defined as a quantile), but their sum doesn't exceed it.
#' @details Jump values are simply isolated prices out of proportion to their neighbors.
#' @examples
#' # create xts time series
#' x_ts <- xts(x=rnorm(1000), order.by=(Sys.time()-3600*(1:1000)))
#' # scrub jump values
#' x_ts <- x_ts[!jump_values(x_ts, thresh_old=1.0)]
jump_values <- function(time_series, vol_window=51, thresh_old=2) {

# calculate simple returns
  diff_series <- diff(time_series)
  diff_series[1, ] <- 0
  colnames(diff_series) <- "diffs"
  diff_series_fut <- lag(diff_series, -1)
  diff_series_fut[nrow(diff_series_fut)] <- 0
  colnames(diff_series_fut) <- "diff_series_fut"

# calculate vo_lat as running quantile
  vo_lat <- runquantile(x=abs(as.vector(diff_series)), k=vol_window, probs=0.9, endrule="constant", align="center")
  vo_lat <- xts(vo_lat, order.by=index(diff_series))
  colnames(vo_lat) <- "volat"
# carry forward non-zero vo_lat values
  vo_lat[vo_lat==0] <- NA
  vo_lat[1] <- 1
  vo_lat <- na.locf(vo_lat)
#  vo_lat <- na.locf(vo_lat, fromLast=TRUE)

# suspect value if abs diffs greater than vo_lat, and if abs sum of diffs less than vo_lat
  sus_pect <- (
    (abs(diff_series) > thresh_old*vo_lat) & 
      (abs(diff_series_fut) > thresh_old*vo_lat) & 
      (abs(diff_series+diff_series_fut) < 2*thresh_old*vo_lat)
  )
  sus_pect[1] <- FALSE
  colnames(sus_pect) <- "suspect"
# cat("Parsing", deparse(substitute(taq_data)), "\n")
# cat("Parsing", strsplit(deparse(substitute(taq_data)), split="[.]")[[1]][4], "on date:", format(to_day), "\tfound", sum(sus_pect), "suspect values\n")
  cat("date:", format(as.Date(index(first(time_series)))), "\tfound", sum(sus_pect), "jump values\n")
  sus_pect
}  # end jump_values



### scrub and aggregate a single day of TAQ data in xts format
# return mid price and volume
#' Scrub a single day of TAQ data and aggregate to OHLC format.
#' 
#' Calculate mid prices from bid and offer TAQ data, scrub them and cbind with volume data, and aggregate to 'period' OHLC data.
#' 
#' @param taq_data \code{TAQ xts} time series.
#' @param vol_window number of data points for estimating running volatility.
#' @param thresh_old level of change above which values are extreme.
#' @param period aggregation period.
#' @return  \code{OHLC xts} time series.
#' @details Valid 'period' character strings include: "minutes", "3 min", "5 min", "10 min", "15 min", "30 min", and "hours".
#' @details The time index is rounded up to the next integer multiple of 'period'.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"), to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random prices
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
scrub_agg <- function(taq_data, vol_window=51, thresh_old=2, period="minutes") {

# convert timezone of index to New_York
  index(taq_data) <- with_tz(index(taq_data), "America/New_York")
# subset data to NYSE trading hours
  taq_data <- taq_data['T09:30:00/T16:00:00', ]
# return NULL if no data
  if (nrow(taq_data)==0)  return(NULL)
  to_day <- as.Date(index(first(taq_data)))

# remove duplicate time stamps using 'duplicated'
  taq_data <- taq_data[!duplicated(index(taq_data)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- taq_data[, 'Ask.Price'] - taq_data[, 'Bid.Price']
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- extreme_values(bid_offer, vol_window=vol_window, thresh_old=thresh_old)
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
  mid_prices[jump_values(mid_prices)] <- NA
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



### scrub and return a single day of TAQ data
#' Scrub a single day of TAQ data in xts format.
#' 
#' @param taq_data \code{TAQ xts} time series.
#' @param vol_window number of data points for estimating running volatility.
#' @param thresh_old level of change above which values are extreme.
#' @return  \code{TAQ xts} time series.
#' @details Jump values are those where neighboring values both exceed the scaled running volatility (defined as a quantile), but their sum doesn't exceed it.
#' @details Jump values are simply isolated prices out of proportion to their neighbors.
#' @examples
#' # create time index of one second intervals for a single day
#' in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"), to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
#' # create xts of random prices
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
#' taq_data <- scrub_TAQ(taq_data, vol_window=11, thresh_old=1)
scrub_TAQ <- function(taq_data, vol_window=51, thresh_old=2) {

# convert time index to New_York
  index(taq_data) <- with_tz(index(taq_data), "America/New_York")
# subset data to NYSE trading hours
  taq_data <- taq_data['T09:30:00/T16:00:00', ]
# return NULL if no data
  if (nrow(taq_data)==0)  return(NULL)
  
# remove duplicate time stamps using duplicated
  taq_data <- taq_data[!duplicated(index(taq_data)), ]
  
# scrub quotes with suspect bid-offer spreads
  bid_offer <- taq_data[, 'Ask.Price'] - taq_data[, 'Bid.Price']
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- extreme_values(bid_offer, vol_window=vol_window, thresh_old=thresh_old)
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
# replace suspect values with NA
  mid_prices[jump_values(mid_prices)] <- NA
  mid_prices <- na.locf(mid_prices)
#  mid_prices <- na.locf(mid_prices, fromLast=TRUE)
  mid_prices <- cbind(mid_prices, taq_data[index(mid_prices), "Volume"])
  mid_prices[is.na(mid_prices[, "Volume"]), "Volume"] <- 0

# aggregate to OHLC minutes data and cumulative volume
  mid_prices <- to.period(x=mid_prices, period="minutes")
# round up times to next minute
  index(mid_prices) <- align.time(x=index(mid_prices), 60)
  mid_prices
}  # end scrub_TAQ



###########
# load and scrub multiple days of data for a single symbol
#' Load, scrub, rbind, and aggregate multiple days of xts time series data for a single symbol.
#' 
#' @param sym_bol \code{character string} representing symbol or ticker.
#' @param data_dir \code{character string} representing directory containing *.RData files.
#' @return  \code{OHLC xts} time series.
#' @details The high frequency data is assumed to be stored in separate directories for each symbol, containing *.RData files, each file with one day of data.  
#' @details Each symbol has its own directory, and the symbol directories reside in 'data_dir'.
#' @details The xts data is loaded, scrubbed, rbinded, and aggregated to OHLC format. 
#' @details The OHLC data is assigned to a variable named 'sym_bol', and then saved to a *.RData file in the 'cwd'.
#' @examples
#' save_OHLC("SPY")
save_OHLC <- function(sym_bol, data_dir="E:/mktdata/sec/", period="minutes") {

# create path to directory containing *.RData files
  file_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
  file_list <- list.files(file_dir)
# create paths to *.RData files
  file_names <- file.path(file_dir, file_list)

# load data into list
  data <- sapply(file_names, function(file_name) {
    cat("loading", sym_bol, "frome file: ", file_name, "\n")
    data_name <- load(file_name)
    get(data_name)
  })

# scrub and aggregate the data
  data <- sapply(data, scrub_agg, period=period)

# recursively "rbind" the list into a single xts
  data <- do_call_rbind(data)
# assign proper column names
  colnames(data) <- sapply(strsplit(colnames(data), split="[.]"), 
                           function(strng) paste(sym_bol, strng[-1], sep="."))

# copy the xts data to a variable with the name 'sym_bol'
  assign(sym_bol, data)

# save the xts data to a file in the 'cwd'
  save(list=eval(sym_bol), file=paste0(sym_bol, ".RData"))

}  # end save_OHLC

