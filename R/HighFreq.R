################################################
###
###  scrubbing and aggregating HFREQ data
###
################################################


###########
# functions


### recursively "rbind" a list xts time series - same as do.call.rbind
do_call_rbind <- function(list_var) {
  # call lapply in a loop to divide list_var by half, binding neighboring elements
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


### identify suspect bid_offer values in univariate xts time series
suspect_bid_offer <- function(bid_offer, agg_vol_window=51, suspect_threshold=2) {

  # calculate vo_lat as running quantile
  vo_lat <- runquantile(x=abs(as.vector(bid_offer)), k=agg_vol_window, probs=0.9, endrule="constant", align="center")
  vo_lat <- xts(vo_lat, order.by=index(bid_offer))
  colnames(vo_lat) <- "volat"
  # carry forward non-zero vo_lat values
  vo_lat[vo_lat==0] <- NA
  vo_lat[1] <- 1
  vo_lat <- na.locf(vo_lat)
  #  vo_lat <- na.locf(vo_lat, fromLast=TRUE)
  
  # find suspect values
  # suspect if bid_offer greater than vo_lat
  sus_pect <- (abs(bid_offer) > 2*suspect_threshold*vo_lat)
  sus_pect[1] <- FALSE
  
  cat("date:", format(as.Date(index(first(bid_offer)))), "\tscrubbed", sum(sus_pect), "suspect bid-offer values\n")
  sus_pect
}  # end suspect_bid_offer



### identify suspect jump values in univariate xts price time series
suspect_jump <- function(price_data, agg_vol_window=51, suspect_threshold=2) {

  # calculate simple returns
  diff_prices <- diff(price_data)
  diff_prices[1, ] <- 0
  colnames(diff_prices) <- "diffs"
  diff_prices_fut <- lag(diff_prices, -1)
  diff_prices_fut[nrow(diff_prices_fut)] <- 0
  colnames(diff_prices_fut) <- "diff_prices_fut"
  
  # calculate vo_lat as running quantile
  vo_lat <- runquantile(x=abs(as.vector(diff_prices)), k=agg_vol_window, probs=0.9, endrule="constant", align="center")
  vo_lat <- xts(vo_lat, order.by=index(diff_prices))
  colnames(vo_lat) <- "volat"
  # carry forward non-zero vo_lat values
  vo_lat[vo_lat==0] <- NA
  vo_lat[1] <- 1
  vo_lat <- na.locf(vo_lat)
  #  vo_lat <- na.locf(vo_lat, fromLast=TRUE)
  
  # find suspect values
  # suspect if abs diffs greater than vo_lat, and if abs sum of diffs less than vo_lat
  sus_pect <- (
    (abs(diff_prices) > suspect_threshold*vo_lat) & 
      (abs(diff_prices_fut) > suspect_threshold*vo_lat) & 
      (abs(diff_prices+diff_prices_fut) < 2*suspect_threshold*vo_lat)
  )
  sus_pect[1] <- FALSE
  colnames(sus_pect) <- "suspect"
  # cat("Parsing", deparse(substitute(taq_data)), "\n")
  # cat("Parsing", strsplit(deparse(substitute(taq_data)), split="[.]")[[1]][4], "on date:", format(to_day), "\tscrubbed", sum(sus_pect), "suspect values\n")
  cat("date:", format(as.Date(index(first(price_data)))), "\tscrubbed", sum(sus_pect), "suspect jump values\n")
  sus_pect
}  # end suspect_jump



### scrub and aggregate a single day of TAQ data in xts format
# return mid price and volume
scrub_agg <- function(taq_data, agg_vol_window=51, suspect_threshold=2) {

  # convert time index to New_York
  index(taq_data) <- with_tz(index(taq_data), "America/New_York")
  # subset data to NYSE trading hours
  taq_data <- taq_data['T09:30:00/T16:00:00', ]
  # return NULL if no data
  if (nrow(taq_data)==0)  return(NULL)
  to_day <- as.Date(index(first(taq_data)))
  
  # remove duplicate time stamps using duplicated
  taq_data <- taq_data[!duplicated(index(taq_data)), ]
  
  # scrub quotes with suspect bid-offer spreads
  bid_offer <- taq_data[, 'Ask.Price'] - taq_data[, 'Bid.Price']
  #  bid_offer <- na.omit(bid_offer)
  sus_pect <- suspect_bid_offer(bid_offer)
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
  mid_prices[suspect_jump(mid_prices)] <- NA
  mid_prices <- na.locf(mid_prices)
  #  mid_prices <- na.locf(mid_prices, fromLast=TRUE)
  mid_prices <- cbind(mid_prices, taq_data[index(mid_prices), "Volume"])
  mid_prices[is.na(mid_prices[, "Volume"]), "Volume"] <- 0
  
  # aggregate to OHLC minutes data and cumulative volume
  mid_prices <- to.period(x=mid_prices, period="minutes")
  # round up times to next minute
  index(mid_prices) <- align.time(x=index(mid_prices), 60)
  mid_prices
}  # end scrub_agg



### scrub and return a single day of TAQ data
scrub_TAQ <- function(taq_data, agg_vol_window=51, suspect_threshold=2) {

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
  sus_pect <- suspect_bid_offer(bid_offer)
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
  mid_prices[suspect_jump(mid_prices)] <- NA
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
load_data <- function(sym_bol) {
  
  # create path to directory with *.RData files
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
  data <- sapply(data, scrub_agg)
  
  # recursively "rbind" the list into a single xts
  data <- do_call_rbind(data)
  
  colnames(data) <- sapply(strsplit(colnames(data), split="[.]"), 
                           function(strng) paste(sym_bol, strng[-1], sep="."))
  
  assign(sym_bol, data)
  
  save(list=eval(sym_bol), file=paste0(sym_bol, ".RData"))
  
}  # end load_data

