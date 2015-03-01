################################################
###
### Demos for managing high frequency data using package 'HighFreq'
###
################################################

rm(list=ls())  # remove all objects

library(HighFreq)

Sys.setenv(TZ="America/New_York")  # Set the time-zone to GMT
setwd("C:/Develop/data")
# search()  # get search path
options(digits.secs=6)
options(digits=5)
options(stringsAsFactors=FALSE)
options(max.print=80)


# set data directories
data_dir <- "E:/mktdata/sec/"
output_dir <- "E:/output/data/"



###########
# load list of instruments

# load list of symbols
sym_bols <- read.csv(file="etf_list_hf.csv")
sym_bols <- sym_bols[[1]]

# define sym_bol
sym_bol <- "SPY"

# load a single day of TAQ data
load(file.path(data_dir, paste0(sym_bol, "/2014.05.02.", sym_bol, ".RData")))
ts_data <- scrub_agg(taq_data=get(sym_bol))
# methods(as.quantmod.OHLC)
# ts_data <- as.quantmod.OHLC(x=ts_data, col.names=c("Open", "High", "Low", "Close", "Volume"))
# is.OHLC(ts_data)
chartSeries(ts_data, name=sym_bol, theme=chartTheme("white"))


###########
# process TAQ data "by hand"

# create path to directory with *.RData files
file_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
file_list <- list.files(file_dir)
# create paths to *.RData files
file_names <- file.path(file_dir, file_list)

# load data into list
ts_data <- sapply(tail(file_names), function(file_name) {
  cat("loading", file_name, "\n")
  data_name <- load(file_name)
  get(data_name)
})
length(ts_data)

# scrub and aggregate the data
ts_data <- lapply(ts_data, scrub_agg)

# flatten list into xts - blows up or takes very long!!!
# ts_data <- do.call(rbind, ts_data)
# recursively "rbind" the list into a single xts
ts_data <- do_call_rbind(ts_data)


# rename the colnames
# colnames(ts_data) <- lapply(strsplit(colnames(ts_data), split="[.]"), 
#                               function(strng) paste(sym_bol, strng[-1], sep="."))
ts_data <- quantmod.OHLC(ts_data)

head(ts_data)
chartSeries(ts_data, name=sym_bol, theme=chartTheme("white"))
chartSeries(ts_data["2008-01-04/2008-01-06"], name=sym_bol, theme=chartTheme("white"))


###########
# process TAQ data using package 'HighFreq': load TAQ data, aggregate to OHLC, and save to file

# process TAQ data for a single symbol
save_OHLC(sym_bol, data_dir=data_dir, output_dir=output_dir, period="15 min")

# process data for list of symbols
sapply(head(sym_bols), save_OHLC, data_dir=data_dir, period="15 min")

# load processed OHLC data for a single symbol
# load(file=paste0(sym_bol, ".RData"))
load(file.path(output_dir, paste0(sym_bol, ".RData")))
# load(file="SPY.RData")
# plot OHLC data
chartSeries(get(sym_bol), name=sym_bol, theme=chartTheme("white"))






########### ignore everything below ###########


###########
# process TAQ data "by hand"

# load one day of TAQ data
load(file.path(data_dir, paste0(sym_bol, "/2014.05.02.", sym_bol, ".RData")))
# load(file.path(data_dir, "SPY/2014.05.02.SPY.RData"))
head(get(sym_bol))


# create path to directory with *.RData files
file_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
file_list <- list.files(file_dir)
# create paths to *.RData files
file_names <- file.path(file_dir, file_list)

# load last six days of data into list
ts_data <- sapply(tail(file_names), function(file_name) {
  cat("loading", file_name, "\n")
  data_name <- load(file_name)
  get(data_name)
})
length(ts_data)

# convert timezone of index to New_York
index(ts_data) <- with_tz(index(ts_data), "America/New_York")
# subset data to NYSE trading hours
daily_data <- daily_data['T09:30:00/T16:00:00', ]

# scrub and aggregate a single day of data
daily_data <- scrub_agg(taq_data=OHLC_data[[3]])
# calculate mid bid-offer prices
mid_prices <- 0.5 * (daily_prices[, 'Bid.Price'] + daily_prices[, 'Ask.Price'])
mid_prices <- na.omit(mid_prices)
colnames(mid_prices) <- "Mid.Price"
chartSeries(daily_data, name=sym_bol, theme=chartTheme("white"))

# scrub and aggregate the data
OHLC_data <- lapply(ts_data, scrub_agg)

# recursively "rbind" the list into a single xts
OHLC_data <- do_call_rbind(OHLC_data)

# rename the colnames
colnames(OHLC_data) <- sapply(strsplit(colnames(OHLC_data), split="[.]"), 
                                 function(strng) paste(sym_bol, strng[-1], sep="."))

head(OHLC_data)
chartSeries(OHLC_data, name=sym_bol, theme=chartTheme("white"))
chartSeries(OHLC_data["2008-01-04/2008-01-06"], name=sym_bol, theme=chartTheme("white"))



### create test data

# create xts time series
x_ts <- xts(x=rnorm(100), order.by=(Sys.time()-3600*(1:100)))
# split time series into daily list
list_xts <- split(x_ts, "days")
# rbind the list back into a time series and compare with the original
identical(x_ts, do_call_rbind(list_xts))


# create time index of one second intervals for a single day
in_dex <- seq(from=as.POSIXct("2015-02-09 09:30:00"), to=as.POSIXct("2015-02-09 16:00:00"), by="1 sec")
# create xts of random prices
x_ts <- xts(cumsum(rnorm(length(in_dex))), order.by=in_dex)
# create vector of random bid-offer prices
bid_offer <- abs(rnorm(length(in_dex)))/10
# create TAQ data using cbind
ts_data <- cbind(x_ts-bid_offer, x_ts+bid_offer)
# add Trade.Price
ts_data <- cbind(ts_data, x_ts+rnorm(length(in_dex))/10)
# add Volume
ts_data <- cbind(ts_data, sample(x=10*(2:18), size=length(in_dex), replace=TRUE))
colnames(ts_data) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
# aggregate to one minute OHLC data
ohlc_data <- scrub_agg(ts_data)
chartSeries(ohlc_data, name=sym_bol, theme=chartTheme("white"))


# mid_prices <- xts(cumsum(rnorm(nrow(mid_prices))), order.by=index(daily_prices))
# sine function with jumps
mid_prices <- xts(sin(22*(1:nrow(daily_prices))/nrow(daily_prices)), order.by=index(daily_prices)) + 2
# mid_prices <- sin(22*(1:dim(daily_prices)[1])/dim(daily_prices)[1])
# prices_scrub <- mid_prices
colnames(mid_prices) <- "Mid.Price"

# add noise
mid_prices[c(1000,3000,5000,7000)] <- 1.1*mid_prices[c(1000,3000,5000,7000)]
mid_prices[c(2000,4000,6000,8000)] <- 0.8*mid_prices[c(1000,3000,5000,7000)]
# diff_prices <- xts(diff_prices, order.by=index(OHLC_data[[6]])[10001:20000])
plot(mid_prices)


# median filter
test.blob <- xts(runmed(x=coredata(mid_prices), 11), order.by=index(mid_prices))
plot(mid_prices, xlab="", ylab="", type='l')
lines(test.blob, col='red', lwd=1)


# calculate stddev, skewness, and quantiles
sd(x=coredata(daily_agg_returns))
skewness(x=coredata(daily_agg_returns))
quantile(x=daily_agg_returns, probs=c(0.05, 0.95))
quantile(x=daily_agg_returns, probs=c(0.1, 0.9))


# plot histograms of daily returns
hist(daily_agg_returns, breaks=200, main="returns", xlab="", ylab="", freq=FALSE)
lines(density(daily_agg_returns), col='red', lwd=1)  # draw density

hist(returns_scrub, breaks=200, main="returns", xlab="", ylab="", freq=FALSE)
lines(density(returns_scrub), col='red', lwd=1)  # draw density

hist(returns_scrub, breaks=300, main="returns", xlab="", ylab="", xlim=c(-0.05, 0.05), freq=FALSE)
lines(density(returns_scrub), col='red', lwd=1)  # draw density

hist(daily_returns, breaks=100, main="returns", xlim=c(-2.0e-4, 2.0e-4), ylim=c(0, 10000), xlab="", ylab="", freq=FALSE)
lines(density(daily_returns), col='red', lwd=1)  # draw density

# title(main=ch.title, line=-1)  # add title

