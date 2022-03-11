################################################
###
### Demos for managing high frequency data using package HighFreq
###
################################################

# Set the time-zone to New_York
Sys.setenv(TZ="America/New_York")
# setwd("C:/Develop/data")
# search()  # get search path
options(digits.secs=6)
options(digits=7)
# suppress spurious timezone warning messages
options(xts_check_TZ=FALSE)
options(stringsAsFactors=FALSE)
options(max.print=80)
rm(list=ls())  # remove all objects

# install HighFreq from github
# install.packages("devtools")
# library(devtools)
devtools::install_github(repo="algoquant/HighFreq")

# install HighFreq from local directory using install.packages()
install.packages(pkgs="C:/Develop/R/HighFreq", repos=NULL, type="source")
install.packages(pkgs="C:/Develop/R/HighFreq", repos=NULL, type="source", lib="C:/Users/Jerzy/Downloads")
install.packages(pkgs="C:/Develop/R/HighFreq", repos=NULL, type="source", lib="C:/Users/Jerzy/Documents/R/win-library/3.2")

# load HighFreq and attach the data
library(HighFreq)
# data(hf_data)

# set data directories
data_dir <- "E:/mktdata/sec/"
output_dir <- "E:/output/data/"

# set data directories
data_dir <- "C:/Develop/data/hfreq/src/"
output_dir <- "C:/Develop/data/hfreq/scrub/"


###########
### process single day of high frequency data

# load list of symbols from file in cwd
# symbolv <- read.csv(file="etf_list_hf.csv")
# symbolv <- symbolv[[1]]
# the file "hf_data.RData" is part of "HighFreq" package, and contains "symbolv"
data("hf_data")

# define symbol
symbol <- "SPY"
interval <- "2013-11-11/2013-11-15"

# load a single day of TAQ data
symbol <- load(file.path(data_dir, 
            paste0(symbol, "/2014.05.02.", symbol, ".RData")))

### scrub a single day of TAQ data (don't aggregate)
taq <- scrub_taq(taq=SPY)

# calculate returns from the secondly TAQ data
returns_running <- run_returns(xtes=taq)

### scrub and aggregate a single day of TAQ data to OHLC
ohlc <- scrub_agg(taq=taq)
chart_Series(ohlc, name=symbol)


###########
### process TAQ data using package HighFreq: load TAQ data, aggregate to OHLC, and save to file

# aggregate TAQ data to 15-min OHLC bar data, for a single symbol, and save to file
save_scrub_agg(symbol, data_dir=data_dir, output_dir=output_dir, period="15 min")

save_taq(symbol, data_dir=data_dir, output_dir=output_dir)

# calculate returns for a single symbol, and save to file
save_rets(symbol, data_dir=data_dir, output_dir=output_dir, period="15 min")

# aggregate data for list of symbols, and save to multiple files
sapply(head(symbolv), save_scrub_agg, data_dir=data_dir, output_dir=output_dir, period="15 min")

# calculate returns for list of symbols, and save to file
sapply(head(symbolv), save_rets, data_dir=data_dir, output_dir=output_dir, period="15 min")

# load processed OHLC data for a single symbol
# load(file=paste0(symbol, ".RData"))
load(file.path(output_dir, paste0(symbol, ".RData")))
# load(file="SPY.RData")
# plot OHLC data
chart_Series(SPY, name=symbol)


###########
### estimating volatility, skewness, and kurtosis using package TTR

# function volatility() from package TTR
volat <- volatility(OHLC=SPY[interval], 
                     calc="yang_zhang", n=10)
volat <- volatility(OHLC=SPY[interval], 
                     calc="rogers_satchell", n=10)
volat <- volatility(OHLC=SPY[interval], 
                     calc="garman_klass", n=10)


# estimating rolling aggregations and moments using package TTR
volat <- runMedian(x=SPY[interval], n=100)
volat <- runSD(x=SPY[interval], n=100)
chart_xts(volat)


###########
### estimating rolling moments using package HighFreq
library(HighFreq)

look_back <- 10

### daily open to close variance and skew
volume_daily <- xts::apply.daily(x=Vo(SPY), FUN=sum)
colnames(volume_daily) <- paste0(symbol, ".Volume")
var_daily <- (6.5*60*60)*xts::apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="roll_variance", method="rogers_satchell")
colnames(var_daily) <- paste0(symbol, ".Var")
skew_daily <- xts::apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="run_skew")
skew_daily <- skew_daily/(var_daily)^(1.5)
colnames(skew_daily) <- paste0(symbol, ".Skew")

# daily Sharpe
sharpe_daily <- xts::apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="run_sharpe")
colnames(sharpe_daily) <- paste0(symbol, ".Sharpe")
chart_Series(sharpe_daily[interval], name=paste(symbol, "Sharpe"))

# simple autocorrelation
returns_running <- run_returns(xtes=SPY)
hurst_daily <- returns_running*lag(returns_running)/returns_running^2
hurst_daily[1, ] <- 0
hurst_daily[is.nan(hurst_daily), ] <- 0
hurst_daily <- xts::apply.daily(x=hurst_daily, FUN=sum)

library(PerformanceAnalytics)
hurst_daily <- xts::apply.daily(x=hurst_daily, FUN=HurstIndex)

hurst_daily <- xts::apply.daily(x=hurst_daily, FUN=sum)

colnames(hurst_daily) <- paste0(symbol, ".autocorr")

x11(width=6, height=4)
interval <- "2013-06-01/"
chart_Series(var_daily[interval], name=paste(symbol, "variance"))
chart_Series(rutils::roll_sum(hurst_daily, look_back=10)[-(1:10)]/10, name=paste(symbol, "Hurst"))
abline(h=0.5, col="blue", lwd=2)

# calculate minutely returns of OHLC data
returns_running <- run_returns(xtes=SPY[interval])

# rolling over endpoints
n_row <- NROW(returns_running)
num_agg <- n_row %/% 10
endpoints <- c(0, n_row-10*num_agg + 10*(0:num_agg))

# rolling average prices
look_back <- 10
prices_rolling <- rutils::roll_sum(SPY[, 1], look_back=look_back)/look_back
colnames(prices_rolling) <- paste0(symbol, ".Prices")
chart_Series(SPY["2013-11-12"], name=paste(symbol, "Prices"))
add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)
legend("top", legend=c("SPY prices", "average prices"), 
       bg="white", lty=c(1, 1), lwd=c(2, 2), 
       col=c("black", "red"), bty="n")


### rolling volume-weighted returns
returns_rolling <- roll_vwap(ohlc=SPY, xtes=returns_running, look_back=look_back)
returns_rolling <- returns_rolling[endpoints, ]

### calculate rolling SPY variance without overnight jumps
var_rolling <- (6.5*60*60)*roll_variance(ohlc=SPY, method="rogers_satchell")
# calculate rolling SPY variance
# var_rolling <- roll_vwap(ohlc=SPY, xtes=var_running, look_back=look_back)
var_rolling <- var_rolling[endpoints, ]

### calculate running and rolling SPY skew
skew_running <- (6.5*60*60)*run_skew(ohlc=SPY)
skew_rolling <- roll_vwap(ohlc=SPY, xtes=skew_running, look_back=look_back)

# calculate intraday seasonality of skew
skewseasonal <- season_ality(skew_running)

interval <- "2014-03-01/"
chart_Series(var_rolling[interval], name=colnames(var_rolling))


### calculate running and rolling Sharpe ratio
sharpe_running <- run_sharpe(ohlc=SPY)
sharpe_rolling <- roll_vwap(ohlc=SPY, xtes=sharpe_running, look_back=look_back)


### calculate rolling Hurst exponents using different methods

# create synthetic OHLC time series of random prices
ohlc <- HighFreq::random_ohlc(indeks=zoo::index(SPY["2009-03"]))

# calculate rolling Hurst for random prices using ratio of price range divided by standard deviation of returns
hurst_rolling <- roll_hurst(ohlc=ohlc, look_back=look_back)

# remove overnight price jumps from SPY in March 2009
ohlc <- remove_jumps(SPY["2009-03"])

# legacy code for identifying overnight price jumps
# find biggest close-to-open price jumps
foo <- log(ohlc[, 1])-rutils::lagxts(log(ohlc[, 4]))
foo <- sort(as.vector(foo))
plot(foo, t="l")
bar <- quantile(foo, probs=c(0.997, 0.003))
bar <- which(foo > bar[1] | foo < bar[2])

# calculate rolling Hurst for SPY in March 2009
hurst_rolling <- roll_hurst(ohlc=ohlc, look_back=look_back)
hurst_rolling <- roll_hurst(ohlc=SPY["2009-03"], look_back=look_back)
chart_Series(hurst_rolling["2009-03-10/2009-03-12"], name=paste("SPY", "hurst_rolling"))

# old version of roll_hurst(): calculate Hurst over end points of SPY
hurst_rolling <- roll_hurst(ohlc=SPY, look_back=11, off_set=0, roll_endpoints=TRUE)
# calculate a series of rolling hurst values using argument off_set
hurst_rolling <- lapply(0:9, roll_hurst, ohlc=SPY, look_back=11, roll_endpoints=TRUE)
hurst_rolling <- rutils::do_call_rbind(hurst_rolling)
# remove daily warmup periods
hurst_rolling <- hurst_rolling["T09:41:00/T16:00:00"]

# calculate rolling Hurst using ratio of range variance estimators
# old version of roll_hurst(): calculate hurst_rolling by applying roll_hurst() over argument off_set
hurst_rolling <- lapply(0:(look_back-1), roll_hurst, ohlc=SPY, look_back=look_back)
hurst_rolling <- rutils::do_call_rbind(hurst_rolling)
# merge and remove duplicate trailing values
duplicatesd <- hurst_rolling[duplicated(index(hurst_rolling)), ]
duplicatess <- rutils::do_call_rbind(lapply(unique(index(duplicatesd)), 
                                            function(duplicates) 
                                              xts(mean(hurst_rolling[duplicates, ]), 
                                                  order.by=duplicates)))
duplicatesd <- index(hurst_rolling) == unique(index(duplicatesd))
hurst_rolling <- hurst_rolling[!duplicatesd, ]
hurst_rolling <- rbind(hurst_rolling, duplicatess)
tail(hurst_rolling, 12)
dim(hurst_rolling)

# remove daily warmup periods
hurst_rolling <- hurst_rolling["T09:41:00/T16:00:00"]

# remove stub periods with low frequency of observations
# calculate frequency of time periods
period_freq <- table(indeks)
period_freq <- structure(as.vector(period_freq), names=names(period_freq))
# remove stub periods
period_freq <- period_freq[period_freq > mean(period_freq)]
period_freq <- period_freq[-length(period_freq)]
hurst_rolling <- hurst_rolling[indeks %in% names(period_freq), ]
chart_Series(x=hurst_rolling["2012-02-13"], 
             name=paste(colnames(hurst_rolling), "10-minute aggregations"))

### rolling autocorrelation
roll_autocorr <- returns_running*lag(returns_running)
roll_autocorr[1, ] <- 0
roll_autocorr <- rutils::roll_sum(roll_autocorr[, 1], look_back=10)/10
roll_autocorr[1:10, ] <- 0
colnames(roll_autocorr) <- paste0(symbol, ".autocorr")
roll_autocorr <- roll_autocorr[endpoints, ]


### design matrix called SPY_design containing columns of aggregations

look_back <- 5
returns_running <- run_returns(xtes=SPY)
returns_rolling <- roll_vwap(ohlc=SPY, xtes=returns_running, look_back=look_back)
colnames(returns_running) <- "returns"
colnames(returns_rolling) <- "returns.WA5"

var_running <- (6.5*60*60)*run_variance(ohlc=SPY)
# var_running <- run_variance(ohlc=SPY, method="rogers_satchell")
# var_diff <- rutils::diffxts(xtes=var_running)
var_rolling <- roll_vwap(ohlc=SPY, xtes=var_running, look_back=look_back)
colnames(var_running) <- "variance"
# colnames(var_diff) <- "variance.diff"
colnames(var_rolling) <- "variance.WA5"

skew_running <- (6.5*60*60)*run_skew(ohlc=SPY)
# skew_diff <- rutils::diffxts(xtes=skew_running)
skew_rolling <- roll_vwap(ohlc=SPY, xtes=skew_running, look_back=look_back)
colnames(skew_running) <- "skew"
# colnames(skew_diff) <- "skew.diff"
colnames(skew_rolling) <- "skew.WA5"

sharpe_running <- run_sharpe(ohlc=SPY)
sharpe_rolling <- roll_vwap(ohlc=SPY, xtes=sharpe_running, look_back=look_back)
colnames(sharpe_running) <- "sharpe_running"
colnames(sharpe_rolling) <- "sharpe_running.WA5"

sharpe_rolling <- roll_sharpe(ohlc=SPY, look_back=look_back)
# sharpe_diff <- rutils::diffxts(xtes=sharpe_rolling)
colnames(sharpe_rolling) <- "sharpe_rolling"
# colnames(sharpe_diff) <- "sharpe.diff"
# x11()
# chart_Series(sharpe_rolling["2013-11-12", ], name=paste("SPY", "sharpe_rolling"))

hurst_rolling <- roll_hurst(ohlc=SPY, look_back=look_back)
# hurst_diff <- rutils::diffxts(xtes=hurst_rolling)
colnames(hurst_rolling) <- "hurst_rolling"
# colnames(hurst_diff) <- "hurst.diff"
chart_Series(hurst_rolling["2009-11-12", ], name=paste("SPY", "hurst_rolling"))

# vol_diff <- rutils::diffxts(xtes=Vo(SPY))
# colnames(vol_diff) <- "volume.diff"
# SPY_design <- cbind(SPY_design, Vo(SPY), vol_diff)

SPY_design <- cbind(returns_running, returns_rolling, var_running, 
                    var_rolling, skew_running, sharpe_running, 
                    sharpe_rolling, hurst_rolling)
colnames(SPY_design)
# SPY_design <- cbind(SPY_design[, -1], SPY_design[, "SPY.Volume"])

# select most significant factors plus interaction terms
SPY_design <- cbind(returns_running, returns_rolling, var_running, skew_running, 
                    hurst_rolling, returns_running*var_running, returns_running*skew_running)
colnames(SPY_design) <- c(colnames(SPY_design)[1:4], "hurst", "rets_var", "rets_skew")

head(SPY_design["2013-11-12", ])
save(SPY_design, file="C:/Develop/data/SPY_design.RData")
load("C:/Develop/data/SPY_design.RData")


### intraday seasonality

# calculate intraday seasonality of volume
volume_seasonal <- season_ality(Vo(SPY[interval]))
colnames(volume_seasonal) <- paste0(symbol, ".volume_seasonal")
seasondata <- volume_seasonal

# calculate intraday seasonality of variance
var_seasonal <- season_ality((6.5*60*60)*run_variance(ohlc=SPY[interval, 1:4], method="rogers_satchell"))
colnames(var_seasonal) <- paste0(symbol, ".var_seasonal")
seasondata <- var_seasonal

# calculate intraday seasonality of variance
# calculate minutely variance of OHLC data scaled to daily units
xtes <- (6.5*60*60)*run_variance(SPY)
# remove overnight variance spikes at "09:31"
xtes <- xtes["T09:32:00/T16:00:00"]
# calculate intraday seasonality of variance
var_seasonal <- season_ality(xtes=xtes)
chart_Series(x=var_seasonal, 
             name=paste(colnames(var_seasonal), "intraday seasonality"))

season_illiquid <- 1e6*sqrt(var_seasonal/volume_seasonal)
foo <- 1e6*ifelse(Vo(SPY)==0, 0, sqrt(run_variance(ohlc=SPY)/Vo(SPY)))
season_illiquid <- season_ality(foo)
season_illiquid <- season_ality(sqrt(run_variance(ohlc=SPY)/Vo(SPY)))
colnames(season_illiquid) <- paste0(symbol, ".season_illiquid")
seasondata <- season_illiquid

# calculate intraday seasonality of skew
skewseasonal <- season_ality(run_skew(ohlc=SPY[interval, 1:4], method="rogers_satchell"))
# skewseasonal <- skewseasonal/(var_seasonal)^(1.5)
colnames(skewseasonal) <- paste0(symbol, ".skewseasonal")
seasondata <- skewseasonal

# calculate intraday seasonality of Hurst exponent
hurst_seasonal <- season_ality(hurst_rolling)
colnames(hurst_seasonal) <- paste0(colnames(hurst_rolling), ".seasonal")
# plot without daily warmup period
chart_Series(x=hurst_seasonal[-(1:10), ], 
             name=paste(colnames(hurst_seasonal), "intraday seasonality"))
# below is for run_hurst() which isn't really true Hurst
hurst_seasonal <- season_ality(run_hurst(ohlc=SPY[interval, 1:4]))
hurst_seasonal <- hurst_seasonal[-NROW(hurst_seasonal)]
colnames(hurst_seasonal) <- paste0(symbol, ".seasonal")
seasondata <- hurst_seasonal

season_autocorr <- returns_running*lag(returns_running)
season_autocorr[1, ] <- 0
season_autocorr <- returns_running*lag(rutils::roll_sum(returns_running, look_back=5)/5)
season_autocorr[1:5, ] <- 0
season_autocorr <- season_ality(season_autocorr)
colnames(season_autocorr) <- paste0(symbol, ".season_autocorr")
seasondata <- rutils::roll_sum(season_autocorr, look_back=5)/5

# daily Hurst exponents
hurst_daily <- xts::apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="run_hurst")
hurst_daily <- xts::apply.daily(x=SPY, FUN=function(x, ...) abs(agg_stats_r(ohlc=x, ...)), calc_bars="run_hurst")
colnames(hurst_daily) <- paste0(symbol, ".Hurst.daily")
chart_Series(roll_sum(hurst_daily, 10)[-(1:10)]/10, name=paste(symbol, "Hurst"))
abline(h=0.5, col="blue", lwd=2)


## aggregate SPY to 10-minute bars
SPY10min <- to.minutes10(SPY)
colnames(SPY10min) <- colnames(SPY)
head(SPY10min)

## plot
seasondata <- seasondata[-(NROW(seasondata))]
x11(width=6, height=4)
chart_Series(x=seasondata, 
             name=paste(colnames(seasondata), "intraday seasonality"))
# or, adjust y-axis range
rutils::chart_xts(seasondata,
          name=paste(colnames(seasondata), "intraday seasonality"), 
          ylim=c(min(seasondata), max(seasondata)/2))
# or, adjust y-axis range by hand
plot_theme <- chart_theme()
plot_theme$format.labels <- "%H:%M"
chobj <- chart_Series(x=seasondata, 
                      name=paste(colnames(seasondata), "intraday seasonality"), 
                      theme=plot_theme, plot=FALSE)
ylim <- chobj$get_ylim()
ylim[[2]] <- structure(c(ylim[[2]][1], ylim[[2]][2]/2), fixed=TRUE)
chobj$set_ylim(ylim)
plot(chobj)


## regressions
datav <- cbind(roll_autocorr, var_rolling)
datav <- cbind(foo, var_rolling)

datav <- cbind(hurst_daily, var_daily)
datav <- cbind(hurst_seasonal, var_seasonal)
datav <- cbind(hurst_daily, skew_daily)

formulav <- as.formula(paste(colnames(datav)[1], 
                             paste(paste(colnames(datav)[-1], 
                                         collapse=" + ")), sep="~"))
formulav <- as.formula(paste(colnames(datav)[1], 
                             paste(paste(colnames(datav)[-1], 
                                         collapse=" + "), "- 1"), sep="~"))
l_m <- lm(formulav, data=datav)
summary(l_m)
# scatterplot of datav
x11()
plot(formulav, data=datav)
# plot of datav
x11()
interval <- "2014-04-07/2014-04-10"
rangev <- range(datav[interval, 1])
chart_Series(datav[interval, 1]/(rangev[2]-rangev[1]), name=paste(symbol, "data"))
rangev <- range(datav[interval, 2])
add_TA(datav[interval, 2]/(rangev[2]-rangev[1]), on=1, col="blue", lwd=2)

# rolling volume-weighted skewness
skew_rolling <- roll_vwap(ohlc=SPY[interval], xtes=run_skew(ohlc=SPY[interval]), look_back=10)
skew_rolling <- skew_rolling/(var_rolling)^(1.5)
skew_rolling[1, ] <- 0
skew_rolling <- xts:::na.locf.xts(skew_rolling)
tail(skew_rolling, 11)


## Create minutely synthetic OHLC time series of random prices

set.seed(1121)  # reset random number generator

# use function random_ohlc()
# ohlc <- HighFreq::random_ohlc()

# create minutely synthetic OHLC time series of random prices by hand
# create time index of one second intervals
indeks <- seq(from=as.POSIXct("2016-01-01 00:00:00"),
              to=as.POSIXct("2016-01-30 00:00:00"), by="1 sec")
# create xts of random prices
volat <- 0.001
xtes <- xts(exp(cumsum(rnorm(length(indeks), sd=volat) - volat^2/2)), order.by=indeks)
colnames(xtes) <- "random"
chart_Series(x=xtes["2016-01-10 09/2016-01-10 10"], name="random prices")
# aggregate to minutes OHLC data
ohlc <- xts::to.period(x=xtes, period="minutes", name="random")
chart_Series(x=ohlc["2016-01-10"], name="random OHLC prices")
# add volume
ohlc <- cbind(ohlc, sample(x=10*(2:18), size=NROW(ohlc), replace=TRUE))
colnames(ohlc)[ 5] <- "random.volume"
tail(ohlc)


## calculate variance estimates and their standard errors using bootstrap and function roll_variance()

# calculate variance estimates using function roll_variance()
var_rolling <- (6.5*60*60)*roll_variance(ohlc)
sum(var_rolling)

# calculate variance using all the different estimators in roll_variance()
methods <- c("close", "garman_klass", "rogers_satchell", "garman_klass_yz", "yang_zhang")
sapply(methods, function(method) {
  sum((6.5*60*60)*roll_variance(ohlc, method=method))
})

# calculate secondly returns from TAQ data
returns_running <- run_returns(xtes=SPY_TAQ)
returns_running <- returns_running - sum(returns_running)/NROW(returns_running) - drop(var(returns_running))/2
indeks <- index(SPY_TAQ)

# calculate standard errors of variance estimators using bootstrap
boot_strap <- sapply(1:100, function(x) {
  boot_sample <- sample(returns_running, replace=TRUE)
  xtes <- xts::xts(exp(cumsum(boot_sample)), order.by=indeks)
# aggregate to minutes OHLC data
  ohlc <- xts::to.period(x=xtes, period="minutes")
# calculate variance estimates
  sapply(methods, function(method) {
    sum(HighFreq::roll_variance(ohlc, look_back=look_back, method=method))
  })  # end sapply
})  # end sapply

# calculate standard errors of variance estimators using TTR and bootstrap
methods <- c("close", "garman.klass", "rogers.satchell", "gk.yz", "yang.zhang")
boot_strap <- sapply(1:100, function(x) {
  boot_sample <- sample(returns_running, replace=TRUE)
  xtes <- xts(exp(cumsum(boot_sample)), order.by=indeks)
  # aggregate to minutes OHLC data
  ohlc <- xts::to.period(x=xtes, period="minutes")
  # calculate variance estimates
  sapply(methods, function(method) {
    sum(na.omit(TTR::volatility(ohlc, n=look_back, N=1, calc=method))^2)
  })  # end sapply
})  # end sapply

head(t(boot_strap))
boot_errors <- apply(boot_strap, MARGIN=1, sd)


# calculate standard errors of variance estimators using bootstrap - parallel version under Windows
library(parallel)  # load package parallel
num_cores <- detectCores() - 1  # number of cores
cluster <- makeCluster(num_cores)  # initialize compute cluster

boot_strap <- parLapply(cluster, 1:1000, 
                        function(x, returns, indeks, methods) {
  boot_sample <- sample(returns, replace=TRUE)
  xtes <- xts::xts(exp(cumsum(boot_sample)), order.by=indeks)
  # aggregate to minutes OHLC data
  ohlc <- xts::to.period(x=xtes, period="minutes")
  # calculate variance estimates
  sapply(methods, function(method) {
    sum((6.5*60*60)*HighFreq::roll_variance(ohlc, method=method))
  })  # end sapply
}, returns=returns_running, indeks=indeks, methods=methods)  # end parLapply

# analyze bootstrapped variance
boot_strap <- do.call(rbind, boot_strap)
head(boot_strap)
sum(is.na(boot_strap))
# means and standard errors from bootstrap
apply(boot_strap, MARGIN=2, 
      function(x) c(mean=mean(x), stderror=sd(x)))

save(boot_strap, file="C:/Develop/data/boot_strap.RData")

stopCluster(cluster)

# calculate variance estimates for SPY
var_rolling <- (6.5*60*60)*roll_variance(SPY)


###########
### benchmark HighFreq::roll_apply() to functions xts:::rollapply.xts(), 
# xts::period.apply(), and PerformanceAnalytics::apply.rolling()

# extract a single day of SPY data
prices <- SPY["2012-02-13"]
look_back <- 11
library(PerformanceAnalytics)  # load package PerformanceAnalytics
# define aggregation function that returns a vector of values
agg_function <- function(xtes)  c(max(xtes[, 1]), min(xtes[, 4]))
agg_function <- function(xtes)  sapply(xtes, mean)
agg_function(prices)
# perform aggregations for vector-valued aggregation function
foo <- xts:::rollapply.xts(prices, FUN=agg_function, width=look_back, align="right", by.column=FALSE)
bar <- HighFreq::roll_apply(prices, agg_fun=agg_function, look_back=look_back)
# perform aggregations for scalar-valued aggregation function
foo <- xts:::rollapply.xts(prices, FUN=sum, width=look_back, align="right")
bar <- HighFreq::roll_apply(prices, agg_fun=sum, look_back=look_back, by_columns=TRUE)

# define end points at 10-minute intervals (SPY is minutely bars)
endpoints <- rutils::calc_endpoints(prices, interval=look_back)
foo <- xts::period.apply(prices[, 1], FUN=sum, INDEX=endpoints)
bar <- HighFreq::roll_apply(prices[, 1], agg_fun=sum, look_back=2, endpoints=endpoints)
foo <- xts::period.apply(prices, FUN=agg_function, INDEX=endpoints)
bar <- HighFreq::roll_apply(prices, agg_fun=agg_function, look_back=2, endpoints=endpoints)

# perform aggregations over length of endpoints
foo_bar <- PerformanceAnalytics::apply.rolling(prices, FUN=agg_function, width=look_back)
head(cbind(foo, bar, foo_bar), 14)

# benchmark the speed of the functionals
library(microbenchmark)
summary(microbenchmark(
  roll_apply=HighFreq::roll_apply(prices, agg_fun=sum, look_back=look_back, by_columns=TRUE),
  rollapply=xts:::rollapply.xts(prices, FUN=sum, width=look_back, align="right"), 
  apply.rolling=apply.rolling(prices, FUN=sum, width=look_back), 
  times=10))[, c(1, 4, 5)]

summary(microbenchmark(
  roll_apply=HighFreq::roll_apply(prices, agg_fun=agg_function, look_back=look_back),
  rollapply=xts:::rollapply.xts(prices, FUN=agg_function, width=look_back, align="right", by.column=FALSE), 
  times=10))[, c(1, 4, 5)]


