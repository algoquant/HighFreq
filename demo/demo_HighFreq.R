################################################
###
### Demos for managing high frequency data using package 'HighFreq'
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
install.packages("devtools")
library(devtools)
install_github(repo="algoquant/HighFreq")

# install HighFreq using install.packages()
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
# process single day of data

# load list of symbols from file in cwd
# sym_bols <- read.csv(file="etf_list_hf.csv")
# sym_bols <- sym_bols[[1]]
# the file "hf_data.RData" is part of "HighFreq" package, and contains "sym_bols"
data("hf_data")

# define sym_bol
sym_bol <- "SPY"
inter_val <- "2013-11-11/2013-11-15"

# load a single day of TAQ data
sym_bol <- load(
  file.path(data_dir, 
            paste0(sym_bol, "/2014.05.02.", sym_bol, ".RData")))

### scrub a single day of TAQ data (don't aggregate)
taq_data <- scrub_TAQ(taq_data=get(sym_bol))

# calculate returns
xts_rets <- calc_rets(xts_data=taq_data)

### scrub and aggregate a single day of TAQ data to OHLC
ohlc_data <- scrub_agg(taq_data=get(sym_bol))
chartSeries(ohlc_data, name=sym_bol, theme=chartTheme("white"))


###########
# process TAQ data using package 'HighFreq': load TAQ data, aggregate to OHLC, and save to file

# aggregate TAQ data for a single symbol, and save to file
save_scrub_agg(sym_bol, data_dir=data_dir, output_dir=output_dir, period="15 min")

save_TAQ(sym_bol, data_dir=data_dir, output_dir=output_dir)

# calculate returns for a single symbol, and save to file
save_rets(sym_bol, data_dir=data_dir, output_dir=output_dir, period="15 min")

# aggregate data for list of symbols, and save to multiple files
sapply(head(sym_bols), save_scrub_agg, data_dir=data_dir, output_dir=output_dir, period="15 min")

# calculate returns for list of symbols, and save to file
sapply(head(sym_bols), save_rets, data_dir=data_dir, output_dir=output_dir, period="15 min")

# load processed OHLC data for a single symbol
# load(file=paste0(sym_bol, ".RData"))
load(file.path(output_dir, paste0(sym_bol, ".RData")))
# load(file="SPY.RData")
# plot OHLC data
chartSeries(get(sym_bol), name=sym_bol, theme=chartTheme("white"))


###########
### estimating volatility, skewness, and kurtosis using package TTR

# function volatility() from package TTR
vol_at <- volatility(OHLC=get(sym_bol)[inter_val], 
                     calc="yang.zhang", n=20)
vol_at <- volatility(OHLC=get(sym_bol)[inter_val], 
                     calc="rogers.satchell", n=20)
vol_at <- volatility(OHLC=get(sym_bol)[inter_val], 
                     calc="garman.klass", n=20)


# estimating rolling aggregations and moments using package TTR
vol_at <- runMedian(x=get(sym_bol)[inter_val], n=100)
vol_at <- runSD(x=get(sym_bol)[inter_val], n=100)
chart_xts(vol_at)


###########
# estimating rolling moments using package HighFreq
library(HighFreq)

# daily open to close variance and skew
var_iance <- apply.daily(x=get(sym_bol), FUN=moment_ohlc)
colnames(var_iance) <- paste(
  strsplit(colnames(get(sym_bol))[1], split="[.]")[[1]][1], 
  "Var", sep=".")
sk_ew <- apply.daily(x=get(sym_bol), FUN=moment_ohlc, mom_fun="skew_ohlc")
sk_ew <- sk_ew/(var_iance)^(1.5)
colnames(sk_ew) <- paste(
  strsplit(colnames(get(sym_bol))[1], split="[.]")[[1]][1], 
  "Skew", sep=".")

re_turns <- calc_rets(xts_data=get(sym_bol))
daily_autocorr <- re_turns[, 1]*lag(re_turns[, 1])
colnames(daily_autocorr) <- paste0(sym_bol, ".autocorr")
daily_autocorr[1, ] <- 0
daily_autocorr <- apply.daily(x=daily_autocorr, FUN=sum)
sum(daily_autocorr)



x11()
inter_val <- "2013-06-01/"
chart_Series(var_iance[inter_val], name=paste(sym_bol, "variance"))


chart_Series(da_ta[inter_val, 1]/(ran_ge[2]-ran_ge[1]), name=paste(sym_bol, "data"))

# chart_xts_panels(cbind(10^5*var_iance, sk_ew))
# chart_xts_panels(cbind(10^5*var_iance["2013-10"], sk_ew["2013-10"]), in_dex=(abs(sk_ew["2013-10"])>1))


# daily close to close volatility



# minutely variance and skew
# minutely returns
re_turns <- calc_rets(xts_data=get(sym_bol)[inter_val])
# rolling volume-weighted returns
roll_returns <- runSum(x=re_turns[, 1], n=20)/20
colnames(roll_returns) <- paste0(sym_bol, ".Rets")
roll_returns[1:20, ] <- 0
foo <- roll_returns[end_points, ]

# rolling volume-weighted variance
var_iance <- roll_moment_ohlc(ohlc=get(sym_bol)[inter_val])
tail(var_iance, 11)
roll_variance <- var_iance[end_points, ]

# rolling volume-weighted autocorrelation
roll_autocorr <- re_turns[, 1]*lag(re_turns[, 1])
roll_autocorr[1, ] <- 0
sum(roll_autocorr)
roll_autocorr <- runSum(x=roll_autocorr[, 1], n=20)/20
roll_autocorr[1:20, ] <- 0
colnames(roll_autocorr) <- paste0(sym_bol, ".autocorr")
roll_autocorr <- roll_autocorr[end_points, ]

# find periods of high variance


da_ta <- cbind(roll_autocorr, roll_variance)

da_ta <- cbind(daily_autocorr, var_iance)
da_ta <- cbind(daily_autocorr, sk_ew)

da_ta <- cbind(foo, roll_variance)

for_mula <- as.formula(paste(colnames(da_ta)[1], 
                             paste(paste(colnames(da_ta)[-1], 
                                         collapse=" + ")), sep="~"))
for_mula <- as.formula(paste(colnames(da_ta)[1], 
                             paste(paste(colnames(da_ta)[-1], 
                                         collapse=" + "), "- 1"), sep="~"))
for_mula
l_m <- lm(for_mula, data=da_ta)
summary(l_m)
# scatterplot of da_ta
x11()
plot(for_mula, data=da_ta)
# plot of da_ta
x11()
inter_val <- "2014-04-07/2014-04-10"
ran_ge <- range(da_ta[inter_val, 1])
chart_Series(da_ta[inter_val, 1]/(ran_ge[2]-ran_ge[1]), name=paste(sym_bol, "data"))
ran_ge <- range(da_ta[inter_val, 2])
add_TA(da_ta[inter_val, 2]/(ran_ge[2]-ran_ge[1]), on=1, col="blue", lwd=2)

# rolling volume-weighted skewness
sk_ew <- roll_moment_ohlc(ohlc=get(sym_bol)[inter_val], mom_fun="skew_ohlc")
sk_ew <- sk_ew/(var_iance)^(1.5)
sk_ew[1, ] <- 0
sk_ew <- na.locf(sk_ew)
tail(sk_ew, 11)


