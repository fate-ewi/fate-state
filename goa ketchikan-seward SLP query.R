setwd("/Users/MikeLitzow/Documents/R/FATE/fate-ewi")
library(ncdf4)
require(chron)
require(zoo)
require(dplyr)

# dowload the latest version of the NCEP/NCAR Reanalysis
# download.file(url="ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/slp.mon.mean.nc", destfile= "/Users/MikeLitzow/Documents/R/NSF-GOA/slp.mon.mean.nc")

# now load the file
nc <- nc_open("/Users/MikeLitzow/Documents/R/NSF-GOA/slp.mon.mean.nc")

# get date info - time is in hours since Jan. 1, 1800!
raw <- ncvar_get(nc, "time")
h <- raw/24
d <- dates(h, origin = c(1,1,1800))
d # we'll just use the entire time series - Jan 1948-present

# examin lat/lons
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

# Ketchikan=N55, E230; Seward=N60, E210

# check Ketchikan coordinates - look right
Kx <- ncvar_get(nc, "lon", start=93, count=1)
Ky <- ncvar_get(nc, "lat", start=15, count=1)
Kx;Ky

# and check Seward coordinates
Sx <- ncvar_get(nc, "lon", start=85, count=1)
Sy <- ncvar_get(nc, "lat", start=13, count=1)
Sx;Sy

Ketchikan <- ncvar_get(nc, "slp", start=c(93,15,1), count=c(1,1,length(d)), verbose = T)
Seward <- ncvar_get(nc, "slp", start=c(85,13,1), count=c(1,1,length(d)), verbose = T)
# and we're interested in the difference between the two
diff <- Ketchikan-Seward
names(diff) <- d # label by date

# quick plot to check
plot(d, diff, type="l")
lines(d, rollmean(diff, k=37, fill=NA), col="red")
   
      abline(h=mean(diff))

# now summarize
yrs <- years(d)
m <- months(d)

# I reckon we should use monthly anomalies to remove the seasonal signal
# this should produce a more informative annual mean
mu <- tapply(diff, m, mean)
sd <- tapply(diff, m, sd)


mu <- rep(mu, length.out=length(d))
sd <- rep(sd, length.out=length(d))
diff.anom <- (diff-mu)/sd

# first, annual means
annual <- tapply(diff, yrs, mean)
plot(names(annual), annual, type="o")

# limit to complete cases
f <- function(x) sum(!is.na(x)) # function to count number of observations
use <- tapply(diff, yrs, f)==12
annual <- annual[use]

# now winter
win <- c("Nov", "Dec", "Jan", "Feb", "Mar") # define winter months
win.yrs <- as.numeric(as.character(yrs)) # as 

# align winter years with year corresponding to January
early.winter <- c("Nov", "Dec")
add <- m %in% early.winter
win.yrs[add] <- win.yrs[add]+1

use <- m %in% win

winter.diff <- diff[use]
win.yrs <- win.yrs[use]

win.mean <- tapply(winter.diff, win.yrs, mean)

# and limit to complete cases
use <- tapply(winter.diff, win.yrs, f)==5
win.mean <- win.mean[use]
# and add a leading NA to fill in 1948
win.mean <- c(NA, win.mean)
names(win.mean)[1] <- 1948

# check with a plot
plot(names(win.mean), win.mean, type="o")

# now spring
spr <- c("Apr", "May", "Jun") # define spring months

use <- m %in% spr

spr.diff <- diff[use]
spr.yrs <- yrs[use]

spr.mean <- tapply(spr.diff, spr.yrs, mean)

# and limit to complete cases
use <- tapply(spr.diff, spr.yrs, f)==3
spr.mean <- spr.mean[use]


# check with a plot
plot(names(spr.mean), spr.mean, type="o")
