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