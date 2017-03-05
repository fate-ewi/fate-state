setwd("/Users/MikeLitzow/Documents/R/FATE/fate-ewi")
library(ncdf4)
require(chron)

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


SLP <- ncvar_get(nc, "slp", start=c(77,11,13), count=c(19,8,length(d)), verbose = T)