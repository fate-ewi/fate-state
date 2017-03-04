setwd("/Users/MikeLitzow/Documents/R/FATE/fate-ewi")
library(ncdf4)
require(rerddap)
require(chron)

ed_search(query = "NCEP/NCAR Reanalysis", which = "grid")
out <- info("esrlNcepRe")

uwind <- griddap(out, time = c("1948-01-01", "2017-02-28"), latitude=60.0,
                 longitude=190.0, fields="uwnd")

# I get different versions of this error with rerddap!

# Error: HTTP Status 500 - Query error: For variable=uwnd axis#1=latitude: 
#"[" was expected at or after position=33, not "6".

# just downloading manually from the erddap site for now

dat <- read.csv("ebs.raw.wind.csv")
head(dat)
str(dat)

splt <- strsplit(as.character(dat$time), "T")

date <- matrix(unlist(splt), ncol=2, byrow=TRUE)

dat$date <-chron(dates.=date[,1], format="y-m-d")

# make vectors of year and month
yr <- years(dat$date)
m <- months(dat$date)

# now convert to wind direction

dat$dir <-(180/pi * atan2(-dat$uwnd, -dat$vwnd))+180
range(dat$dir)

# now make columns of 1s and 0s to indicate if the wind is blowing NW or SE

dat$NW <- dat$SE <- 0

for(i in 1:nrow(dat)){
  # i <- 2
  if(dat$dir[i] >=105 & dat$dir[i] <=165) dat$SE[i] <- 1
  if(dat$dir[i] >=285 & dat$dir[i] <=345) dat$NW[i] <- 1
}

# function to calculate the proportion of days with a particular wind direction
f <- function(x) sum(x)/sum(!is.na(x)) 

# get monthly sums of proportion of days with wind from each direction
prop.SE <- tapply(dat$SE, list(yr,m), f)
prop.NW <- tapply(dat$NW, list(yr,m), f)

# Danielson et al. 2012 GRL recommend seasons of Oct-Apr and May-Sept

MaySepNW <- prop.NW[,5:9] 
MaySepSE <- prop.SE[,5:9]

# now get summer means
sumNW <- rowMeans(MaySepNW)
plot(names(sumNW), sumNW, type="o")
sumSE <- rowMeans(MaySepSE)
plot(names(sumSE), sumSE, type="o")

# and winter...
OctAprNW <- prop.NW[2:nrow(prop.NW), 1:4]
#add in Oct Nov Dec from previous year
OctAprNW <- cbind(OctAprNW, prop.NW[1:(nrow(prop.NW)-1),10:12]) 
colMeans(OctAprNW, na.rm=T) # check how the months compare - pretty similar

# get means
winNW <- rowMeans(OctAprNW)
plot(names(winNW), winNW, type="o") # more of a coherent trend than the summer TS!

OctAprSE <- prop.SE[2:nrow(prop.SE), 1:4]
#add in Oct Nov Dec from previous year
OctAprSE <- cbind(OctAprSE, prop.SE[1:(nrow(prop.SE)-1),10:12]) 
colMeans(OctAprSE, na.rm=T) # check how the months compare - pretty similar

# get means
winSE <- rowMeans(OctAprSE)
plot(names(winSE), winSE, type="o") 

# and save as a .csv file
out <- cbind(winSE, winNW, sumSE[2:length(sumSE)], sumNW[2:length(sumNW)])
colnames(out)[3:4] <- c("sumSE", "sumNW")
write.csv(out, "ebs.winds.csv")

# out of curiosity...
cor(out, use="p") # quite independent of each other!
