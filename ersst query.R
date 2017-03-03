setwd("/Users/MikeLitzow/Documents/R/FATE")
library(ncdf4)
library(maps)
library(mapdata)
library(chron)
library(fields)
require(zoo)
require(ape)

# dowload the latest version of ERSSTv4 
# download.file(url="ftp://ftp.cdc.noaa.gov/Datasets/noaa.ersst/sst.mnmean.v4.nc", destfile= "/Users/MikeLitzow/Documents/R/NSF-GOA/sst.mnmean.v4.nc")

# now load the file
nc <- nc_open("/Users/MikeLitzow/Documents/R/NSF-GOA/sst.mnmean.v4.nc")

# get info
nc

# view dates (middle of month):
ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)

# assign dates
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# we will start with Jan 1950 and end with most recent month
d[c(1153,length(d))] # check - those are correct
d <- d[1153:length(d)] # restrict to the desired dates

# save months and years for use later on
m <- months(d)
yrs <- years(d)

# and set a couple functions for standardizing below
f1 <- function(x) tapply(x, m, mean)
f2 <- function(x) tapply(x, m, sd)

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

##############################################################################
# Western GOA
# 54-62N and 198-216E
x <- ncvar_get(nc, "lon", start=100, count=10)
y <- ncvar_get(nc, "lat", start=14, count=5)
x;y # check

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(100,14,1153), count=c(10,5,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out Bristol Bay and far offshore areas
blank <- c("N56E198", "N58E198", "N60E198","N62E198",  "N56E200", "N58E200", "N58E202", "N54E216",
           "N56E216", "N54E214", "N56E214", "N54E212", "N54E210", "N54E208")
SST[,blank] <- NA

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(40,66))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# check with some plots!
par(las=1)
wgoa.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(wgoa.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("WGOA SST & 13-mo rolling mean")

# calculate sd and ar(1) on rolling windows=25% of total length
w <- round(0.25*nrow(SST))

# create objects to store rolling mean, sd, and ar(1) values
wgoa.rmean <- wgoa.rsd <- wgoa.rar <- rep(NA, length.out=nrow(SST))
SST.mean <- rowMeans(SST.anom,na.rm=T)
for(i in 1:(nrow(SST)-w)){
 # i <- 1
  wgoa.rmean[i+round(w/2)] <- mean(SST.mean[i:(i+(w-1))])
  wgoa.rsd[i+round(w/2)] <- sd(SST.mean[i:(i+(w-1))])
  wgoa.rar[i+round(w/2)] <- ar(SST.mean[i:(i+(w-1))], aic=F, order.max=1)$ar
}

# scale them all for comparison!
wgoa.rmean <- scale(wgoa.rmean)
wgoa.rsd <- scale(wgoa.rsd)
wgoa.rar <- scale(wgoa.rar)

# and plot
rmean.plot <- ts(wgoa.rmean, start=c(1950,1), frequency=12)
rsd.plot <- ts(wgoa.rsd, start=c(1950,1), frequency=12)
rar.plot <- ts(wgoa.rar, start=c(1950,1), frequency=12)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
      xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("bottomright", c("Mean", "Temp. SD", "Temp. AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("WGOA SST")

# now plot spatial metrics

wgoa.sp.sd <- apply(SST.anom, 1, sd, na.rm=T)

# now calculate distance matrix for spatial autocorrelation
# using Moran's I in package ape
dd <- rdist.earth(cbind(lon, lat), miles = F)
dd.inv <- 1/dd # we weight by the inverse of distance
diag(dd.inv) <- 0

wgoa.sp.ar <- NA # object for saving Moran scores

for(i in 1:nrow(SST.anom)){
# i <- 1
mi <- Moran.I(SST.anom[i,], dd.inv, na.rm = T)
wgoa.sp.ar[i] <- mi$observed
}

# scale and plot
mean.plot <- ts(rollmean(scale(SST.mean), 37, fill=NA), start=c(1950,1), frequency=12)
sp.sd.plot <- ts(rollmean(scale(wgoa.sp.sd), 37, fill=NA), start=c(1950,1), frequency=12)
sp.ar.plot <- ts(rollmean(scale(wgoa.sp.ar), 37, fill=NA), start=c(1950,1), frequency=12)


plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean", "Sp. SD", "Sp. cor."), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("WGOA SST - 37 mo. rolling means")

pdf("WGOA SST plots.pdf", 10.5, 8)
par(mfrow=c(2,2), las=1, mar=c(4,4,2,1), oma=c(0,0,2,0), mgp=c(2.5,1,0))

SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x,y,z, col=tim.colors(64), zlim=c(2,17), xlim=c(160,240), ylim=c(40,66), ylab="N latitude"
           , xlab="Longitude", xaxt="n")
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
axis(1, at=seq(160,240,by=20), labels=c("160E", "180", "160W", "140W", "120W"))
mtext("A) Mean SST pattern", adj=0)

wgoa.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(wgoa.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("B) Monthly anomalies & 13-month rolling mean", adj=0)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("bottomright", c("Mean SST", "Temporal SD", "Temporal AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
mtext("C) Temporal EWI (rolling window = 25% time series length)", adj=0)

plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean SST", "Spatial SD", "Spatial correlation"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
mtext("D) Spatial EWI (smoothed with 37-month rolling mean)", adj=0)
mtext("Western Gulf of Alaska", outer=T, cex=2)
dev.off()

# rename and save the various TS
SST.dat <- data.frame(wgoa.mean=rowMeans(SST.anom, na.rm=T), wgoa.roll.mean.25percent=wgoa.rmean, 
                      wgoa.roll.sd=wgoa.rsd, wgoa.roll.ar=wgoa.rar, 
                      wgoa.spatial.sd=wgoa.sp.sd, wgoa.spatial.cor=wgoa.sp.ar)
head(SST.dat)

####
# now calculate seasonal means
# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
spr <- c("Apr", "May", "Jun")

# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]
plot(names(SST.win), SST.win, type="b")

# finally, add a NA for 1950 to align with spring!
SST.win <- c(NA, SST.win)
names(SST.win)[1] <- 1950

# calculate spring mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% spr # logical vector that's true for spring months
SST.mean <- SST.mean[use] # select spring means only
spr.yrs <- yrs[use] # restrict years vector to spring months
SST.spr <- tapply(SST.mean, spr.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
spr.count <- tapply(SST.mean, spr.yrs, f3)==3 # logical vector to identify complete winters!
# there's an oddity with spring months - we can have a year for which spr.count==NA
# fix that here
spr.count[is.na(spr.count)] <- F

SST.spr <- SST.spr[spr.count] # restrict to complete cases
plot(names(SST.spr), SST.spr, type="b")

# now ensure that win and spr time series are the same length
if(length(SST.win)==length(SST.spr)) SST.spr <- SST.spr else SST.spr <- c(SST.spr,NA)

seas.SST <- data.frame(wgoa.NDJFM=SST.win, wgoa.AMJ=SST.spr)
##########################################################################################
##########################################################################################
# now the Eastern GOA
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

# 52-60N and 218-232E
x <- ncvar_get(nc, "lon", start=110, count=8)
y <- ncvar_get(nc, "lat", start=15, count=5)
x;y # check

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(110,15,1153), count=c(8,5,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out the far offshore areas
blank <- c("N52E218", "N54E218", "N56E218", "N52E220", "N54E220", "N52E222")
SST[,blank] <- NA


# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(40,66))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# check with some plots!
par(las=1)
egoa.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(egoa.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("EGOA SST & 13-mo rolling mean")

# calculate sd and ar(1) on rolling windows=25% of total length
w <- round(0.25*nrow(SST))

# create objects to store rolling mean, sd, and ar(1) values
egoa.rmean <- egoa.rsd <- egoa.rar <- rep(NA, length.out=nrow(SST))
SST.mean <- rowMeans(SST.anom,na.rm=T)
for(i in 1:(nrow(SST)-w)){
  # i <- 1
  egoa.rmean[i+round(w/2)] <- mean(SST.mean[i:(i+(w-1))])
  egoa.rsd[i+round(w/2)] <- sd(SST.mean[i:(i+(w-1))])
  egoa.rar[i+round(w/2)] <- ar(SST.mean[i:(i+(w-1))], aic=F, order.max=1)$ar
}

# scale them all for comparison!
egoa.rmean <- scale(egoa.rmean)
egoa.rsd <- scale(egoa.rsd)
egoa.rar <- scale(egoa.rar)

# and plot
rmean.plot <- ts(egoa.rmean, start=c(1950,1), frequency=12)
rsd.plot <- ts(egoa.rsd, start=c(1950,1), frequency=12)
rar.plot <- ts(egoa.rar, start=c(1950,1), frequency=12)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("bottomright", c("Mean", "Temp. SD", "Temp. AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("EGOA SST")

# now plot spatial metrics


egoa.sp.sd <- apply(SST.anom, 1, sd, na.rm=T)

# now calculate distance matrix for spatial autocorrelation
# using Moran's I in package ape
dd <- rdist.earth(cbind(lon, lat), miles = F)
dd.inv <- 1/dd # we weight by the inverse of distance
diag(dd.inv) <- 0

egoa.sp.ar <- NA # object for saving Moran scores

for(i in 1:nrow(SST.anom)){
  # i <- 1
  mi <- Moran.I(SST.anom[i,], dd.inv, na.rm = T)
  egoa.sp.ar[i] <- mi$observed
}

# scale and plot
mean.plot <- ts(rollmean(scale(SST.mean), 37, fill=NA), start=c(1950,1), frequency=12)
sp.sd.plot <- ts(rollmean(scale(egoa.sp.sd), 37, fill=NA), start=c(1950,1), frequency=12)
sp.ar.plot <- ts(rollmean(scale(egoa.sp.ar), 37, fill=NA), start=c(1950,1), frequency=12)


plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean", "Sp. SD", "Sp. cor."), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("EGOA SST - 37 mo. rolling means")

# plot as a pdf to summarize info across areas
pdf("EGOA SST plots.pdf", 10.5, 8)
par(mfrow=c(2,2), las=1, mar=c(4,4,2,1), oma=c(0,0,2,0), mgp=c(2.5,1,0))

SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x,y,z, col=tim.colors(64), zlim=c(2,17), xlim=c(160,240), ylim=c(40,66), ylab="N latitude"
           , xlab="Longitude", xaxt="n")
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
axis(1, at=seq(160,240,by=20), labels=c("160E", "180", "160W", "140W", "120W"))
mtext("A) Mean SST pattern", adj=0)

egoa.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(egoa.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("B) Monthly anomalies & 13-month rolling mean", adj=0)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("bottomright", c("Mean SST", "Temporal SD", "Temporal AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
mtext("C) Temporal EWI (rolling window = 25% time series length)", adj=0)

plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean SST", "Spatial SD", "Spatial correlation"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
mtext("D) Spatial EWI (smoothed with 37-month rolling mean)", adj=0)
mtext("Eastern Gulf of Alaska", outer=T, cex=2)
dev.off()


# rename and save the various TS
egoa.mean <- rowMeans(SST.anom, na.rm=T)
egoa.roll.mean.25percent <- egoa.rmean
egoa.roll.sd <- egoa.rsd
egoa.roll.ar <- egoa.rar
egoa.spatial.sd <- egoa.sp.sd
egoa.spatial.cor <- egoa.sp.ar

SST.dat <- cbind(SST.dat, egoa.mean, egoa.roll.mean.25percent, egoa.roll.sd,
                  egoa.roll.ar, egoa.spatial.sd, egoa.spatial.cor)
head(SST.dat)

####
# now calculate seasonal means
# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
spr <- c("Apr", "May", "Jun")

# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]
plot(names(SST.win), SST.win, type="b")

# finally, add a NA for 1950 to align with spring!
SST.win <- c(NA, SST.win)
names(SST.win)[1] <- 1950

# calculate spring mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% spr # logical vector that's true for spring months
SST.mean <- SST.mean[use] # select spring means only
spr.yrs <- yrs[use] # restrict years vector to spring months
SST.spr <- tapply(SST.mean, spr.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
spr.count <- tapply(SST.mean, spr.yrs, f3)==3 # logical vector to identify complete winters!
# there's an oddity with spring months - we can have a year for which spr.count==NA
# fix that here
spr.count[is.na(spr.count)] <- F

SST.spr <- SST.spr[spr.count] # restrict to complete cases
plot(names(SST.spr), SST.spr, type="b")

# now ensure that win and spr time series are the same length
if(length(SST.win)==length(SST.spr)) SST.spr <- SST.spr else SST.spr <- c(SST.spr,NA)

egoa.NDJFM <- SST.win
egoa.AMJ <- SST.spr

seas.SST <- cbind(seas.SST, egoa.NDJFM, egoa.AMJ)

##############################################################
#######################################################################
##############################################################################
# EBS
# 54-62 deg. N, 186-202 deg. E
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=9)
y <- ncvar_get(nc, "lat", start=14, count=5)
x; y


# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,14,1153), count=c(9,5,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out as needed
blank <- c("N54E186", "N54E188", "N54E190",  "N54E196", "N54E198", "N54E200", "N54E202",
           "N56E186", "N56E188", "N58E186",  "N56E200", "N56E202")
SST[,blank] <- NA

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(40,66))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# check with some plots!
par(las=1)
ebs.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(ebs.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("EBS SST & 13-mo rolling mean")

# calculate sd and ar(1) on rolling windows=25% of total length
w <- round(0.25*nrow(SST))

# create objects to store rolling mean, sd, and ar(1) values
ebs.rmean <- ebs.rsd <- ebs.rar <- rep(NA, length.out=nrow(SST))
SST.mean <- rowMeans(SST.anom,na.rm=T)
for(i in 1:(nrow(SST)-w)){
  # i <- 1
  ebs.rmean[i+round(w/2)] <- mean(SST.mean[i:(i+(w-1))])
  ebs.rsd[i+round(w/2)] <- sd(SST.mean[i:(i+(w-1))])
  ebs.rar[i+round(w/2)] <- ar(SST.mean[i:(i+(w-1))], aic=F, order.max=1)$ar
}

# scale them all for comparison!
ebs.rmean <- scale(ebs.rmean)
ebs.rsd <- scale(ebs.rsd)
ebs.rar <- scale(ebs.rar)

# and plot
rmean.plot <- ts(ebs.rmean, start=c(1950,1), frequency=12)
rsd.plot <- ts(ebs.rsd, start=c(1950,1), frequency=12)
rar.plot <- ts(ebs.rar, start=c(1950,1), frequency=12)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("bottomright", c("Mean", "Temp. SD", "Temp. AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("EBS SST")

# now plot spatial metrics


ebs.sp.sd <- apply(SST.anom, 1, sd, na.rm=T)

# now calculate distance matrix for spatial autocorrelation
# using Moran's I in package ape
dd <- rdist.earth(cbind(lon, lat), miles = F)
dd.inv <- 1/dd # we weight by the inverse of distance
diag(dd.inv) <- 0

ebs.sp.ar <- NA # object for saving Moran scores

for(i in 1:nrow(SST.anom)){
  # i <- 1
  mi <- Moran.I(SST.anom[i,], dd.inv, na.rm = T)
  ebs.sp.ar[i] <- mi$observed
}

# scale and plot
mean.plot <- ts(rollmean(scale(SST.mean), 37, fill=NA), start=c(1950,1), frequency=12)
sp.sd.plot <- ts(rollmean(scale(ebs.sp.sd), 37, fill=NA), start=c(1950,1), frequency=12)
sp.ar.plot <- ts(rollmean(scale(ebs.sp.ar), 37, fill=NA), start=c(1950,1), frequency=12)


plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean", "Sp. SD", "Sp. cor."), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("EBS SST - 37 mo. rolling means")

  pdf("EBS SST plots.pdf", 10.5, 8)
  par(mfrow=c(2,2), las=1, mar=c(4,4,2,1), oma=c(0,0,2,0), mgp=c(2.5,1,0))
  
  SST.mean <- colMeans(SST)
  z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image.plot(x,y,z, col=tim.colors(64), zlim=c(2,17), xlim=c(160,240), ylim=c(40,66), ylab="N latitude"
             , xlab="Longitude", xaxt="n")
  map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
  axis(1, at=seq(160,240,by=20), labels=c("160E", "180", "160W", "140W", "120W"))
  mtext("A) Mean SST pattern", adj=0)
  
  ebs.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
  plot(ebs.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
  sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
  lines(sm, col="red", lwd=1.5)
  abline(h=0)
  mtext("B) Monthly anomalies & 13-month rolling mean", adj=0)
  
  plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
       xlab="", ylab="Anomaly")
  lines(rsd.plot, lwd=1.5, col="#009E73")
  lines(rar.plot, lwd=1.5, col="#0072B2")
  abline(h=0, col="grey")
  legend("bottomright", c("Mean SST", "Temporal SD", "Temporal AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
  mtext("C) Temporal EWI (rolling window = 25% time series length)", adj=0)
  
  plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
       xlab="", ylab="Anomaly")
  lines(sp.sd.plot, lwd=1.5, col="#009E73")
  lines(sp.ar.plot, lwd=1.5, col="#0072B2")
  abline(h=0, col="grey")
  legend("topleft", c("Mean SST", "Spatial SD", "Spatial correlation"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
  mtext("D) Spatial EWI (smoothed with 37-month rolling mean)", adj=0)
  mtext("Eastern Bering Sea", outer=T, cex=2)
  dev.off()
  
# now rename and save the TS
ebs.mean <- rowMeans(SST.anom, na.rm=T)
ebs.roll.mean.25percent <- ebs.rmean
ebs.roll.sd <- ebs.rsd
ebs.roll.ar <- ebs.rar
ebs.spatial.sd <- ebs.sp.sd
ebs.spatial.cor <- ebs.sp.ar

SST.dat <- cbind(SST.dat, ebs.mean, ebs.roll.mean.25percent, ebs.roll.sd,
                 ebs.roll.ar, ebs.spatial.sd, ebs.spatial.cor)

head(SST.dat)

####
# now calculate seasonal means
# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
spr <- c("Apr", "May", "Jun")

# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]
plot(names(SST.win), SST.win, type="b")

# finally, add a NA for 1950 to align with spring!
SST.win <- c(NA, SST.win)
names(SST.win)[1] <- 1950

# calculate spring mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% spr # logical vector that's true for spring months
SST.mean <- SST.mean[use] # select spring means only
spr.yrs <- yrs[use] # restrict years vector to spring months
SST.spr <- tapply(SST.mean, spr.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
spr.count <- tapply(SST.mean, spr.yrs, f3)==3 # logical vector to identify complete winters!
# there's an oddity with spring months - we can have a year for which spr.count==NA
# fix that here
spr.count[is.na(spr.count)] <- F

SST.spr <- SST.spr[spr.count] # restrict to complete cases
plot(names(SST.spr), SST.spr, type="b")

# now ensure that win and spr time series are the same length
if(length(SST.win)==length(SST.spr)) SST.spr <- SST.spr else SST.spr <- c(SST.spr,NA)

ebs.NDJFM <- SST.win
ebs.AMJ <- SST.spr

seas.SST <- cbind(seas.SST, ebs.NDJFM, ebs.AMJ)

#########################
##############################################################
####################################################################
##############################################################################
# NCC

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

# 42-50 deg. N, 234-242 deg. E
x <- ncvar_get(nc, "lon", start=117, count=3)
y <- ncvar_get(nc, "lat", start=20, count=5)
x; y 


# get required sst data
SST <- ncvar_get(nc, "sst", start=c(117,20,1153), count=c(3,5,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out non-coastal waters
#blank <- c("N32E236", "N34E236", "N36E236", "N32E238","N34E238", "N32E240", "N32E244")

#SST[,blank] <- NA

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(180,246), ylim=c(20,56))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# check with some plots!
par(las=1)
ncc.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(ncc.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("NCC SST & 13-mo rolling mean")

# calculate sd and ar(1) on rolling windows=25% of total length
w <- round(0.25*nrow(SST))

# create objects to store rolling mean, sd, and ar(1) values
ncc.rmean <- ncc.rsd <- ncc.rar <- rep(NA, length.out=nrow(SST))
SST.mean <- rowMeans(SST.anom,na.rm=T)
for(i in 1:(nrow(SST)-w)){
  # i <- 1
  ncc.rmean[i+round(w/2)] <- mean(SST.mean[i:(i+(w-1))])
  ncc.rsd[i+round(w/2)] <- sd(SST.mean[i:(i+(w-1))])
  ncc.rar[i+round(w/2)] <- ar(SST.mean[i:(i+(w-1))], aic=F, order.max=1)$ar
}

# scale them all for comparison!
ncc.rmean <- scale(ncc.rmean)
ncc.rsd <- scale(ncc.rsd)
ncc.rar <- scale(ncc.rar)

# and plot
rmean.plot <- ts(ncc.rmean, start=c(1950,1), frequency=12)
rsd.plot <- ts(ncc.rsd, start=c(1950,1), frequency=12)
rar.plot <- ts(ncc.rar, start=c(1950,1), frequency=12)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("bottomright", c("Mean", "Temp. SD", "Temp. AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("NCC SST")

# now plot spatial metrics
ncc.sp.sd <- apply(SST.anom, 1, sd, na.rm=T)

# now calculate distance matrix for spatial autocorrelation
# using Moran's I in package ape
dd <- rdist.earth(cbind(lon, lat), miles = F)
dd.inv <- 1/dd # we weight by the inverse of distance
diag(dd.inv) <- 0

ncc.sp.ar <- NA # object for saving Moran scores

for(i in 1:nrow(SST.anom)){
  # i <- 1
  mi <- Moran.I(SST.anom[i,], dd.inv, na.rm = T)
  ncc.sp.ar[i] <- mi$observed
}

# scale and plot
mean.plot <- ts(rollmean(scale(SST.mean), 37, fill=NA), start=c(1950,1), frequency=12)
sp.sd.plot <- ts(rollmean(scale(ncc.sp.sd), 37, fill=NA), start=c(1950,1), frequency=12)
sp.ar.plot <- ts(rollmean(scale(ncc.sp.ar), 37, fill=NA), start=c(1950,1), frequency=12)


plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean", "Sp. SD", "Sp. cor."), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("NCC SST - 37 mo. rolling means")

pdf("NCC SST plots.pdf", 10.5, 8)
par(mfrow=c(2,2), las=1, mar=c(4,4,2,1), oma=c(0,0,2,0), mgp=c(2.5,1,0))

SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x,y,z, col=tim.colors(64), zlim=c(2,17), xlim=c(180,246), ylim=c(20,56), ylab="N latitude"
           , xlab="W longitude", xaxt="n")
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
axis(1, at=seq(180,240,by=10), labels=seq(180,120,by=-10))
mtext("A) Mean SST pattern", adj=0)

ncc.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(ncc.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("B) Monthly anomalies & 13-month rolling mean", adj=0)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("bottomright", c("Mean SST", "Temporal SD", "Temporal AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
mtext("C) Temporal EWI (rolling window = 25% time series length)", adj=0)

plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean SST", "Spatial SD", "Spatial correlation"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
mtext("D) Spatial EWI (smoothed with 37-month rolling mean)", adj=0)
mtext("Northern California Current", outer=T, cex=2)
dev.off()


# now rename and save the TS
ncc.mean <- rowMeans(SST.anom, na.rm=T)
ncc.roll.mean.25percent <- ncc.rmean
ncc.roll.sd <- ncc.rsd
ncc.roll.ar <- ncc.rar
ncc.spatial.sd <- ncc.sp.sd
ncc.spatial.cor <- ncc.sp.ar

SST.dat <- cbind(SST.dat, ncc.mean, ncc.roll.mean.25percent, ncc.roll.sd,
                 ncc.roll.ar, ncc.spatial.sd, ncc.spatial.cor)

colnames(SST.dat)

####
# now calculate seasonal means
# define seasons for CCE - this is based on Black et alia's (2015 GCB) definition 
# of winter and summer seasons for upwelling
win <- c("Jan", "Feb","Mar")
sum <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
# get a vector of years for each observation
yrs <- years(names(SST.mean))
SST.win <- tapply(SST.mean, yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, yrs, f3)==3 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]
plot(names(SST.win), SST.win, type="b")

# calculate summer mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% sum # logical vector that's true for summer months
SST.mean <- SST.mean[use] # select summer means only
sum.yrs <- years(names(SST.mean)) # restrict years vector to summer months
SST.sum <- tapply(SST.mean, sum.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
sum.count <- tapply(SST.mean, sum.yrs, f3)==7 # logical vector to identify complete winters!

SST.sum <- SST.sum[sum.count] # restrict to complete cases
plot(names(SST.spr), SST.spr, type="b")

# now ensure that win and sum time series are the same length
if(length(SST.win)==length(SST.sum)) SST.sum <- SST.sum else SST.sum <- c(SST.sum,NA)

ncc.JFM <- SST.win
ncc.AMJJASO <- SST.sum

seas.SST <- cbind(seas.SST, ncc.JFM, ncc.AMJJASO)

#######################################
###############################
###################
# SCC

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

# 32-40 deg. N, 234-242 deg. E
x <- ncvar_get(nc, "lon", start=118, count=5)
y <- ncvar_get(nc, "lat", start=25, count=5)
x; y 

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(118,25,1153), count=c(5,5,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out non-coastal waters
blank <- c("N32E234", "N34E234", "N32E236", "N34E236", "N36E234", "N32E238")

SST[,blank] <- NA

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(180,246), ylim=c(20,56))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# check with some plots!
par(las=1)
scc.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(scc.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("SCC SST & 13-mo rolling mean")

# calculate sd and ar(1) on rolling windows=25% of total length
w <- round(0.25*nrow(SST))

# create objects to store rolling mean, sd, and ar(1) values
scc.rmean <- scc.rsd <- scc.rar <- rep(NA, length.out=nrow(SST))
SST.mean <- rowMeans(SST.anom,na.rm=T)
for(i in 1:(nrow(SST)-w)){
  # i <- 1
  scc.rmean[i+round(w/2)] <- mean(SST.mean[i:(i+(w-1))])
  scc.rsd[i+round(w/2)] <- sd(SST.mean[i:(i+(w-1))])
  scc.rar[i+round(w/2)] <- ar(SST.mean[i:(i+(w-1))], aic=F, order.max=1)$ar
}

# scale them all for comparison!
scc.rmean <- scale(scc.rmean)
scc.rsd <- scale(scc.rsd)
scc.rar <- scale(scc.rar)

# and plot
rmean.plot <- ts(scc.rmean, start=c(1950,1), frequency=12)
rsd.plot <- ts(scc.rsd, start=c(1950,1), frequency=12)
rar.plot <- ts(scc.rar, start=c(1950,1), frequency=12)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("bottomright", c("Mean", "Temp. SD", "Temp. AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("SCC SST")

# now plot spatial metrics
scc.sp.sd <- apply(SST.anom, 1, sd, na.rm=T)

# now calculate distance matrix for spatial autocorrelation
# using Moran's I in package ape
dd <- rdist.earth(cbind(lon, lat), miles = F)
dd.inv <- 1/dd # we weight by the inverse of distance
diag(dd.inv) <- 0

scc.sp.ar <- NA # object for saving Moran scores

for(i in 1:nrow(SST.anom)){
  # i <- 1
  mi <- Moran.I(SST.anom[i,], dd.inv, na.rm = T)
  scc.sp.ar[i] <- mi$observed
}

# scale and plot
mean.plot <- ts(rollmean(scale(SST.mean), 37, fill=NA), start=c(1950,1), frequency=12)
sp.sd.plot <- ts(rollmean(scale(scc.sp.sd), 37, fill=NA), start=c(1950,1), frequency=12)
sp.ar.plot <- ts(rollmean(scale(scc.sp.ar), 37, fill=NA), start=c(1950,1), frequency=12)


plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean", "Sp. SD", "Sp. cor."), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n", cex=0.8)
mtext("SCC SST - 37 mo. rolling means")

pdf("SCC SST plots.pdf", 10.5, 8)
par(mfrow=c(2,2), las=1, mar=c(4,4,2,1), oma=c(0,0,2,0), mgp=c(2.5,1,0))

SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x,y,z, col=tim.colors(64), zlim=c(2,17), xlim=c(180,246), ylim=c(20,56), ylab="N latitude"
           , xlab="W longitude", xaxt="n")
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
axis(1, at=seq(180,240,by=10), labels=seq(180,120,by=-10))
mtext("A) Mean SST pattern", adj=0)

scc.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(scc.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("B) Monthly anomalies & 13-month rolling mean", adj=0)

plot(rmean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(rmean.plot, rsd.plot, rar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(rsd.plot, lwd=1.5, col="#009E73")
lines(rar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean SST", "Temporal SD", "Temporal AR(1)"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
mtext("C) Temporal EWI (rolling window = 25% time series length)", adj=0)

plot(mean.plot, type="l", lwd=1.5, col="#E69F00", ylim=range(mean.plot, sp.sd.plot, sp.ar.plot, na.rm=T),
     xlab="", ylab="Anomaly")
lines(sp.sd.plot, lwd=1.5, col="#009E73")
lines(sp.ar.plot, lwd=1.5, col="#0072B2")
abline(h=0, col="grey")
legend("topleft", c("Mean SST", "Spatial SD", "Spatial correlation"), text.col=c("#E69F00", "#009E73", "#0072B2"), bty="n")
mtext("D) Spatial EWI (smoothed with 37-month rolling mean)", adj=0)
mtext("Southern California Current", outer=T, cex=2)
dev.off()


# now rename and save the TS
scc.mean <- rowMeans(SST.anom, na.rm=T)
scc.roll.mean.25percent <- scc.rmean
scc.roll.sd <- scc.rsd
scc.roll.ar <- scc.rar
scc.spatial.sd <- scc.sp.sd
scc.spatial.cor <- scc.sp.ar

SST.dat <- cbind(SST.dat, scc.mean, scc.roll.mean.25percent, scc.roll.sd,
                 scc.roll.ar, scc.spatial.sd, scc.spatial.cor)

colnames(SST.dat)
write.csv(SST.dat, "sst data by system all months.csv")

####
# now calculate seasonal means
# define seasons for CCE - this is based on Black et alia's (2015 GCB) definition 
# of winter and summer seasons for upwelling
win <- c("Jan", "Feb","Mar")
sum <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
# get a vector of years for each observation
yrs <- years(names(SST.mean))
SST.win <- tapply(SST.mean, yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, yrs, f3)==3 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]
plot(names(SST.win), SST.win, type="b")

# calculate summer mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% sum # logical vector that's true for summer months
SST.mean <- SST.mean[use] # select summer means only
sum.yrs <- years(names(SST.mean)) # restrict years vector to summer months
SST.sum <- tapply(SST.mean, sum.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
sum.count <- tapply(SST.mean, sum.yrs, f3)==7 # logical vector to identify complete winters!

SST.sum <- SST.sum[sum.count] # restrict to complete cases
plot(names(SST.spr), SST.spr, type="b")

# now ensure that win and sum time series are the same length
if(length(SST.win)==length(SST.sum)) SST.sum <- SST.sum else SST.sum <- c(SST.sum,NA)

scc.JFM <- SST.win
scc.AMJJASO <- SST.sum

seas.SST <- cbind(seas.SST, scc.JFM, scc.AMJJASO)

write.csv(seas.SST, "seasonal SST by system.csv")
