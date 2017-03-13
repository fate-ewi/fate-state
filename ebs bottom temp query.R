setwd("/Users/MikeLitzow/Documents/R/seagrant #2")
require(mgcv)
require(zoo)
require(mgcv)
require(mice)
require(Hmisc)

# load haul data for 1982-2015
BShauls <- read.csv("hauls2015.csv")

# now load 2013-2016 data and pull out 2016 haul data
setwd("/Users/MikeLitzow/Documents/R/FATE")
xtra <- read.csv("ebs2013_2016.csv")

# stations sampled every year
st <- read.csv("Bering Sea trawl stations sampled every year 1982-2016.csv")
# copying them in here for reference....
st[,2]

# [1] A-02   A-03   A-04   A-05   B-02   B-03   B-04   B-05   B-06   B-07   C-01   C-02   C-03   C-04
# [15] C-05   C-06   C-07   C-08   C-09   C-18   D-01   D-02   D-04   D-05   D-06   D-07   D-08   D-09
# [29] D-10   D-18   E-01   E-02   E-03   E-05   E-06   E-07   E-08   E-09   E-10   E-12   E-18   E-19
# [43] E-21   E-22   F-01   F-02   F-03   F-04   F-05   F-06   F-07   F-08   F-09   F-10   F-11   F-12
# [57] F-13   F-14   F-18   F-19   F-21   F-22   F-24   F-25   G-01   G-02   G-03   G-04   G-05   G-06
# [71] G-07   G-08   G-09   G-10   G-11   G-12   G-13   G-14   G-15   G-18   G-19   G-20   G-21   G-22
# [85] G-23   G-25   G-26   GF1918 GF2019 GF2120 GF2221 H-01   H-02   H-03   H-04   H-05   H-06   H-07
# [99] H-08   H-09   H-10   H-11   H-12   H-13   H-14   H-15   H-16   H-18   H-19   H-20   H-21   H-22
# [113] H-23   H-24   H-25   H-26   HG1918 HG2019 I-01   I-02   I-03   I-04   I-05   I-06   I-07   I-08
# [127] I-09   I-10   I-11   I-12   I-14   I-15   I-16   I-18   I-19   I-20   I-21   I-22   I-23   I-24
# [141] I-26   IH1918 IH2019 IH2120 IH2221 J-01   J-02   J-03   J-04   J-05   J-06   J-07   J-08   J-09
# [155] J-10   J-11   J-13   J-14   J-18   J-19   J-20   J-21   J-22   J-23   J-24   J-25   J-26   JI1918
# [169] JI2019 JI2120 JI2221 K-01   K-02   K-03   K-04   K-05   K-06   K-07   K-08   K-09   K-10   K-11
# [183] K-12   K-14   K-18   K-19   K-20   K-21   K-22   K-23   K-24   K-25   K-26   K-27   L-01   L-02
# [197] L-03   L-04   L-05   L-06   L-07   L-08   L-09   L-18   L-19   L-20   L-21   L-22   L-23   L-24
# [211] L-25   L-26   L-27   L-28   L-29   M-02   M-03   M-04   M-05   M-06   M-07   M-19   M-20   M-21
# [225] M-22   M-23   M-24   M-25   M-26   M-27   M-28   M-29   N-02   N-03   N-04   N-06   N-07   N-21
# [239] N-22   N-23   N-24   N-25   N-26   N-27   N-28   N-29   O-02   O-03   O-04   O-21   O-22   O-23
# [253] O-24   O-25   O-26   O-27   O-28   O-29   O-30   O-31   P-21   P-22   P-23   P-24   P-26   P-27
# [267] P-28   P-29   P-30   P-31   Q-21   Q-22   Q-23   Q-26   Q-27   Q-28   Q-29   Q-30   R-23   R-25
# [281] R-26   R-27   R-28   R-29   R-30

head(xtra)

# get average sampling date and temperature for each
# first, need to convert DATETIME to day of year
xtra$DATETIME <- as.character(xtra$DATETIME)
xt <- strptime(xtra$DATETIME, "%m/%d/%y %H:%M")#specifies abbreviated month and year without century
xtra$JULIAN <- xt$yday

# get date and bottom temp for each station in the new years...
new.date <- tapply(xtra$JULIAN, list(xtra$YEAR, xtra$STATION), mean)
new.temp <- tapply(xtra$BOT_TEMP, list(xtra$YEAR, xtra$STATION), mean)

# restrict to the stations that were sampled every year
keep <- colnames(new.date) %in% st[,2]
new.date <- new.date[,keep]
new.temp <- new.temp[,keep]

# now get the same info for the older hauls
str(BShauls)

# get Julian date
xt <- strptime(as.character(BShauls$START_TIME), "%d-%b-%y")#specifies abbreviated month and year without century
BShauls$JULIAN <- xt$yday

old.date <- tapply(BShauls$JULIAN, list(BShauls$YEAR, BShauls$STATIONID), mean)
old.temp <- tapply(BShauls$GEAR_TEMPERATURE, list(BShauls$YEAR, BShauls$STATIONID), mean)

# restrict to the stations that were sampled every year
keep <- colnames(old.date) %in% st[,2]
old.date <- old.date[,keep]
old.temp <- old.temp[,keep]

# check!
identical(colnames(old.date), colnames(new.date))
dim(old.date)
dim(new.date)

# check again - looks good!
identical(old.date[rownames(old.date)==2015,], new.date[rownames(new.date)==2015,]) #T
identical(old.temp[rownames(old.temp)==2015,], new.temp[rownames(new.temp)==2015,]) #T

# now combine the 2016 data with older data
sample.day <- rbind(old.date, new.date[4,])
b.temp <- rbind(old.temp, new.temp[4,])

rownames(sample.day) <- rownames(b.temp) <- 1982:2016

sum(is.na(sample.day)) # no missing values!
sum(is.na(b.temp)) #415 missing values for temperature!

# just to remove any possible effect of mising values, I'll fill in with multiple imputation

#quickly estimate missing temps with mice

# 1. Determine "best" predictor variables for each column:
x <- b.temp
r <- rcorr(x)$r

# Choose 20 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=20))
diag(pred) <- FALSE

Ni <- 20 # number of imputations
imp <- mice(x,Ni, imp="norm", predictor = pred)

#extract the complete imputed data set!
xx <- complete(imp, "long")

means <- matrix(nrow = nrow(x), ncol = ncol(x))
for(i in 1:ncol(x)){
  c <- i+2
  means[,i]<- tapply(xx[,c],xx[,2], mean)
}

colnames(means) <- colnames(x)
rownames(means) <- 1982:2016

# so "means" is the average of the imputations for the missing values

# get mean sampling date and bottom temp for each year
mean.date <- rowMeans(sample.day)
mean.temp <- rowMeans(means)

# quick check to make sure that imputations are realistic!
old.mean <- rowMeans(old.temp, na.rm=T)
new.mean <- rowMeans(new.temp, na.rm=T)
plot(old.mean, mean.temp[names(mean.temp) %in% 1982:2015]) # looks fine!
plot(new.mean, mean.temp[names(mean.temp) %in% 2013:2016]) # looks fine!

# now fit a gam to elucidate the effect of timing of sampling on temperature
# limiting edf (smoothing) with knots=4
mod <- gam(mean.temp ~ s(mean.date, k=4))
summary(mod)
plot(mod, se=T, residuals = T, pch=1, rug=F)

# adjust temp as residuals + intercept
res <- residuals(mod, type="response")

adj.temp <- output$p.coeff + res

# quick comparison of the raw and adjusted temps
plot(1982:2016, mean.temp, type="o", col="blue")
lines(1982:2016, adj.temp, type="o", col= "dark green") # only makes a sizeable difference for 1999

plot(names(adj.temp), adj.temp, type="o")
# just....wow
adj.temp[length(adj.temp)]-max(adj.temp[1:(length(adj.temp)-1)]) # 0.77 degrees above previous max

# and save
ebs.bottom.temp <- data.frame(mean.sample.day = mean.date, raw.temp = mean.temp,
                              adjusted.temp = adj.temp)
write.csv(ebs.bottom.temp, "ebs.bottom.temp.csv")
