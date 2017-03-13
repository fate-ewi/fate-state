setwd("/Users/MikeLitzow/Documents/R/NSF-GOA")
require(mgcv)
require(zoo)
require(lubridate)
require(dplyr)

# gak1 data can be downloaded at:
# http://www.ims.uaf.edu/gak1/data/TimeSeries/gak1.dat

data <- read.csv("gak1 2-17-17.csv")
head(data)

# appears to be a bad date in the 2013 data - asked Seth Danielson about this (2.17.17)
# dropping for now
## NOW THAT CRUISE DATE HAS BEEN CORRECTED!
# data[data$cruise == "LO128",]

# data <- filter(data, cruise!="LO128")


setwd("/Users/MikeLitzow/Documents/R/FATE/fate-ewi")

data$year <- floor(data$dec.year)

data$day <- yday(date_decimal(data$dec.year))


head(data)

sal.set <- as.data.frame(tapply(data$sigma.t, list(data$dec.year, data$depth), mean))
colnames(sal.set) <- c("d0", "d10", "d20", "d30", "d50", "d75", "d100", "d150", "d200", "d250")
cor(sigma.set, use="pa")
cor(sal.set, use="pa")

sal.set$year <- floor(as.numeric(rownames(sal.set))) 
sal.set$day <- yday(date_decimal(as.numeric(rownames(sal.set))))

# now get FMA salinity
salFMA <- filter(sal.set, day>31 & day <=120)

sal10mu <- tapply(salFMA$d10, salFMA$year, mean, na.rm=T)
sal20mu <- tapply(salFMA$d20, salFMA$year, mean, na.rm=T)
sal30mu <- tapply(salFMA$d30, salFMA$year, mean, na.rm=T)

cor(cbind(sal10mu, sal20mu, sal30mu)) # all the same!

 !is.na(salFMA$d30) # no missing values for 30 m salinity!
sal30day <- tapply(salFMA$day, salFMA$year, mean, na.rm=T)
sal30n <- tapply(salFMA$d30, salFMA$year, function(x) sum(!is.na(x)))

plot(names(sal30mu), sal30mu, type="o") # getting fresher!

mod <- gam(sal30mu ~ s(sal30day, k=4))
summary(mod) # best fit is linear!

adj.sal30 <- residuals(mod, type="response") + summary(mod)$p.coeff

# and check
plot(names(sal30mu), sal30mu, type="o", col="blue")
lines(names(adj.sal30), adj.sal30, type="o", col= "dark green")
# only makes a small difference in a couple years

sal.out <- data.frame(raw.sal30=sal30mu, adj.sal30=adj.sal30, mean.day=sal30day)

write.csv(sal.out, "GAK1 30m salinity.csv")

