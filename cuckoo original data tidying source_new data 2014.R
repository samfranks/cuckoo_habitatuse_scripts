##########################################################
#
#  Cuckoo ORIGINAL data - tidying source file, adding new variables (distance moved between points, bearing, movement groups, and movement types)
#
#  Samantha Franks
#	15 Nov 2013
# 4 Dec 2013
# 23 Dec 2013
# 2 Sep 2014 - re-run with Argos Movebank data (pre-filtered) and new birds
#
#
##########################################################

library(rgdal)
library(maptools)
library(shape)
library(splancs)
library(geosphere)
library(adehabitat) # 
library(sp) # converts R dataframes into spatially explicit dta
library(rgeos)
library(raster)
library(adehabitatHR)
library(chron)
library(plyr)


###-----------------------------------------------###
#         Load cuckoo data & tidying
###-----------------------------------------------###

# cuckoo coordinate locations are in Lat-Long and UTM coordinates and WGS84 datum

setwd("C:/Users/samf/Documents/Git/cuckoos/data")

##############  USES ARGOS MOVEBANK PRE-FILTERED DATA
d <- read.csv("BTO Cuckoo migration study_20140902.csv", header=T)

d2 <- with(d, data.frame(name=individual.local.identifier, id=event.id, timestamp, long=location.long, lat=location.lat, loc.class=argos.lc, error.major=argos.semi.major, error.minor=argos.semi.minor))

d2$long <- as.numeric(as.character(d2$long))

d2$name <- revalue(d2$name, c("lloyd"="Lloyd", "Rob Roy"="Roy", "Jake"="Patch"))


x <- do.call(rbind,(strsplit(as.character(d2$timestamp), " ")))

dat.datetime <- chron(x[,1], x[,2], format=c(dates="d/m/y", times="h:m:s"))

d2 <- data.frame(d2, datetime=dat.datetime, year=years(dat.datetime), month=as.numeric(months(dat.datetime)), day=days(dat.datetime))

d2 <- d2[order(d2$name, d2$datetime),]

# removes locations outside Afro-Palearctic (reference points in North America plus any really weird points in the Atlantic)
d3 <- subset(d2, long > -30)

d4 <- subset(d3, loc.class=="1" | loc.class=="2" | loc.class=="3") # subset by best location classes

d5 <- subset(d4, name!="New Forest 1")
d5 <- droplevels(d5)

###-----------------------------------------------###
#         Add distance, bearing, mgroup and mtype to original data
###-----------------------------------------------###


newdat <- split(d5, list(d5$name)) # split dataset by individual bird

trimdata <- function(dat) { # distance, bearing, etc function
  
  # for each element of data, do the following
  
  # order dat according to location_date - not needed if dataset has already been ordered in Excel by location_date
  
  
  ### --- DISTANCE & BEARING CALCULATION --- ###
  
  # calculate distances and bearings   
  # distbearmvmt <- function(input){ # requires 2 variables called newlongs & newlats or long & lat or centroidlong & centroidlat
  # use newlongs & newlats to calculate distances based on resampled points
  # use long & lat to calculate distances based on original measured location point
  
  input <- dat
  
  dist.cuckoo <- input[,c("long","lat")]
  coordinates(dist.cuckoo) <- c("long", "lat")
  proj4string(dist.cuckoo) <- CRS("+init=epsg:4326")
  
  d1 <- dist.cuckoo@coords
  d2 <- dist.cuckoo@coords[-1,] # create new matrix minus the first row so that d2 starts at the second observation
  last.d2 <- matrix(c(0,0), nrow=1, dimnames=list(1, c("long","lat"))) # create a placeholder last row in the second distance matrix
  d2 <- rbind(d2, c(0,0))
  
  dist <- distCosine(d1,d2)/1000 # distance between points, in km
  bear <- bearing(d1,d2) # bearing between points
  dist <- c(NA,dist[-length(dist)])
  bear <- c(NA,bear[-length(bear)])
  
  distbear <- data.frame(input, distance=dist, bearing=bear)
  
  
  
  ### --- MOVEMENT GROUPS --- ###
  
  mgroup <- rep(NA,nrow(distbear))
  
  for (n in 1:nrow(distbear)){
    if (is.na(distbear$distance[n])) {
      mgroup[n] <- 1
    } else if (distbear$distance[n] <= 30) {
      mgroup[n] <- mgroup[n-1]
    } else {mgroup[n] <- mgroup[n-1] + 1}
  }
  
  ### --- MOVEMENT TYPES --- ###
  
  mtype <- c("C", rep(NA, nrow(distbear)-1))
  
  for (n in 2:(nrow(distbear)-1)){
    if (mgroup[n] != mgroup[n-1] & mgroup[n] != mgroup[n+1]) {
      mtype[n] <- "M"
    } else {
      mtype[n] <- "S"
    }
  }
  
  completedata <- data.frame(distbear,mgroup,mtype)
  return(completedata)
  
}

original <- lapply(newdat, trimdata) # apply function to calculate distance, bearing, mgroup, mtype across all elements of dataset (ie. across each individual bird's data)

fulloriginal <- do.call(rbind,original)

setwd("C:/Users/samf/Documents/Git/cuckoos/data")
write.csv(fulloriginal, "allbirdfinal - Movebank original data new birds 20140903.csv", row.names=F)

setwd("C:/Users/samf/Documents/Git/cuckoos/data/Movebank original data new birds + distance, movement groups, etc")
for (i in 1:length(original)){
  write.csv(original[[i]], paste(original[[i]]$name[1],".csv", sep=""), row.names=FALSE)
}

# # ### Histogram of cuckoo distance movements
# # 
# # hist(fulloriginal$distance)
# # 
# par(mfrow=c(2,2))
# hist(subset(fulloriginal, distance > 20)$distance, main=c("d > 20km"), xlab=c("Distance (km)"))
# hist(subset(fulloriginal, distance > 200 & distance <= 500)$distance, main=c("100km < d < 500km"), xlab=c("Distance (km)"))
# hist(subset(fulloriginal, distance > 500)$distance, main=c("d > 500km"), xlab=c("Distance (km)"))
# hist(subset(fulloriginal, distance > 20 & distance <= 200)$distance, main=c("20km < d < 200km"), breaks=20, xlab=c("Distance (km)"))
# 
# par(mfrow=c(2,2), mar=c(4,2,2,1))
# hist(subset(fulloriginal, distance > 0 & distance <= 100)$distance, main=c("0km < d < 100km"), xlab=c("Distance (km)"))
# hist(subset(fulloriginal, distance > 10 & distance <= 100)$distance, main=c("10km < d < 100km"), xlab=c("Distance (km)"))
# hist(subset(fulloriginal, distance > 20 & distance <= 50)$distance, main=c("20km < d < 100km"), xlab=c("Distance (km)"))
# hist(subset(fulloriginal, distance < 10)$distance, main=c("d < 10km"), xlab=c("Distance (km)"))
# 
# par(mfrow=c(2,2))
# hist(subset(fulloriginal, distance > 20 & distance <= 200)$distance, main=c("20km < d < 200km"), xlab=c("Distance (km)"), breaks=10)
# hist(subset(fulloriginal, distance > 20 & distance <= 200)$distance, main=c("20km < d < 200km"), xlab=c("Distance (km)"), breaks=20)
# hist(subset(fulloriginal, distance > 20 & distance <= 100)$distance, main=c("20km < d < 100km"), xlab=c("Distance (km)"), breaks=10)
# hist(subset(fulloriginal, distance > 20 & distance <= 100)$distance, main=c("20km < d < 100km"), xlab=c("Distance (km)"), breaks=20)
