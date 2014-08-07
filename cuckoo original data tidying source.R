##########################################################
#
#  Cuckoo ORIGINAL data - tidying source file, adding new variables (distance moved between points, bearing, movement groups, and movement types)
#
#  Samantha Franks
#	15 Nov 2013
# 4 Dec 2013
# 23 Dec 2013
#
#
##########################################################

library(NCStats)
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


###-----------------------------------------------###
#         Load cuckoo data & tidying
###-----------------------------------------------###

# cuckoo coordinate locations are in Lat-Long and UTM coordinates and WGS84 datum

setwd("C:/Users/samf/Documents/Git/cuckoos/data")
#d <- read.csv("BTO Cuckoo migration study - raw data + new variables 15112013.csv", header=T) # cuckoo movement

d <- read.csv("cuckoos raw data 20131118_reformatted_20131213.csv", header=T)

# cuckoo <- read.csv("cuckoo movements.csv", header=T) # cuckoo movement
# 
# reqcuckoo <- with(cuckoo, data.frame(timestamp, julian=julian.date, year, month, day, id, lat, long, duration, days.since.tagging, name, country, continent, movement.number))

# $filter refers to the values identified as best in the transmission cycle by the Douglas filter applied in Movebank: TRUE values are the best locations
# $loc.class refers to the Argos location class: from best to worst, 3,2,1,0,A,B,Z. Alpha location classes have no associated errors, so remove. Location class=0 have gigantic major axis errors, so also remove these. Use only location classes 1:3
# error.major = error radius along the major axis of the ellipse
# error.minor = error radius along the minor axis of the ellipse
# test location refers to points that were used as either reference points or tested outside Chris Hewson's house in Cambridge (long/lat = 0.125/52.135)

# d2 <- with(d, data.frame(name=individual.local.identifier, id=event.id, filter=visible, tcycle=Transmission.cycle, timestamp, julian, year, month, day, long=location.long, lat=location.lat, loc.class=argos.lc, error.major=argos.semi.major, error.minor=argos.semi.minor))

d2 <- with(d, data.frame(name, id, tcycle=transmission.cycle, timestamp, location_date, julian, year, month, day, hour, long=longitude, lat=latitude, loc.class=location_class, error.major=semi_major_axis, error.minor=semi_minor_axis, ell.orientation=ellipse_orientation, sensor1=X_1, sensor2=X_2, sensor3=X_3, sensor4=X_4,test_location))

# d2$filter <- as.factor(d2$filter) # define the filter as a factor rather than logical

d2$tcycle <- as.factor(d2$tcycle)

# removes locations outside Afro-Palearctic (reference points in North America plus any really weird points in the Atlantic)
d3 <- Subset(d2, long > -30)

d4 <- Subset(d3, loc.class=="1" | loc.class=="2" | loc.class=="3") # subset by best location classes

d5 <- Subset(d4, test_location=="N")

d6 <- Subset(d5, name!="" & name!="115592" & name!="115601" & name!="121792")

d7 <- Subset(d6, error.major>0 & error.minor>0)

unique_locdates <- which(duplicated(d7$location_date)==FALSE)

d8 <- d7[unique_locdates,]

d9 <- na.omit(d8) # complete cases only

###-----------------------------------------------###
#         Add distance, bearing, mgroup and mtype to original data
###-----------------------------------------------###


byname <- d9$name

newdat <- split(d9, list(byname)) # split dataset by individual bird

trimdata <- function(dat) { # distance, bearing, etc function
  
  # for each element of data, do the following
  
  # order dat according to location_date - not needed if dataset has already been ordered in Excel by location_date
  # datorder <- dat[order(dat$location_date),]
  
  smallerror <- Subset(dat, error.major <= 5000)
  
  
  ### --- DISTANCE & BEARING CALCULATION --- ###
  
  # calculate distances and bearings   
  # distbearmvmt <- function(input){ # requires 2 variables called newlongs & newlats or long & lat or centroidlong & centroidlat
  # use newlongs & newlats to calculate distances based on resampled points
  # use long & lat to calculate distances based on original measured location point
  
  input <- smallerror
  
  dist.cuckoo <- input[,c("long","lat")]
  coordinates(dist.cuckoo) <- c("long", "lat")
  proj4string(dist.cuckoo) <- CRS("+proj=longlat +datum=WGS84")
  
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
write.csv(fulloriginal, "allbirdfinal - original data.csv", row.names=F)

setwd("C:/Users/samf/Documents/Git/cuckoos/data/original data + distance, movement groups, etc")
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
