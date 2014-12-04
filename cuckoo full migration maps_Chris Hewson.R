##########################################################
#
# Cuckoo migration tracks for Chris Hewson
#
# 1. plot stopovers from breeding location to winter solstice each individual on same plot
#     - different colours for SE/SW strategy (group Balkans as SE)
#     - stopover point size relative to between first and last transmission from same location
# 
# 2. plot repeated journeys for individuals with > 1 track on same plot
#     - different colours for each individual
#
# 12 Feb 2014
# 28 Nov 2014 - updated to include 2014 cohort of birds
#	Samantha Franks
#	
#
##########################################################


######-------------------   NOTES -----------------##########
# longs and lats from Movebank are in WGS84 GCS

rm(list=ls())

### Load packages

library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
#library(shapefiles) # allows you to read and write ESRI shapefiles
library(rgeos)
library(raster)
library(geosphere)
library(rworldmap) # getMap function
library(rworldxtra)
library(reshape)
library(plyr)
library(dismo)

options(scipen=6)

### -------- Working directories --------- ###

Mac <- FALSE

if (!Mac) parentwd <- parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")

datawd <- paste(parentwd, "data", sep="/")
scriptswd <- paste(parentwd, "scripts", sep="/")
outputwd <- paste(parentwd, "output", "Chris Hewson maps", sep="/")

##################################################################
#
#   MAP 1: STOPOVERS & MIGRATION TRACKS
#
##################################################################

###-----------------------------------------------------------###
#         1. CREATE MAP BASE LAYERS
###-----------------------------------------------------------###


### ------- Europe/Africa national boundaries layer ------- ###

dat.world <- getMap(resolution="high")

Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa" | dat.world$REGION=="Asia"),]
#Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia"),]
Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]

Eur.Afr <- spTransform(Eur.Afr,CRS("+proj=longlat +datum=WGS84")) # transform map to exactly the same projection as cuckoo newlongs/newlats data

Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]

clipEurAfr <- extent(-21,60,-20,64)
EurAfrmap <- crop(Eur.Afr,clipEurAfr)


###-----------------------------------------------------------###
#         2. SET UP CUCKOO DATA
###-----------------------------------------------------------###

setwd(datawd)

### which dataset to use
# cuckoo <- read.csv("Chris Hewson maps_all cuckoos best location data_05022014.csv", header=T) # data used for map production Feb 2014
d <- read.csv("BTO Cuckoo migration study_Movebank_20141128.csv", header=T) # data to use for map production 28/11/2014, includes 2014 cohort, uses data downloaded straight from Movebank as .csv without any other formatting

# remove Wistman, Karma, New Forest 1 (these don't get out of UK)
rm.names <- c("Wistman","Karma","New Forest 1")

# run source code to add dates, calculate distances, movement groups and movement types
# needs to be passed an object called "d"
# produces object called "fulloriginal"
source(paste(scriptswd, "cuckoo Movebank data tidying source_new data 2014.R", sep="/"))

# run source code to calculate stopover duration
# needs to be passed an object called "present.date", which has a variable called datetime (of class chron), and which is ordered according to name, then by date
# outputs an object called present.LOS (the dataset with a new column with stopover duration calculated)
# also outputs an object called LOStable.ref, which is a quick reference table with each individual's mgroup and corresponding stopover duration
present.date <- fulloriginal[order(fulloriginal$name, fulloriginal$datetime),]
source(paste(scriptswd, "source code to add stopover duration data for analysis.R", sep="/"))

setwd(datawd)
# add some extra data about migration strategy
strategy.dat <- read.csv("cuckoo migratory strategy and Sahara crossing success 2014.csv", header=T)
strategy.dat <- subset(strategy.dat, select=c("name2","year.tagged","migratory.strategy","number.of.years","successful.Sahara.crossing"))
strategy.dat <- rename(strategy.dat, c("name2"="name", "migratory.strategy"="strat","number.of.years"="no.yrs", "successful.Sahara.crossing"="Sahara.success"))
# strategy.dat <- subset(strategy.dat, name!="Idemili" & name!="Karma" & name!="Wistman")

# merge migration strategy info with bird data
fulldat <- merge(present.LOS, strategy.dat, by="name")

# cuckoo2 <- read.csv("all cuckoos_best locations_20140331.csv", header=T)
# lloyd2 <- subset(cuckoo2, individual.local.identifier=="Lloyd")
# lloyd2 <- droplevels(lloyd2)


###--- Autumn only ---###

juliandate <- as.numeric()

# cannot get julian() to work as a function on each row and spit out the right value, so for now is written as a slow loop!
for (i in 1:nrow(fulldat)) {
  juliandate[i] <- julian(x=months(fulldat$datetime[i]), d=days(fulldat$datetime[i]), y=years(fulldat$datetime[i]), c(day=1, month=1, year=years(fulldat$datetime[i])))
}

fulldat2 <- cbind(fulldat, juliandate)

autumnonly <- subset(fulldat2, juliandate >= 131 & juliandate <=365 )
# autumnonly <- Subset(fulldat, julian >= 131 & julian <=356)
# Clip dataset to only southward migration breeding May 10 to Dec 21 (julian >=131 and <=356)

###--- remove stopover points in Africa in May (Tor 15) ---###

remove <- which(autumnonly$name=="Tor" & autumnonly$mgroup==15)
autumnonly <- autumnonly[-remove,]

###--- remove very long stopover points in Spain (Jake 42 and 43) ---###

remove <- which(autumnonly$name=="Jake" & autumnonly$mgroup >= 42)
autumnonly <- autumnonly[-remove,]

# change Whortle's 2014 autumn migration strategy to SE and name associated to Whortle2
autumnonly[autumnonly$name=="Whortle" & autumnonly$year=="2014","strat"] <- "SE"

autumnonly$name <- factor(autumnonly$name, levels=c(levels(autumnonly$name), "Whortle 2"))
autumnonly[autumnonly$name=="Whortle" & autumnonly$year=="2014", "name"] <- "Whortle 2"

###########################
###########################

###--- No migratory movements ---###

stopovers <- subset(autumnonly, mtype!="M")

###--- Pull out individuals with multiple years & subset data ---###

multiyeartable <- ifelse(table(autumnonly$name, autumnonly$year) > 0, 1, 0)
multiyearnames <- c(names(which(apply(multiyeartable, 1, sum) > 1)), "Whortle", "Whortle 2")

multiyear <- autumnonly[grep(paste(multiyearnames, collapse="|"), autumnonly$name),]
multiyear <- droplevels(multiyear)
multiyear.stopovers <- subset(multiyear, mtype!="M")
multiyear.stopovers <- droplevels(multiyear.stopovers)


###-----------------------------------------------------------###
#         3. CREATE SPATIAL LINES CONNECTING YEAR *ALL* STOPOVERS
###-----------------------------------------------------------###

map2 <- TRUE

if (map2 == TRUE) {
  
  # split dataset
  byname <- multiyear$name
  forlines <- split(multiyear, list(byname))
}

if (map2 == FALSE) {
  
  # split dataset
  byname <- autumnonly$name
  forlines <- split(autumnonly, list(byname))
}

allbirdlines <- list()

for (n in 1:length(forlines)) {
  
  dat <- forlines[[n]]
  dat <- droplevels(dat)
  dat$year <- as.factor(dat$year)
  
  yearlevel <- levels(dat$year)
  print(paste(levels(dat$name), length(yearlevel)))
  
  if (length(yearlevel) == 1) {
    
    line1 <- Line(cbind(dat$long, dat$lat))
    birdlines <- Lines(list(line1), ID = levels(dat$name))
    
  } else if (length(yearlevel) == 2) {
    
    year1 <- dat[dat$year == yearlevel[1],]
    year2 <- dat[dat$year == yearlevel[2],]
    
    line1 <- Line(cbind(year1$long, year1$lat))
    line2 <- Line(cbind(year2$long, year2$lat))
    birdlines <- Lines(list(line1,line2), ID = levels(dat$name))
    
  } else if (length(yearlevel) == 3) {
    
    year1 <- dat[dat$year == yearlevel[1],]
    year2 <- dat[dat$year == yearlevel[2],]
    year3 <- dat[dat$year == yearlevel[3],]
    
    line1 <- Line(cbind(year1$long, year1$lat))
    line2 <- Line(cbind(year2$long, year2$lat))
    line3 <- Line(cbind(year3$long, year3$lat))    
    birdlines <- Lines(list(line1,line2,line3), ID = levels(dat$name))
    
  } else {
    
    year1 <- dat[dat$year == yearlevel[1],]
    year2 <- dat[dat$year == yearlevel[2],]
    year3 <- dat[dat$year == yearlevel[3],]
    year4 <- dat[dat$year == yearlevel[4],]
    
    line1 <- Line(cbind(year1$long, year1$lat))
    line2 <- Line(cbind(year2$long, year2$lat))
    line3 <- Line(cbind(year3$long, year3$lat))
    line4 <- Line(cbind(year4$long, year4$lat)) 
    birdlines <- Lines(list(line1,line2,line3,line4), ID = levels(dat$name))
  }
  
  allbirdlines[[n]] <- birdlines
  
}

line.routes <- SpatialLines(allbirdlines)
proj4string(line.routes) <- CRS("+proj=longlat +datum=WGS84")

if (!map2) {
  mig.routes <- SpatialLinesDataFrame(line.routes, data.frame(unique(autumnonly[,c("name","strat")]), row.names=sapply(slot(line.routes, 'lines'), function(i) slot(i, 'ID'))))
}

if (map2) {
  mig.routes <- SpatialLinesDataFrame(line.routes, data.frame(unique(multiyear[,c("name","strat")]), row.names=sapply(slot(line.routes, 'lines'), function(i) slot(i, 'ID'))))
}

###-----------------------------------------------------------###
#         4. CREATE MAIN STOPOVER LOCATION POINTS EACH CUCKOO
###-----------------------------------------------------------###

###########
#
# Extract the centroid of the MCP for stopovers >= 2 transmissions
# Use dataframe(stopovers), which is only C and S mtypes
#
###########

if (map2 == TRUE) {
  
  # split dataset
  byname <- multiyear.stopovers$name
  newdat <- split(multiyear.stopovers, list(byname))
}

if (map2 == FALSE) {
  
  # split dataset
  byname <- stopovers$name
  newdat <- split(stopovers, list(byname))
}


# for each movement group of each cuckoo with >= 3 points, draw a polygon around the mgroup's points and extract the centroid
# if 2 points, then draw a line and calculate the mid-point
# choose points for MCP such that 5% of outliers are excluded

mig.stops <- list() # create list of all centroid points (or just single points) for each cuckoo's southward migration

for (n in 1:length(newdat)) {
  
  dat <- newdat[[n]]
  dat <- droplevels(dat)
  dat$mgroup <- as.factor(dat$mgroup)
  
  centroids <- matrix(data=NA, nrow=length(levels(dat$mgroup)), ncol=3)
    
  coordinates(dat) <- c("long","lat")
  proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
  
  for (i in 1:length(levels(dat$mgroup))){
  
    mgrouplevel <- levels(dat$mgroup)[i]
    usedat <- dat[dat$mgroup==mgrouplevel,]
    
    if (nrow(usedat) >= 3) {
      
      # place "centroid" point of stopover in the centre of the bbox of the spatial points
      centroidlong <- bbox(usedat)["long","min"] + ((bbox(usedat)["long","max"]-bbox(usedat)["long","min"])/2)
      centroidlat <- bbox(usedat)["lat","min"] + ((bbox(usedat)["lat","max"]-bbox(usedat)["lat","min"])/2)
      
#       r1 <- rbind(coordinates(usedat),coordinates(usedat[1,]))
#       mg.poly <- Polygon(r1) # this fails to make a closed polygon for Martin mgroup 11 for some reason, gives ridiculous centroid
      centroids[i,] <- cbind(centroidlong, centroidlat, as.numeric(mgrouplevel))
      
    } else if (nrow(usedat) == 2) {
      
      midpt <- midPoint(coordinates(usedat)[1,],coordinates(usedat)[2,])
      centroids[i,] <- cbind(midpt[1], midpt[2], as.numeric(mgrouplevel))
      
    } else {
      
      centroids[i,] <- cbind(coordinates(usedat), as.numeric(mgrouplevel))
    }
  }
  
  centroids <- data.frame(centroids)
  centroids <- cbind(rep(levels(dat$name), nrow(centroids)), centroids)
  colnames(centroids) <- c("name", "long","lat","mgroup")
  
  centroids2 <- merge(centroids, unique(dat@data[,c("name","mgroup","LOS","year","strat")]), by=c("name","mgroup"))  
  centroids2 <- centroids2[order(centroids2$mgroup),]
    
  centroids2$mgroup <- as.factor(centroids2$mgroup)
  
  mig.stops[[n]] <- centroids2
  
}

# ###-----------------------------------------------------------###
# #         5. ADD STOPOVER DURATION & YEAR TO LOCATION POINTS
# ###-----------------------------------------------------------###
# 
# # for each cuckoo
# 
# fullroute <- list()
# 
# for (n in 1:length(newdat)) {
#   
#   dat <- newdat[[n]]
#   dat <- droplevels(dat)
#   dat$mgroup <- as.factor(dat$mgroup)
#   dat$year <- as.factor(dat$year)
#   
#   # for each level of mgroup:
#   #   1. sum duration.days
#   #   2. extract year
#   
#   year.LOS <- data.frame(matrix(NA, nrow=nrow(mig.stops[[n]]), ncol=2, dimnames=list(c(), c("year","LOS"))))
#   
#   for (i in 1:length(levels(dat$mgroup))){
#     
#     mgrouplevel <- levels(dat$mgroup)[i]
#     usedat <- dat[dat$mgroup==mgrouplevel,]
#     usedat <- droplevels(usedat)
#     
#     year <- levels(usedat$year)
#     
#     LOS <- sum(usedat$duration.days, na.rm=TRUE)
#     
#     year.LOS[i,] <- c(year, LOS)
#         
#   }
#   
#   fullroute[[n]] <- data.frame(mig.stops[[n]], year.LOS)
#   fullroute[[n]]$LOS <- as.numeric(fullroute[[n]]$LOS)
#   fullroute[[n]]$year <- as.factor(fullroute[[n]]$year)
#   
# }


###-----------------------------------------------------------###
#         6. ADD EXTRA VARIABLES ON MIG STRATEGY, LOS BIN
###-----------------------------------------------------------###

# re-concatenate individual bird list elements into a single dataframe with all individuals

allbirds <- do.call(rbind, mig.stops) # rbind all list elements

# allbirds$LOS <- as.numeric(allbirds$LOS)
# allbirds$year <- as.factor(allbirds$year)

if (map2 == FALSE) {
  # remove stopovers < 3 days
  remove <- which(allbirds$LOS < 3 & allbirds$mgroup!="1")
  longstops <- allbirds[-remove,]
}

if (map2 == TRUE) {
#   longstops <- allbirds
  remove <- which(allbirds$LOS < 3 & allbirds$mgroup!="1")
  longstops <- allbirds[-remove,]
}

# ###--- assign migratory strategy (map 1 only) ---###
# 
# if (map2 == FALSE) {
#   
#   migstrategy <- function(name) {
#     if (name=="BB" | name=="Chance" | name=="Chris" | name=="David" | name=="Indy" | name=="Iolo" | name=="Kasper" | name=="Livingstone" | name=="Lloyd" | name=="Martin" | name=="Mungo" | name=="Patch" | name=="Roy" | name=="Sussex" | name=="Sussex" | name=="Tor" | name=="Wallace" | name=="Waller") {
#       strategy <- "SE"
#     } else {
#       strategy <- "SW"
#     }
#     return(strategy)
#   }
#   
#   mig.strategy <- sapply(longstops$name, migstrategy)
#   
# }


###--- convert LOS into binned categorical variable ---###
###--- create a cex vector for LOS ---###

groupLOS <- function(LOSvariable){
  
#   if (LOSvariable <= 5) {
#     bin <- 0.3
#   } else if (LOSvariable > 5 & LOSvariable <= 10) {
#     bin <- 0.5
#   } else if (LOSvariable > 10 & LOSvariable <= 30) {
#     bin <- 0.7
#   } else {
#     bin <- 0.9
#   }
#   return(bin)
  
#   if (LOSvariable <= 7) {
#     bin <- 0.3
#   } else if (LOSvariable > 7 & LOSvariable <= 14) {
#     bin <- 0.6
#   } else if (LOSvariable > 14 & LOSvariable <= 35) {
#     bin <- 0.9
#   } else {
#     bin <- 1.2
#   }
#   return(bin)
  
#   if (LOSvariable <= 14) {
#     bin <- 0.5
#   } else {
#     bin <- 1.2
#   }
#   return(bin)
  
  #for (i in 1:nrow(longstops)) {
    
    #LOSvariable <- longstops$LOS[i]
    
    if (LOSvariable <= 14) {
      bin <- 0.3
    } else if (LOSvariable > 14 & LOSvariable <= 30) {
      bin <- 0.7
    } else if (is.na(LOSvariable)) {
      bin <- 0
    } else {
      bin <- 1.1
    }
    return(bin)
  
}

# if any capture locations have NA for LOS, then set LOS to 0
longstops$LOS <- ifelse(is.na(longstops$LOS), 0, longstops$LOS)
LOS.bin <- sapply(na.omit(longstops$LOS), groupLOS)

###--- FULL DATAFRAME ---###

if (map2 == FALSE) {
  cuckoo.pts <- data.frame(longstops, LOS.bin)
}

if (map2 == TRUE) {
  cuckoo.pts <- data.frame(longstops, LOS.bin)
}

###--- convert long/lat to coordinates in WGS84 ---###

coordinates(cuckoo.pts) <- c("long","lat")
proj4string(cuckoo.pts) <- CRS("+proj=longlat +datum=WGS84")

###--- extract country / continent names ---###

geoinfo <- over(cuckoo.pts, EurAfrmap)[,c("REGION","ADMIN")]
colnames(geoinfo) <- c("continent", "country")
geoinfo <- droplevels(geoinfo)

toplot.pts <- data.frame(cuckoo.pts, geoinfo)

###--- convert long/lat to coordinates in WGS84 ---###

coordinates(toplot.pts) <- c("long","lat")
proj4string(toplot.pts) <- CRS("+proj=longlat +datum=WGS84")

###--- points to remove because dead bird, in sea etc ---###

if (map2 == FALSE) {
  
  # TO CHECK
  # BB mgroup 29 - on Wadden Sea island, ok - toplot.pts[8,]
  # Martin mgroup 1 and 2 - capture occasion but in sea. . .not sure - toplot.pts[185:186,]
  # Waller mgroup 1 - middle of North Sea, remove
  
  remove <- which(toplot.pts$name=="Waller" & toplot.pts$mgroup==1)
  toplot.pts <- toplot.pts[-remove,]
  
}

###--- add continent/country labels for NA points --###

if (map2 == FALSE) {
  which(is.na(toplot.pts$country))
  
  toplot.pts[which(is.na(toplot.pts$continent)), "continent"] <- "Europe"
  toplot.pts[which(is.na(toplot.pts$country)), "country"] <- c("Netherlands", "United Kingdom", "United Kingdom")
  
}

# if (map2 == TRUE) {
#   which(is.na(toplot.pts$country))
#   
#   toplot.pts[which(is.na(toplot.pts$continent)), "continent"] <- "Europe"
#   toplot.pts[which(is.na(toplot.pts$country)), "country"] <- c("Netherlands")
#   
# }
  

###--- remove stopover points in UK that are breeding & > 10d ---###

remove <- which(toplot.pts$country=="United Kingdom")
toplot.pts <- toplot.pts[-remove,]

if (!map2) {
  ref <- merge(autumnonly, toplot.pts@data, by=c("name","mgroup","LOS","strat","year"))
}

if (map2) {
  ref <- merge(multiyear, toplot.pts@data, by=c("name","mgroup","LOS","strat","year"))
}

########## EMSWORTHY MAY BE CODED AS INCORRECT STRATEGY, CHECK WITH CHRIS (starts SW, but moves onto SE from Italy, track makes a bend over Tunisia)
########## WHORTLE CHANGES STRATEGY IN YEAR 2 (2014) - starts as SW bird, moves onto SE

###-----------------------------------------------------------###
#         7. PLOT MAPS
###-----------------------------------------------------------###

##################### MAP 1 ########################

if (map2 == FALSE) {
  
  ###--- create a colour vector for migratory strategy ---###
  
  colstrategy <- factor(toplot.pts@data$strat)
  colstrategy <- revalue(colstrategy, c("SW"="blue", "SE"="red"))
  colstrategy.lines <- factor(mig.routes@data$strat)
  colstrategy.lines <- revalue(colstrategy.lines, c("SW"="blue", "SE"="red"))
  
  setwd(outputwd)
  
  tiff("all cuckoos southward migration and LOS_final.tiff",res=300,height=2200,width=2000,units="px")
  
  par(mar=c(0,0,0,0))
  plot(EurAfrmap, border="grey60")
  
  
#   x <- subset(toplot.pts, name=="lloyd")
#   y <- subset(mig.routes, name=="lloyd")
#   colstrategy <- factor(x@data$strat)
#   colstrategy <- revalue(colstrategy, c("SW"="blue", "SE"="red"))
#   colstrategy.lines <- factor(y@data$strat)
#   colstrategy.lines <- revalue(colstrategy.lines, c("SW"="blue", "SE"="red"))
#   
#   # plot major (> 3 day) stopovers
#   plot(x, pch=16, col=as.character(colstrategy), cex=toplot.pts$LOS.bin, add=T)
#   
#     # plot migration routes
#   plot(y, col=as.character(colstrategy.lines), lwd=0.6, add=T)
  
  # plot major (> 3 day) stopovers
  plot(toplot.pts, pch=16, col=as.character(colstrategy), cex=toplot.pts$LOS.bin, add=T)
  
  
  # plot migration routes
  plot(mig.routes, col=as.character(colstrategy.lines), lwd=0.6, add=T)
  
  legend(35, 62, legend=c("< 14 days", "14-30 days", "> 30 days", "SE strategy", "SW strategy"), pch=16, cex=0.7, pt.cex=c(0.3, 0.7, 1.1, 0, 0), text.col=c("black","black","black","red","blue"), bty="n")
  
  dev.off()
  
}

##################### MAP 2 ########################

if (map2 == TRUE) {
  
  ###--- create a colour vector for each individual ---###
  
  # 2 blacks for 2 Whortle years
  colours <- c("blue", "deeppink", "forestgreen", "cyan3", "darkorchid1","red","orange","olivedrab3","peachpuff3","black", "black")
  
  colstrategy <- factor(toplot.pts@data$name, labels=colours)
  colstrategy.lines <- factor(mig.routes@data$name, labels=colours)
  
  setwd(outputwd)
  
  tiff("multiyear cuckoos southward migration and LOS_final.tiff",res=300,height=2200,width=2000,units="px")
  
  par(mar=c(0,0,0,0))
  plot(EurAfrmap, border="grey60")
  
  # plot major (> 3 day) stopovers
  plot(toplot.pts, pch=16, col=as.character(colstrategy), cex=toplot.pts$LOS.bin, add=T)
  
  # plot migration routes
  plot(mig.routes, col=as.character(colstrategy.lines), lwd=0.9, add=T)
  
  #   legend(37, 65, legend=c("< 14 days", "14-30 days", "> 30 days", levels(toplot.pts@data$name)), pch=16, cex=0.7, pt.cex=c(0.3, 0.7, 1.1, 0, 0, 0, 0, 0,0,0,0,0,0), text.col=c("black", "black", "black", colours), bty="n")
  
  legend(37, 60, legend=c("< 14 days", "14-30 days", "> 30 days"), pch=16, cex=0.7, pt.cex=c(0.3, 0.7, 1.1), text.col=c("black", "black", "black"), bty="n")
  
  
  legend(-15, 5, legend=c(levels(toplot.pts@data$name)[-which(levels(toplot.pts@data$name)=="Whortle 2")]), pch=16, cex=1, pt.cex=c(0, 0, 0, 0, 0,0,0,0,0,0), text.col=c(colours), bty="n")
  
  
  dev.off()  
}

# 
# ############## MAP 2 (version 2, separate maps/bird ########
# 
# if (map2 == TRUE) {
#   
#   ###--- create a colour vector for each individual ---###
#   
#   colours <- c("blue", "deeppink", "forestgreen", "cyan3", "darkorchid1")
#   
#   colstrategy <- factor(toplot.pts@data$name, labels=colours)
#   colstrategy.lines <- factor(mig.routes@data$name, labels=colours)
#   
#   setwd(outputwd)
#   
#   tiff("multiyear cuckoos southward migration and LOS_v2.tiff",res=300,height=2000,width=2400,units="px")
#   
#   lo <- layout(rbind(c(1,2,3),c(4,5,6)))
#   #layout.show(lo)
#   par(mar=c(0,0,0,0))
#   
#   for (i in 1:length(levels(toplot.pts$name))){
#     
#     name <- levels(toplot.pts$name)[i]
#     
#     plot(EurAfrmap, border="grey60")
#     plot(toplot.pts[toplot.pts$name==name,], pch=16, col=colours[i], cex=toplot.pts$LOS.bin, add=T)
#     
#     text(-5,-5, name, font=2, cex=1.5)
#   }
#   
#   
#   legend(80, 30, legend=c("< 14 days", "14-30 days", "> 30 days"), pch=16, cex=1.5, pt.cex=c(0.7, 1.1, 1.5), text.col=c("black", "black", "black"), xpd=NA, bty="n")
#   
#   dev.off()
#   
#     
# }




######################
# 
# windows(10,10)
# par(mfrow=c(5,2))
# hist(allbirds[allbirds$LOS < 10,"LOS"])
# hist(allbirds[allbirds$LOS < 20 & allbirds$LOS > 10,"LOS"])
# hist(allbirds[allbirds$LOS < 30 & allbirds$LOS > 20,"LOS"])
# hist(allbirds[allbirds$LOS < 40 & allbirds$LOS > 30,"LOS"])
# hist(allbirds[allbirds$LOS < 50 & allbirds$LOS > 40,"LOS"])
# hist(allbirds[allbirds$LOS < 60 & allbirds$LOS > 50,"LOS"])
# hist(allbirds[allbirds$LOS < 70 & allbirds$LOS > 60,"LOS"])
# hist(allbirds[allbirds$LOS < 80 & allbirds$LOS > 70,"LOS"])
# hist(allbirds[allbirds$LOS < 90 & allbirds$LOS > 80,"LOS"])
# hist(allbirds[allbirds$LOS > 90,"LOS"])

######################################################################################################################################################################################################

##################################################################
#
#     MAP 2: TRACKS FOR INDIVIDUALS WITH MULTIPLE YEARS
#
##################################################################

