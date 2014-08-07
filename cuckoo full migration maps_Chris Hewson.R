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
#	Samantha Franks
#	
#
##########################################################

rm(list=ls())

### Load packages

library(adehabitat) # 
library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
#library(shapefiles) # allows you to read and write ESRI shapefiles
library(rgeos)
library(raster)
library(geosphere)
library(adehabitatHR)
library(rworldmap) # getMap function
library(rworldxtra)
library(NCStats)

### -------- Working directories --------- ###

datawd <- c("C:/Users/samf/Documents/Git/cuckoos/data")

outputwd <- ("C:/Users/samf/Documents/Git/cuckoos/output/Chris Hewson maps")

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

cuckoo <- read.csv("Chris Hewson maps_all cuckoos best location data_05022014.csv", header=T)

cuckoo2 <- read.csv("all cuckoos_best locations_20140331.csv", header=T)
lloyd2 <- subset(cuckoo2, individual.local.identifier=="Lloyd")
lloyd2 <- droplevels(lloyd2)

strategy.dat <- read.csv("cuckoo migratory strategy.csv", header=T)

###### ADD MOVEMENT TYPE + MOVEMENT GROUP ######

### --- MOVEMENT GROUPS --- ###

mgroup <- rep(NA,nrow(cuckoo))

for (n in 1:nrow(cuckoo)){
  if (is.na(cuckoo$distance[n])) {
    mgroup[n] <- 1
  } else if (cuckoo$distance[n] <= 25) {
    mgroup[n] <- mgroup[n-1]
  } else {mgroup[n] <- mgroup[n-1] + 1}
}

### --- MOVEMENT TYPES --- ###

mtype <- rep(NA, nrow(cuckoo))

for (n in 1:nrow(cuckoo)){
  if (is.na(cuckoo$distance[n])) {
    mtype[n] <- "C"
  } else if (mgroup[n] != mgroup[n-1] & mgroup[n] != mgroup[n+1]) {
    mtype[n] <- "M"
  } else {
    mtype[n] <- "S"
  }
}

completedata <- data.frame(cuckoo,mgroup,mtype)

newcuckoo <- with(completedata, data.frame(id, name, event.id, timestamp, year, month, day, hour, julian, duration.days, duration.hrs, long.WGS84, lat.WGS84, direction, distance, mgroup, mtype))

newcuckoo$year <- as.factor(newcuckoo$year)
newcuckoo$mgroup <- as.factor(newcuckoo$mgroup)

###--- remove birds with no mig strategy

# remove Karma (no migratory strategy)
newcuckoo <- Subset(newcuckoo, name!="Karma")

###--- Autumn only ---###

autumnonly <- Subset(newcuckoo, julian >= 131 & julian <=356)
# Clip dataset to only southward migration breeding May 10 to Dec 21 (julian >=131 and <=356)

###########################
########################### test Lloyd

lloyd <- subset(autumnonly, name=="Lloyd")
lloyd2 <- droplevels(lloyd2)

coordinates(lloyd2) <- c("long.WGS84","lat.WGS84")
coordinates(lloyd2) <- c("location.long","location.lat")

proj4string(lloyd2) <- CRS("+proj=longlat +datum=WGS84")

plot(EurAfrmap)
plot(lloyd2, col="blue", pch=16, add=T)

###########################
###########################

###--- No migratory movements ---###

stopovers <- Subset(autumnonly, mtype!="M")

###--- Pull out individuals with multiple years ---###

ind.by.year <- as.data.frame.matrix(table(autumnonly$name, autumnonly$year))
ind.by.year2 <- data.frame(name=rownames(ind.by.year), ind.by.year)
rownames(ind.by.year2) <- c(1:29)

multiyear <- Subset(autumnonly, name=="BB" | name=="Chance" | name=="Chris" | name=="David" | name=="Lyster")

multiyear.stopovers <- Subset(multiyear, mtype!="M")

###-----------------------------------------------------------###
#         3. CREATE SPATIAL LINES CONNECTING YEAR *ALL* STOPOVERS
###-----------------------------------------------------------###

map2 <- FALSE

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
  
  yearlevel <- levels(dat$year)
  
  if (length(yearlevel) == 1) {
    
    line1 <- Line(cbind(dat$long.WGS84, dat$lat.WGS84))
    birdlines <- Lines(list(line1), ID = levels(dat$name))
    
  } else if (length(yearlevel) == 2) {
    
    year1 <- dat[dat$year == yearlevel[1],]
    year2 <- dat[dat$year == yearlevel[2],]
    
    line1 <- Line(cbind(year1$long.WGS84, year1$lat.WGS84))
    line2 <- Line(cbind(year2$long.WGS84, year2$lat.WGS84))
    birdlines <- Lines(list(line1,line2), ID = levels(dat$name))
    
  } else {
    
    year1 <- dat[dat$year == yearlevel[1],]
    year2 <- dat[dat$year == yearlevel[2],]
    year3 <- dat[dat$year == yearlevel[3],]
    
    line1 <- Line(cbind(year1$long.WGS84, year1$lat.WGS84))
    line2 <- Line(cbind(year2$long.WGS84, year2$lat.WGS84))
    line3 <- Line(cbind(year3$long.WGS84, year3$lat.WGS84))    
    birdlines <- Lines(list(line1,line2,line3), ID = levels(dat$name))
    
  }
  
  allbirdlines[[n]] <- birdlines
  
}

line.routes <- SpatialLines(allbirdlines)
proj4string(line.routes) <- CRS("+proj=longlat +datum=WGS84")

mig.routes <- SpatialLinesDataFrame(line.routes, data.frame(strategy.dat, row.names=levels(autumnonly$name)))


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
  
  centroids <- matrix(data=NA, nrow=length(levels(dat$mgroup)), ncol=3)
    
  coordinates(dat) <- c("long.WGS84","lat.WGS84")
  proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
  
  for (i in 1:length(levels(dat$mgroup))){
  
    mgrouplevel <- levels(dat$mgroup)[i]
    usedat <- dat[dat$mgroup==mgrouplevel,]
    
    if (nrow(usedat) >= 3) {
      r1 <- rbind(coordinates(usedat),coordinates(usedat[1,]))
      mg.poly <- Polygon(r1)
      centroids[i,] <- cbind(mg.poly@labpt[1], mg.poly@labpt[2], as.numeric(mgrouplevel))
      
    } else if (nrow(usedat) == 2) {
      
      midpt <- midPoint(coordinates(usedat)[1,],coordinates(usedat)[2,])
      centroids[i,] <- cbind(midpt[1], midpt[2], as.numeric(mgrouplevel))
      
    } else {
      
      centroids[i,] <- cbind(coordinates(usedat), as.numeric(mgrouplevel))
    }
  }
  
  centroids <- data.frame(centroids)
  centroids <- cbind(rep(levels(dat$name), nrow(centroids)), centroids)
  colnames(centroids) <- c("name", "long.WGS84","lat.WGS84","mgroup")
  centroids$mgroup <- as.factor(centroids$mgroup)
  
  mig.stops[[n]] <- centroids
  
}
 
###-----------------------------------------------------------###
#         5. ADD STOPOVER DURATION & YEAR TO LOCATION POINTS
###-----------------------------------------------------------###

# for each cuckoo

fullroute <- list()

for (n in 1:length(newdat)) {
  
  dat <- newdat[[n]]
  dat <- droplevels(dat)
  dat$mgroup <- as.factor(dat$mgroup)
  dat$year <- as.factor(dat$year)
  
  # for each level of mgroup:
  #   1. sum duration.days
  #   2. extract year
  
  year.LOS <- data.frame(matrix(NA, nrow=nrow(mig.stops[[n]]), ncol=2, dimnames=list(c(), c("year","LOS"))))
  
  for (i in 1:length(levels(dat$mgroup))){
    
    mgrouplevel <- levels(dat$mgroup)[i]
    usedat <- dat[dat$mgroup==mgrouplevel,]
    usedat <- droplevels(usedat)
    
    year <- levels(usedat$year)
    
    LOS <- sum(usedat$duration.days, na.rm=TRUE)
    
    year.LOS[i,] <- c(year, LOS)
        
  }
  
  fullroute[[n]] <- data.frame(mig.stops[[n]], year.LOS)
  fullroute[[n]]$LOS <- as.numeric(fullroute[[n]]$LOS)
  fullroute[[n]]$year <- as.factor(fullroute[[n]]$year)
  
}


###-----------------------------------------------------------###
#         6. ADD EXTRA VARIABLES ON MIG STRATEGY, LOS BIN
###-----------------------------------------------------------###

# re-concatenate individual bird list elements into a single dataframe with all individuals

allbirds <- do.call(rbind, fullroute) # rbind all list elements

allbirds$LOS <- as.numeric(allbirds$LOS)
allbirds$year <- as.factor(allbirds$year)

if (map2 == FALSE) {
  # remove stopovers < 3 days
  remove <- which(allbirds$LOS < 3 & allbirds$mgroup!=1)
  longstops <- allbirds[-remove,]
}

if (map2 == TRUE) {
  longstops <- allbirds
}

###--- assign migratory strategy (map 1 only) ---###

if (map2 == FALSE) {
  
  migstrategy <- function(name) {
    if (name=="BB" | name=="Chance" | name=="Chris" | name=="David" | name=="Indy" | name=="Iolo" | name=="Kasper" | name=="Livingstone" | name=="Lloyd" | name=="Martin" | name=="Mungo" | name=="Patch" | name=="Roy" | name=="Sussex" | name=="Sussex" | name=="Tor" | name=="Wallace" | name=="Waller") {
      strategy <- "SE"
    } else {
      strategy <- "SW"
    }
    return(strategy)
  }
  
  mig.strategy <- sapply(longstops$name, migstrategy)
  
}


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
  
    if (LOSvariable <= 14) {
      bin <- 0.3
    } else if (LOSvariable > 14 & LOSvariable <= 30) {
      bin <- 0.7
    } else {
      bin <- 1.1
    }
    return(bin)
  
}

LOS.bin <- sapply(longstops$LOS, groupLOS)

###--- FULL DATAFRAME ---###

if (map2 == FALSE) {
  cuckoo.pts <- data.frame(longstops, LOS.bin, mig.strategy)
}

if (map2 == TRUE) {
  cuckoo.pts <- data.frame(longstops, LOS.bin)
}

###--- convert long/lat to coordinates in WGS84 ---###

coordinates(cuckoo.pts) <- c("long.WGS84","lat.WGS84")
proj4string(cuckoo.pts) <- CRS("+proj=longlat +datum=WGS84")

###--- extract country / continent names ---###

geoinfo <- over(cuckoo.pts, EurAfrmap)[,c("REGION","ADMIN")]
colnames(geoinfo) <- c("continent", "country")
geoinfo <- droplevels(geoinfo)

toplot.pts <- data.frame(cuckoo.pts, geoinfo)

###--- convert long/lat to coordinates in WGS84 ---###

coordinates(toplot.pts) <- c("long.WGS84","lat.WGS84")
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

if (map2 == TRUE) {
  which(is.na(toplot.pts$country))
  
  toplot.pts[which(is.na(toplot.pts$continent)), "continent"] <- "Europe"
  toplot.pts[which(is.na(toplot.pts$country)), "country"] <- c("Netherlands")
  
}
  

###--- remove stopover points in UK that are breeding & > 10d ---###

remove <- which(toplot.pts$country=="United Kingdom")
toplot.pts <- toplot.pts[-remove,]


###-----------------------------------------------------------###
#         7. PLOT MAPS
###-----------------------------------------------------------###

##################### MAP 1 ########################

if (map2 == FALSE) {
  
  ###--- create a colour vector for migratory strategy ---###
  
  colstrategy <- factor(toplot.pts@data$mig.strategy, labels=c("blue","red"))
  colstrategy.lines <- factor(mig.routes@data$strategy, labels=c("blue", "red"))
  
  setwd(outputwd)
  
  tiff("all cuckoos southward migration and LOS_final.tiff",res=300,height=2200,width=2000,units="px")
  
  par(mar=c(0,0,0,0))
  plot(EurAfrmap, border="grey60")
  
  # plot major (> 3 day) stopovers
  plot(toplot.pts, pch=16, col=as.character(colstrategy), cex=toplot.pts$LOS.bin, add=T)
  
  # plot migration routes
  plot(mig.routes, col=as.character(colstrategy.lines), lwd=0.6, add=T)
  
  legend(35, 62, legend=c("< 14 days", "14-30 days", "> 30 days", "SW strategy", "SE strategy"), pch=16, cex=0.7, pt.cex=c(0.3, 0.7, 1.1, 0, 0), text.col=c("black","black","black","red","blue"), bty="n")
  
  dev.off()
  
}

##################### MAP 2 ########################

if (map2 == TRUE) {
  
  ###--- create a colour vector for each individual ---###
  
  colours <- c("blue", "deeppink", "forestgreen", "cyan3", "darkorchid1")
  
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
  
  legend(37, 65, legend=c("< 14 days", "14-30 days", "> 30 days", levels(toplot.pts@data$name)), pch=16, cex=0.7, pt.cex=c(0.3, 0.7, 1.1, 0, 0, 0, 0, 0), text.col=c("black", "black", "black", colours), bty="n")
  
  dev.off()
  
}


############## MAP 2 (version 2, separate maps/bird ########

if (map2 == TRUE) {
  
  ###--- create a colour vector for each individual ---###
  
  colours <- c("blue", "deeppink", "forestgreen", "cyan3", "darkorchid1")
  
  colstrategy <- factor(toplot.pts@data$name, labels=colours)
  colstrategy.lines <- factor(mig.routes@data$name, labels=colours)
  
  setwd(outputwd)
  
  tiff("multiyear cuckoos southward migration and LOS_v2.tiff",res=300,height=2000,width=2400,units="px")
  
  lo <- layout(rbind(c(1,2,3),c(4,5,6)))
  #layout.show(lo)
  par(mar=c(0,0,0,0))
  
  for (i in 1:length(levels(toplot.pts$name))){
    
    name <- levels(toplot.pts$name)[i]
    
    plot(EurAfrmap, border="grey60")
    plot(toplot.pts[toplot.pts$name==name,], pch=16, col=colours[i], cex=toplot.pts$LOS.bin, add=T)
    
    text(-5,-5, name, font=2, cex=1.5)
  }
  
  
  legend(80, 30, legend=c("< 14 days", "14-30 days", "> 30 days"), pch=16, cex=1.5, pt.cex=c(0.7, 1.1, 1.5), text.col=c("black", "black", "black"), xpd=NA, bty="n")
  
  dev.off()
  
    
}
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

