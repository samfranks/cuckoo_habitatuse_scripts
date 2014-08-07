##########################################################
#
# Source code to extract cuckoo locations from different times of year
#
# Samantha Franks
# 7 Mar 2014
#
##########################################################

###--- SUMMARY ---###

# this code takes cuckoo location points for a single dataset (single .csv file) and subsets the data according to desired times of year and season
# ALSO
# - identifies final stopovers made during southward migration in Europe
# - adds migration strategy
# - adds age
# - adds success of Sahara crossing
# can be used with either original or resampled location points

###--- PARENT DATA REQUIRED ---###

### data are fed into the source from the parent script

# original <- TRUE or FALSE - is original (non-resampled) data being used?
# dataset <- x, where x is the name of the dataset you want to extract info from
# jdate1 <- a, where a is the minimum julian date bracketing your period of interest
# jdate2 <- b, where b is the maximum julian date bracketing your period of interest
# countryname <- y, where y is the name of the country/countries you want to extract info from
# continentname <- z, where z is the name of the continent you want to extract info from
# select.mtype <- m, where m is the type of movement to use or not use (M=migration, C=capture, S=stopover)

######### NOTES ABOUT THE DATASETS #########

###--- Projection and datum information ---###
# resampled cuckoo locations - epsg: 3395 - +proj=merc +ellps=wGS84 +datum=WGS84
# original cuckoo locations - epsg: 4326 - +proj=longlat +ellps=wGS84 +datum=WGS84
# CORINE data - epsg: 3035 - LAEA CRS

###--- README ---###

# name = individual id
# id = id value of location point
# tcycle = value assigned to a transmission cycle
# timestamp = UNIX translated timestamp for a point
# location_date = original UNIX timestamp
# julian = julian date of point
# year, month, day, hour = time information about point
#   - includes only points between May 10 (julian 131) up until the date a bird left Europe
#   - subsetted for autumn migration points only in "cuckoo migration, stopovers and land cover.R"
# loc.class = location class (quality) of the point: 3 best, 1 worst
# error.major/minor,orientation = information about the error ellipse - semi-major/minor axes + orientation
# sensor data
# test_location = whether the point was a test location or not (a reference point in North America or outside Chris Hewson's house in Cambridge - should all be N by this point, having removed any Y from the dataset)
# distance, bearing = distance and bearing moved between original location points (calculated from sequential points in the original dataset, and propagated with the observation through the bootstrapping procedure)
#   - "cuckoo original data tidying source.R"
# mgroup, mtype = movement group and type, determined from distance/bearings of the original data - at this point in the analysis, should not include any M (migration) movements, nor mgroups with < 2 transmission cycles (subsetted out in cuckoo migration, stopovers and land cover maps.R), also should be only points in Europe
# newlongs/newlats (resampled) = epsg: 3395 +proj=merc +ellps=wGS84 +datum=WGS84
# long/lat (original) = epsg: 4326 - +proj=longlat +ellps=wGS84 +datum=WGS84

###--- LOAD PACKAGES ---###

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

##################################################
##################################################

###-----------------------------------------------------------###
#     PARENT INFORMATION ON VARIABLES (passed from PARENT SCRIPT)
###-----------------------------------------------------------###

### NOTE: if using as source code, then can leave out LOAD DATASET section and use the data passed from parent script.  Loop should be set up in the parent script to run through all datasets in a folder

original <- FALSE # whether original or resampled data
# dataset <- x, where x is the name of the dataset you want to extract info from

# Clip dataset to only southward migration breeding May 10 to Dec 21 (julian >=131 and <=365)
jdate1 <- 131 #minimum julian date bracketing period of interest
jdate2 <- 365 #maximum julian date bracketing period of interest

continentname <- "Europe" #name of the continent you want to extract info from

# countryname <- y, where y is the name of the country/countries you want to extract info from

select.mtype <- "M"

##################################################
##################################################


###-----------------------------------------------------------###
#         SET WORKING DIRECTORIES
###-----------------------------------------------------------###

parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")

datawd <- ("C:/Users/samf/Documents/Git/cuckoos/data")
outputwd <- ("C:/Users/samf/Documents/Git/cuckoos/output")

resampwd <- ("/resampled data with uncertainty 100 bootstraps")
origwd <- ("/original data + distance, movement groups, etc")

orig.data.outputwd <- ("/original data + extra grouping variables")
resamp.data.outputwd <- ("/resampled data + extra grouping variables")

##################################################
##################################################

###-----------------------------------------------------------###
#         LOAD DATASET - left out if used as source code
###-----------------------------------------------------------###

# change to resampwd or origwd depending on whether corine values are being extracted for resampled or original data

if (original) setwd(paste(datawd,origwd,sep=""))
if (!original) setwd(paste(datawd,resampwd,sep=""))

cuckoofiles <- list.files()

# a <- 3 #test

# if this is being used as source code, then "dataset" should be the name of the dataset passed to this code
dataset <- read.csv(cuckoofiles[a], header=T)

##################################################
##################################################


###-----------------------------------------------------------###
#         LOAD MAP BASE LAYERS
###-----------------------------------------------------------###

### ------- Europe/Africa national boundaries layer ------- ###

dat.world <- getMap(resolution="high")

Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]

Eur.Afr <- spTransform(Eur.Afr,CRS("+init=epsg:4326")) # transform map to exactly the same projection as cuckoo newlongs/newlats data

Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]

###-----------------------------------------------------------###
#         EXTRACT GEOINFO AND ADD TO DATASET
###-----------------------------------------------------------###

# long/lat (original) are in CRS("+proj=longlat +datum=WGS84"), or CRS("+init=epsg:4326")
# long/lat (resampledl) are in CRS("+init=epsg:3395")
# newlongs/newlats are in CRS("+init=epsg:3395")
# centroidlong/lat are in CRS("+init=epsg:3395")

if (original) {
  coordinates(dataset) <- c("long","lat")
  proj4string(dataset) <- CRS("+init=epsg:4326")
}

if (!original) {
  coordinates(dataset) <- c("newlongs","newlats")
  proj4string(dataset) <- CRS("+init=epsg:3395")
  dataset <- spTransform(dataset,CRS("+init=epsg:4326")) # transform newlongs/newlats from World Mercator/WGS84 to longlat WGS84
}

### over() retrieves values of Eur.Afr at points described by the dataset
# extracts continent and administrative country names from Eur.Afr
geoinfo <- over(dataset, Eur.Afr)[,c("REGION","ADMIN")]
colnames(geoinfo) <- c("continent", "country")
geoinfo <- droplevels(geoinfo)

### add a new factor for continent/country NAs, which should be points at sea
atsea.pts <- which(is.na(geoinfo$continent))

continent.factors <- c(levels(geoinfo$continent), "at sea")
levels(geoinfo$continent) <- continent.factors
geoinfo$continent[atsea.pts] <- "at sea"

country.factors <- c(levels(geoinfo$country),"at sea")
levels(geoinfo$country) <- country.factors
geoinfo$country[atsea.pts] <- "at sea"

### add geographic info to dataset and redefine projection

if (original) {
  fulldataset <- data.frame(dataset, geoinfo)
  coordinates(fulldataset) <- c("long","lat")
  proj4string(fulldataset) <- CRS("+init=epsg:4326")
}

if (!original) {
  fulldataset <- data.frame(dataset, geoinfo)
  coordinates(fulldataset) <- c("newlongs","newlats")
  proj4string(fulldataset) <- CRS("+init=epsg:4326")
}

### create continent- or country-specific geoinfo datasets

# if Europe, then also include at sea points

if (continentname == "Europe") {
  geodataset <- subset(fulldataset, continent==continentname | continent=="at sea")
  geodataset@data <- droplevels(geodataset@data)
}

if (continentname != "Europe") {
  geodataset <- subset(fulldataset, continent==continentname)
  geodataset@data <- droplevels(geodataset@data)
}

if ("countryname" %in% ls()) {
  geodataset <- subset(fulldataset, country==countryname)
}
  
###-----------------------------------------------------------###
#         SUBSET DATA OF INTEREST (DATES, MTYPES, STOPOVERS WITH > 2 TRANSMISSION CYCLES)
###-----------------------------------------------------------###

# SUBSET 1 - SEASON (DATES)
subset1 <- subset(geodataset, julian >= jdate1 & julian <= jdate2)
subset1@data <- droplevels(subset1@data)

# SUBSET 2 - MOVEMENT TYPE (STOPOVERS ONLY)
subset2 <- subset(subset1, mtype!=select.mtype)
subset2@data <- droplevels(subset2@data)

# SUBSET 3 - MGROUPS WITH # of TCYCLES < 2
subset3 <- subset2
#subset3$mgroup <- as.factor(subset3$mgroup)

for (i in 1:length(levels(as.factor(subset2$mgroup)))) {
  
  mgroup.number <- as.integer(levels(as.factor(subset2$mgroup))[i])
  
  temp <- subset(subset2, mgroup==mgroup.number)
  temp@data <- droplevels(temp@data)
  
  if (length(levels(as.factor(temp$tcycle))) < 2) {
    subset3 <- subset(subset3, mgroup!=mgroup.number)
    subset3@data <- droplevels(subset3@data)
  } else {
    subset3 <- subset3
  }
}

### ADD MORE SUBSETTING VARIABLES IF DESIRED

nextdat <- subset3

###-----------------------------------------------------------###
#         ADD LAST EUROPEAN STOPOVER VARIABLE
###-----------------------------------------------------------###

yearlevels <- levels(as.factor(nextdat$year))
maxmgroup.byyear <- rep(NA, length(yearlevels))

for (i in 1:length(yearlevels)){
  maxmgroup.byyear[i] <- max(nextdat@data[nextdat@data$year==yearlevels[i], "mgroup"]) # extract max mgroup number for each level of year
}

laststop <- factor(levels=c("N","Y"))
laststop[which(as.factor(nextdat@data$mgroup) %in% maxmgroup.byyear)] <- "Y"
laststop[which(!as.factor(nextdat@data$mgroup) %in% maxmgroup.byyear)] <- "N"
with.laststop <- data.frame(nextdat,laststop)

###-----------------------------------------------------------###
#         ADD MIGRATION STRATEGY, AGE, and SAHARA SUCCESS
###-----------------------------------------------------------###

setwd(datawd)
othervar <- read.csv("cuckoo migratory strategy and Sahara crossing success.csv", header=TRUE)

strategy <- rep(othervar[othervar$name==levels(with.laststop$name),"migratory.strategy"], nrow(with.laststop))

age <- rep(othervar[othervar$name==levels(with.laststop$name),"age.at.capture"], nrow(with.laststop))

Sahara.success <- factor(levels=c("N","Y"))

if (is.na(othervar[othervar$name==levels(with.laststop$name),"year.failed.Sahara.crossing"]))
  Sahara.success[1:nrow(with.laststop)] <- "Y"

if(!(is.na(othervar[othervar$name==levels(with.laststop$name),"year.failed.Sahara.crossing"]))) {
  Sahara.success[which((with.laststop$year) %in% othervar[othervar$name==levels(with.laststop$name),"year.failed.Sahara.crossing"])] <- "N"
  Sahara.success[which(!(with.laststop$year) %in% othervar[othervar$name==levels(with.laststop$name),"year.failed.Sahara.crossing"])] <- "Y"
}

stopoversite <- factor(levels=c("N","Y"))
stopoversite[1:length(with.laststop$mgroup)] <- "Y"

# if only one "breeding site"
if (is.na(othervar[othervar$name==levels(with.laststop$name),"breed.mgroup2"])) {
  stopoversite[which((with.laststop$mgroup) %in% othervar[othervar$name==levels(with.laststop$name), "breed.mgroup1"])] <- "N"
  #stopoversite[which(!(with.laststop$mgroup) %in% othervar[othervar$name==levels(with.laststop$name),"breed.mgroup1"])] <- "Y"
} 

# if 2 "breeding sites"
if (is.na(othervar[othervar$name==levels(with.laststop$name),"breed.mgroup3"])) {
  stopoversite[which((with.laststop$mgroup) %in% othervar[othervar$name==levels(with.laststop$name), "breed.mgroup1"])] <- "N"
  stopoversite[which((with.laststop$mgroup) %in% othervar[othervar$name==levels(with.laststop$name), "breed.mgroup2"])] <- "N"
}

# if 3 "breeding sites"
if (!(is.na(othervar[othervar$name==levels(with.laststop$name),"breed.mgroup3"]))) {
  stopoversite[which((with.laststop$mgroup) %in% othervar[othervar$name==levels(with.laststop$name), "breed.mgroup1"])] <- "N"
  stopoversite[which((with.laststop$mgroup) %in% othervar[othervar$name==levels(with.laststop$name), "breed.mgroup2"])] <- "N"
  stopoversite[which((with.laststop$mgroup) %in% othervar[othervar$name==levels(with.laststop$name), "breed.mgroup3"])] <- "N"
}

newdataset <- data.frame(with.laststop,stopoversite,strategy,age,Sahara.success)

###-----------------------------------------------------------###
#         OUTPUT NEW DATASETS
###-----------------------------------------------------------###

### outputs are without projection information
# original data (longs, lats) are in epsg: 4326
# resampled data (newlongs/newlats) are in epsg: 4326
# resampled data (longs/lats) are in epsg: 3395

if (original) setwd(paste(datawd,orig.data.outputwd,sep=""))
if (!original) setwd(paste(datawd,resamp.data.outputwd,sep=""))

write.csv(newdataset,(paste(levels(newdataset$name),".csv",sep="")), row.names=F)