##########################################################
#
#  Cuckoo resampled AND original location point land cover value extractions
#  Plots of proportional habitat use
#
#	Samantha Franks
#	28 Nov 2013
# 17 Dec 2013 - a) run with corrected mgroups + mtypes from rerun resampling data; 2) now also includes code for to run originl data
# 24 Dec 2013 - now removes mgroups with points from only a single tcycle
#
#
##########################################################


library(adehabitat) # 
library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
#library(shapefiles) # allows you to read and write ESRI shapefiles
library(rgeos)
library(raster)
library(geosphere)
library(adehabitatHR)
library(rworldmap)
library(rworldxtra)
library(NCStats)

if (!("Europeraster" %in% ls())) { ### run prep code to load rasters if not in workspace already
  
  ### -------- Working directories --------- ###
  
  # GISwd <- c("C:/Users/samf/Documents/R/projects/cuckoos/ArcGIS/all original layers WGS84 export")
  corinewd <- c("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/data/corine raster 100x100")
  
  ###-----------------------------------------------------------###
  #         Load base layers
  ###-----------------------------------------------------------###
  
  ### CORINE layer ###
  
  # projection and datum information for all maps
  # epsg: 3395 - +proj=merc +ellps=wGS84 +datum=WGS84
  # epsg: 4326 - +proj=latlong +ellps=wGS84 +datum=WGS84
  ### epsg: 3035 - LAEA CRS, used for the Corine Land Cover 2006 data
  
  corine.crs <- CRS("+init=epsg:3035")
  
  ### ------- CORINE land cover raster file ------- ###
  setwd(corinewd)
  r <- raster("g100_06.tif")
  
  Europebox <- extent(2483500,5890375,1276625,4286750) # bounding box extent for relevant parts of Europe (no Scandanavia or Turkey)
  Europeraster <- crop(r,Europebox)
  
  ### ------- Europe/Africa national boundaries layer ------- ###
  
  dat.world <- getMap(resolution="high")
  
  Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
  Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
  Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]
  
  Eur.Afr <- spTransform(Eur.Afr,CRS("+proj=longlat +datum=WGS84")) # transform map to exactly the same projection as cuckoo newlongs/newlats data
  
  Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
  Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]
  
  ### -------- Protected area network -------- ###
  
}

###------------------------------------------------------------###
#         EXTRACT CORINE VALUES FOR LOCATION POINTS
###------------------------------------------------------------###

### TO RUN ORIGINAL or RESAMPLED data, change below line

original <- TRUE # or FALSE if resampled

# create list to hold datasets with extracted corine values for each cuckoo

corine.values <- list()

### --- Set working directories --- ###

setwd("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos")
parentwd <- getwd()

datawd <- ("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/data")
outputwd <- ("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/output")

resampwd <- ("/resampled data with uncertainty 100 bootstraps")
origwd <- ("/original data + distance, movement groups, etc")

# for loop for each cuckoo file starts here
for (a in 1:31) {
  
  ### --- RESAMPLED or ORIGINAL cuckoo data: Load and add geoinfo --- ###
  
  # change to resampwd or origwd depending on whether corine values are being extracted for resampled or original data
  
  if (original) setwd(paste(datawd,origwd,sep=""))
  if (!original) setwd(paste(datawd,resampwd,sep=""))
  cuckoofiles <- list.files()
  
  cuckoo <- read.csv(cuckoofiles[a], header=T)
  
  # long/lat are in CRS("+proj=longlat +datum=WGS84")
  # newlongs/newlats are in CRS("+init=epsg:3395")
  # centroidlong/lat are in CRS("+init=epsg:3395")
  
  # ---------------------- RESAMPLED DATA -----------------------#
  # transform newlongs/newlats from World Mercator/WGS84 to longlat WGS84
  if (!original) {
    coordinates(cuckoo) <- c("newlongs","newlats")
    proj4string(cuckoo) <- CRS("+init=epsg:3395") 
    cuckoo <- spTransform(cuckoo,CRS("+proj=longlat +datum=WGS84")) 
  }
  
  # ---------------------- ORIGINAL DATA -----------------------#
  if (original) {
    coordinates(cuckoo) <- c("long","lat")
    proj4string(cuckoo) <- CRS("+proj=longlat +datum=WGS84")
  }
  
  ### over() retrieves values of Eur.Afr at points described by cuckoo
  # extracts continent and administrative country names from Eur.Afr
  geoinfo <- over(cuckoo, Eur.Afr)[,c("REGION","ADMIN")]
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
  
  ### add geographic info to cuckoo dataset and redefine projection of newlongs/newlats
  fullcuckoo <- data.frame(cuckoo, geoinfo)
  
  # ---------------------- RESAMPLED DATA -----------------------#
  if (!original) {
    coordinates(fullcuckoo) <- c("newlongs","newlats")
    proj4string(fullcuckoo) <- CRS("+proj=longlat +datum=WGS84")
 }
  
  # ---------------------- ORIGINAL DATA -----------------------#
  if (original) {
    coordinates(fullcuckoo) <- c("long","lat")
    proj4string(fullcuckoo) <- CRS("+proj=longlat +datum=WGS84")
  }
  
  
  
  ### create Europe only (+ at sea points) and Africa only data sets
  Eurcuckoo <- fullcuckoo[which(fullcuckoo$continent=="Europe" | fullcuckoo$continent=="at sea"),]
  Afrcuckoo <- fullcuckoo[which(fullcuckoo$continent=="Africa"),]
  
  ###------------------------------------------------------------###
  #         Trim European data - autumn stopovers only
  ###------------------------------------------------------------###
  
  # subset out only stopovers in Europe on southward migration
  
  ### RESAMPLED or ORIGINAL CUCKOO DATA
  
  Eurcuckoo$mgroup <- as.factor(Eurcuckoo$mgroup)
  
  #stopovers <- Eurcuckoo[which(Eurcuckoo$mtype=="C" | Eurcuckoo$mtype=="S" | is.na(Eurcuckoo$mtype)),] # includes capture (breeding) location, stopovers, as well as last known location (if in Europe)
  #stopovers <- stopovers[which(stopovers$mgroup!=1),] # this line removes the capture occasion
  
  autumn <- Eurcuckoo[which(Eurcuckoo$julian >= 131),] # all julian dates after May 10 (removes spring stopovers)
  autumn@data <- droplevels(autumn@data)
  
  autumn$mgroup <- as.factor(autumn$mgroup)
  
  # can add other subsetting variables here, and change whatever "dat" is called
  
  # for each level of mgroup, check the number of levels of tcycle - if < 2, then Subset data to remove that mgroup
  
  newautumn <- autumn
  
  if (original) {
    # if original data being used, then Subset mtype to be C,S, or N/A
    newautumn <- subset(newautumn, newautumn$mtype=="C" | newautumn$mtype=="S" | is.na(newautumn$mtype))
  }
  
  for (i in 1:length(levels(autumn$mgroup))) {
    temp <- newautumn[newautumn$mgroup==levels(autumn$mgroup)[i],]
    temp@data <- droplevels(temp@data)
    if (length(levels(as.factor(temp$tcycle))) < 2) {
      newautumn <- newautumn[which(newautumn$mgroup!=levels(autumn$mgroup)[i]),]
      newautumn@data <- droplevels(newautumn@data)
    } else {
      newautumn <- newautumn
    }
  }
  
  dat <- newautumn
  
  # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
  
  newdat <- spTransform(dat, CRS=corine.crs)
  
  ###-------------------------------------------------------------###
  #   Land-cover point value extractions
  ###-------------------------------------------------------------###
  
  corine.values[[a]] <- data.frame(newdat, corine.values=extract(Europeraster, newdat@coords))
    
}

# set directory of csv output dataset with corine values
if (!original) setwd(paste(datawd,"/resampled data corine extracted values",sep=""))
if (original) setwd(paste(datawd,"/original data corine extracted values",sep=""))


for (i in 1:31){
   write.csv(corine.values[[i]], paste(corine.values[[i]]$name[1],".csv", sep=""), row.names=FALSE)
}

######## ADD SOURCE CODE HERE: "add variables to corine values. . . 

setwd(paste(parentwd,"/scripts",sep=""))
source("add variables to corine values - age, breeding mgroup, migration strategy, and Sahara crossing.R")

# write output of source file - elements of list corinewithgroupvar - each to .csv file
# CHANGE folder to ORIGINAL or RESAMPLED data folder in the parent script

if (!original) setwd(paste(datawd,"/resampled data corine + extra grouping variables",sep=""))
if (original) setwd(paste(datawd,"/original data corine + extra grouping variables",sep=""))

for (i in 1:31){
  write.csv(corinewithgroupvar[[i]],(paste(levels(corinewithgroupvar[[i]]$name),".csv",sep="")), row.names=F)
}

################################

###---------------------------------------------------------###
#
#  Barplots of land-cover types used
#
###---------------------------------------------------------###

### --- FUNCTION TO PLOT CUCKOO HABITAT USE, MULTIPLE SITES POOLED --- ###

habsummary <- function(dataset,landcover){ # datasets required are point land cover values, clclegend data (landcover), and par settings for specific # of stopovers
  
  dataset <- droplevels(dataset)
  dataset$mgroup <- as.factor(as.character(dataset$mgroup))
  
  # calculate proportion of each habitat type used
  prophab <- prop.table(table(dataset$corine.name))
  
  # produce colour labels for barplot according to RGB colour codes of CLC legend
  landcolours <- rgb(landcover[which(landcover$LABEL4 %in% levels(dataset$corine.name)),c("rgb1","rgb2","rgb3")], max=255)
  
  # plot the proportional habitat use
  
  par(mar=c(9,4,2,1), oma=c(0,0,2,0)) # specific par settings for number of plots to produce
  
  x <- barplot(prophab, names.arg=c(""), ylim=c(0,max(prophab)+0.05), col=landcolours, cex.axis=0.8, las=2)
  text(x, par("usr")[3] - 0.01, srt = 45, adj = 1, xpd = TRUE, cex=0.7, labels = levels(dataset$corine.name))
  
  #   label for breeding site or stopovers based on whether breedsite=Y or N
  if (levels(dataset$breedsite)=="Y") {
    mtext(paste(dataset$name[1], ", breeding habitat use in ", dataset$country[1],sep=""), cex=0.8,font=2,line=1)
  } else {
    mtext(paste(dataset$name[1], ", overall stopover habitat use in Europe", sep=""), cex=0.8,font=2,line=1)
  }
  
}

### --- FUNCTION TO PLOT CUCKOO HABITAT USE, SINGLE SITES, MULTIPLE PLOTS PER PAGE --- ###

habsummary.stopovers <- function(dataset,landcover){ # datasets required are point land cover values, clclegend data (landcover), and par settings for specific # of stopovers
  
  dataset <- droplevels(dataset)
  dataset$mgroup <- as.factor(as.character(dataset$mgroup))
  
  # calculate proportion of each habitat type used
  prophab <- prop.table(table(dataset$corine.name))
  
  # produce colour labels for barplot according to RGB colour codes of CLC legend
  landcolours <- rgb(landcover[which(landcover$LABEL4 %in% levels(dataset$corine.name)),c("rgb1","rgb2","rgb3")], max=255)
  
  # plot the proportional habitat use
  
  x <- barplot(prophab, names.arg=c(""), ylim=c(0,max(prophab)+0.05), col=landcolours, cex.axis=0.8, las=2)
  text(x, par("usr")[3] - 0.03, srt = 45, adj = 1, xpd = TRUE, cex=0.7, labels = levels(dataset$corine.name))
  
  # if dataset has a single stopover, then label is for that specific stopover; otherwise, label is for Europe generally
  mtext(paste(dataset$name[1], ", habitat use in ", dataset$country[1], " ",dataset$year[1],"\nLast stopover? = ", dataset$laststop[1], sep=""), cex=0.5,font=2,line=1)
  
}

### --- SET-UP DATA --- ###

### import corine legend
setwd(datawd)
clclegend <- read.csv("clc_legend.csv", header=T)

### Set RESAMPLED or ORIGINAL data
if (!original) dataid <- c("resampled data")
if (original) dataid <- c("original data")

### Begin loop through bird datasets

# leave out Idemili ([8])
for (a in 1:30) {
  ### select cuckoo dataset
  setwd(paste(datawd,"/",dataid," corine + extra grouping variables",sep=""))
  
  corinefiles <- list.files()[-8]
  dat <- read.csv(corinefiles[a], header=T)
  
  ### identify last stopovers in Europe
  # identify last stopovers before Sahara crossing
  # for each level of year, identify max(mgroup)
  
  # convert mgroup to integer, via character
  dat$mgroup <- as.integer(as.character(dat$mgroup))
  
  yearlevels <- levels(as.factor(dat$year))
  maxmgroup.byyear <- rep(NA, length(yearlevels))
  
  for (i in 1:length(yearlevels)){
    maxmgroup.byyear[i] <- max(dat[dat$year==yearlevels[i], "mgroup"]) # extract max mgroup number for each level of year
  }
  
  laststop <- factor(levels=c("N","Y"))
  laststop[which(dat$mgroup %in% maxmgroup.byyear)] <- "Y"
  laststop[which(!dat$mgroup %in% maxmgroup.byyear)] <- "N"
  dat <- data.frame(dat,laststop)
  
  dat <- Subset(dat, !is.na(corine.values))
  
  ### add names of corine land types to cuckoo dataset
  corine.name <- rep(NA,nrow(dat))
  corine.nameLABEL1 <- rep(NA,nrow(dat))
  corine.nameLABEL2 <- rep(NA,nrow(dat))
  
  for (i in 1:nrow(dat)){
    corine.name[i] <- as.character(clclegend[which(clclegend$GRID_CODE %in% dat[i,"corine.values"]), "LABEL4"])
    corine.nameLABEL1[i] <- as.character(clclegend[which(clclegend$GRID_CODE %in% dat[i,"corine.values"]), "LABEL1"])
    corine.nameLABEL2[i] <- as.character(clclegend[which(clclegend$GRID_CODE %in% dat[i,"corine.values"]), "LABEL2"])
  }
  
  newland <- data.frame(dat, corine.name, corine.nameLABEL1, corine.nameLABEL2)
  
  # order land cover levels in dataset so they appear in order of grid codes in barplot
  orderedcorine <- clclegend$LABEL4[which(clclegend$LABEL4 %in% levels(newland$corine.name))] 
  newland$corine.name <- factor(newland$corine.name, orderedcorine, ordered=TRUE)
  
  write.csv(newland,paste(levels(newland$name),".csv",sep=""), row.names=F)

}      ### REMOVE this bracket if you want to output barplots of habitat use. Leave bracket if you only want to produce datasets.
  
  ### --- LAND-USE DATASETS & PLOTS --- ###
  
  setwd(paste(outputwd,"/Europe habitat use/", dataid,sep=""))
  
  # UK breeding site
  breeding <- Subset(newland, breedsite=="Y")
  
  tiff(paste(levels(breeding$name)," - UK breeding - ", dataid, ".tiff ", sep=""),res=300,height=2000,width=2400,units="px")
  
  habsummary(breeding,clclegend)
  
  dev.off()
  
  # all sites in Europe, stopovers only, pooled - no breeding site or other UK movements that are not an autumn stopover
  stopovers <- Subset(newland, breedsite == "N")
  
  if (nrow(stopovers)>0) {
  
  tiff(paste(levels(stopovers$name)," - all stopovers pooled - ", dataid, ".tiff", sep=""),res=300,height=2000,width=2400,units="px")
  
  habsummary(stopovers,clclegend)
  
  dev.off()
  
  
  
  # all stopovers, individually
  # convert mgroup to level and split data by stopover (level of mgroup)
  stopovers$mgroup <- as.factor(as.character(stopovers$mgroup))
  stopovers$mgroup <- reorder(stopovers$mgroup)
  bymgroup <- stopovers$mgroup
  splitstopovers <- split(stopovers, list(bymgroup))
  
  
  tiff(paste(levels(newland$name)," - individual stopovers - ", dataid, ".tiff", sep=""),res=300,height=2000,width=3000,units="px")
  
  if (length(splitstopovers) < 8){
    par(mfrow=c(2, ceiling(length(splitstopovers)/2)), mar=c(9,4,3,2), oma=c(0,0,2,1))
  } else {
    par(mfrow=c(3, ceiling(length(splitstopovers)/3)), mar=c(9,4,3,2), oma=c(0,0,2,1))
    
  }
  
  for (i in 1:length(splitstopovers)){
    
    habsummary.stopovers(splitstopovers[[i]],clclegend)
  }
  
  dev.off()
  
  }
  
}
