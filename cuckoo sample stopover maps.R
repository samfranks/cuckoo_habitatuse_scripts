##########################################################
#
#	Cuckoo sample stopover maps
#
# Sample stopover maps for cuckoo presentations
#
#	Samantha Franks
# 16 April 2014
# 10 Oct 2014 - more maps for Chris' Korea talk, cuckoo slides
#
##########################################################

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

#rm(list=ls())

### Load packages

library(adehabitat) # 
library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
#library(shapefiles) # allows you to read and write ESRI shapefiles
library(rgeos)
library(raster)
library(rasterVis)
library(geosphere)
library(adehabitatHR)
library(GISTools)

#### LOAD WORKSPACE ####
#load("~/Git/cuckoos/PA corine extraction.RData")
#load("~/cuckoos/PA corine temp precip extraction.RData") # load workspace on cluster
load("~/Git/cuckoos/PA epsg3035 corine temp precip extraction.RData") # load workspace on PC

#################################################
#
####      SET WORKING DIRECTORIES     ####
#
#################################################

parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output/sample maps", sep="")

resampwd <- ("/resampled data + extra grouping variables")
origwd <- ("/original data + extra grouping variables")

corinewd <- paste(parentwd, "/data/corine raster 100x100/v17", sep="")
elevationwd <- paste(parentwd, "/data/elevation data", sep="")


# if (!cluster) {
#   GISwd <- c("C:/Users/samf/Documents/GIS/cuckoos/protected areas")
# }
# 
# if (cluster) {
#   GISwd <- c("/users1/samf/cuckoos")
# }

corine.crs <- CRS("+init=epsg:3035")
epsg4326 <- CRS("+init=epsg:4326")

setwd(parentwd)

###-----------------------------------------------------------###
#         LOAD BASE LAYERS - CORINE & PROTECTED AREA RASTERS
###-----------------------------------------------------------###

setwd(corinewd)
r <- raster("g100_06.tif")

# Protected Area shapefile (PA.epsg3035) already loaded as part of the workspace, with correct projection for corine layer (epsg 3035)

###=====================================================================###
#====                  SET-UP DATA & STOPOVER POLYGONS                ====
###=====================================================================###

### SET SPECIFICATIONS FOR OUTPUT - which bird and which stopover (year and country); some birds may have multiple stops within a country in a given year
birdname <- "Chance" # name of cuckoo to map
whichstop <- "2013Germany" # yearcountry
original <- FALSE

whichyear <- as.numeric(substring(whichstop,1,4))  # 4 character year
whichcountry <- substring(whichstop, 5) # country name

###------------------ LOAD DATA, SUBSET, ADD PROJECTION INFO, TRANSFORM ------------------

if (original) setwd(paste(datawd,origwd,sep=""))
if (!original) setwd(paste(datawd,resampwd,sep=""))

dataset <- read.csv(paste(birdname, ".csv", sep=""), header=T)

# convert mgroup to factor
dataset$mgroup <- as.factor(dataset$mgroup)

dataset2 <- data.frame(dataset)

# subset dataset by stopover sites (no breeding sites)
dataset3 <- subset(dataset2, stopoversite=="Y")
dataset3 <- droplevels(dataset3)

# long/lat (original) are in CRS("+init=epsg:4326")
# long/lat (resampled) are in CRS("+init=epsg:3395")  
# newlongs/newlats (resampled) are in CRS("+init=epsg:4326")
# Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)


# ---------------------- RESAMPLED DATA -----------------------#
if (!original) {
  coordinates(dataset) <- c("newlongs","newlats")
  proj4string(dataset) <- CRS("+init=epsg:4326") 
  dataset.epsg3035 <- spTransform(dataset, CRS=corine.crs) 
}

# ---------------------- ORIGINAL DATA -----------------------#
if (original) {
  coordinates(dataset) <- c("long","lat")
  proj4string(dataset) <- CRS("+init=epsg:4326")
  dataset.epsg3035 <- spTransform(dataset, CRS=corine.crs)
  
}

#### ================= CREATE SUBSET OF DESIRED DATA FOR BIRD BY YEAR and COUNTRY =============== ####
newdat <- dataset.epsg3035[which(dataset.epsg3035$year==whichyear & dataset.epsg3035$country==whichcountry),]

# convert mgroup to factor
newdat@data$mgroup <- as.factor(newdat@data$mgroup)
newdat@data <- droplevels(newdat@data)

if (original) plotorig <- newdat
if (!original) plotresamp <- newdat


###------------------------------------------------------------###
#   CREATE DATASETS FOR ANALYSIS (include extra grouping variables) USING SOURCE CODE
###------------------------------------------------------------###


  
#################################################
#         MAPS
#################################################
  
# create new directory for specific bird to put output maps in, and setwd() as this new directory

gClip <- function(toclip, clipwith){
  clipwith <- extent(newdat)
  clipby <- as(extent(clipwith), "SpatialPolygons")
  gIntersection(toclip, clipby, byid = T)
}

setwd(outputwd)

mlarge <- 70000
msmall <- 20000

croppedraster.large <- crop(r,extent(plotresamp) + mlarge) # crop raster to +/- 10000m (10km) to each of N/S, E/W of stopover extent
croppedraster.small <- crop(r,extent(plotresamp) + msmall) # crop raster to +/- 10000m (10km) to each of N/S, E/W of stopover extent

croppedPA <- PA.epsg3035[plotresamp,]

tiff("sample map large scale.tiff", res=150,height=800,width=800,units="px")
### plot the large scale map showing entire PA coverage in area
plot(croppedraster.large)
plot(croppedPA, add=T, lwd=1.2, col=rgb(1,130,185,150,max=255))
#plot(croppedPA, lwd=1.2, col=rgb(1,130,185,150,max=255))
plot(plotresamp, cex=0.5, pch=16, col="black", add=T)
plot(plotorig, cex=0.5, pch=16, col="purple4", add=T)
# extentlarge.vector <- as.vector(extent(croppedraster.large))
# map.scale(extentlarge.vector[1]+10000, extentlarge.vector[3]+10000, 20000, "km", 5, sfcol="black")
dev.off()

PA.clipped <- gClip(croppedPA, extent(newdat)+msmall)

tiff("sample map small scale.tiff", res=150,height=800,width=800,units="px")
### plot the small scale map
plot(croppedraster.small)
plot(PA.clipped, add=T, lwd=1.2, col=rgb(1,130,185,150,max=255))
#plot(croppedPA, lwd=1.2, col=rgb(1,130,185,150,max=255))
plot(plotresamp, cex=1, pch=16, col="black", add=T)
plot(plotorig, cex=1, pch=16, col="purple4", add=T)
dev.off()

###-------------------------------------------------------###
#     EUROPE MAP, ALL STOPOVERS 
###-------------------------------------------------------###

# # export Europe autumn stopovers to a tiff file, with different colours for different years
# 
# if (length(levels(as.factor(newdat@data$year)))==1) {
#   colpoints <- factor(newdat@data$year, labels=c("black"))  
# } else if (length(levels(as.factor(newdat@data$year)))==2) {
#   colpoints <- factor(newdat@data$year, labels=c("darkmagenta","black"))  
# } else {
#   colpoints <- factor(newdat@data$year, labels=c("royalblue","darkmagenta","black"))  
# }
# 
# tiff(paste(levels(newdat$name)," - Europe.tiff"),res=300,height=2000,width=2400,units="px")
# 
# plot(Europeraster)
# plot(newdat, pch=16, cex=0.8, col=as.character(colpoints), add=T)
# coords <- cbind(2700000,2700000)
# textpts <- SpatialPoints(coords)
# proj4string(textpts) <- corine.crs
# text(textpts, levels(newdat$name))
# 
# dev.off()

###-----------------------------------------------###
#   CALCULATE MCP FOR RESAMPLED STOPOVER LOCATIONS
###-----------------------------------------------###

# for each movement group of each cuckoo, draw a MCP around the mgroup's points using mcp() to estimate home range
# choose points for MCP such that 5% of outliers are excluded
newdat@data <- droplevels(newdat@data)
stoppoints <- newdat[,c("mgroup")]
MCP <- mcp(stoppoints, percent=100)

### --- CROP Europeraster to box around stopover area --- ###

# clip the raster to a bounding box with an extent that buffers each MCP by ~ 5km and plot the clipped rasters with the stopover points on it, and export to a TIFF

# for each polygon in MCP, create a vector of the min(newlongs),max(newlongs),min(newlats),max(newlats), and use this to crop the raster Europeraster
# clip the raster to a bounding box with an extent that buffers the MCP by ~ 10km

MCPbox <- list()

# ??? can convert below loop to a function if desired

bufferrad <- 40000

for (i in 1:length(mcp)){
  polycoords <- mcp@polygons[[i]]@Polygons[[1]]@coords
  MCPbox[[i]] <- c(min(polycoords[,1])-bufferrad, max(polycoords[,1])+bufferrad, min(polycoords[,2])-bufferrad, max(polycoords[,2])+bufferrad)
  MCPbox
}

clipMCP <- function(MCPminmax){
  MCPbounds <- extent(MCPminmax)
  return(MCPbounds)
}

mcp.bbox <- lapply(MCPbox,clipMCP)

### --- EXTRACT STOPOVER MCP CENTROIDS --- ###

extractcentroid <- function(spdf){ #function to extract centroid values from the polygons (list components) of a SpatialPolygonsDataFrame
  centroid <- spdf@labpt
  return(centroid)
}

extractid <- function(spdf){ #function to extract id name (mgroup)
  id <- spdf@ID
  return(id)
}

centroids <- t(sapply(mcp@polygons,extractcentroid)) # represents centroids of mcps describing the stopovers
stopovers <- t(sapply(mcp@polygons,extractid)) # represents centroids of MCPs describing the stopovers

centroids <- as.data.frame(centroids)
stopovers <- as.factor(stopovers)

cent.stop <- data.frame(stopovers,centroids)
colnames(cent.stop) <- c("mgroup","centroidlong","centroidlat")

bird.name <- rep(levels(newdat$name), nrow(cent.stop))

bird.year <- rep(NA,length(mcp.bbox))
last.stopover <- factor(levels=c("Y","N"))

for (i in 1:length(mcp.bbox)){
  bird.year[i] <- newdat@data[which(newdat@data$mgroup==mcp$id[i])[1],"year"]
  last.stopover[i] <- newdat@data[which(newdat@data$mgroup==mcp$id[i])[1],"laststop"]
}

centroid.stopovers[[a]] <- data.frame(name=bird.name, year=bird.year, laststop=last.stopover, cent.stop)

# plot the new cropped rasters for each stopover movement, with the relevant stopover points and MCP on it
# export plots to a TIFF


for (i in 1:length(mcp.bbox)) {
  
  i<-1
  
  croppedraster <- crop(r,mcp.bbox[[i]])
  
  tiff(paste(levels(newdat$name), " - movement ", mcp$id[i], "raster showing Chance 29 random PAs.tiff", sep=""),res=150,height=800,width=800,units="px")
  
  plot(croppedraster)
  
  #plot(PA.LochLomond, add=T, lwd=1.2, col=rgb(1,130,185,100, max=255))
  #plot(PA.LochLomond, lwd=2, col=rgb(1,130,185,80,max=255))
  
  plot(PA.DEU.Chance, add=T, lwd=1.2, col=rgb(1,130,185,150,max=255))
  #plot(PA.DEU.Chance, lwd=1.2, col=rgb(1,130,185,100, max=255))
  
  polygontoplot <- SpatialPolygons(list(MCP@polygons[[i]])) # MCP around stopover points
  plot(polygontoplot, col=rgb(0,0,0,75,max=255), lwd = 2, add=T)
  
  polygontoplot <- SpatialPolygons(list(randomMCP@polygons[[i]])) # MCP around stopover points
  plot(polygontoplot, col=rgb(255,255,255,200,max=255), lwd = 2, add=T)
  
  #plot(mergePA, add=T, border="blue", col= "blue", ID="29") 
  
  plot(newdat[which(newdat$mgroup==mcp$id[i]),], pch=16, cex=0.4, col="black", add=T)
  
  coords <- cbind(MCPbox[[i]][1]+500,MCPbox[[i]][4]-3000)
  textpts <- SpatialPoints(coords)
  proj4string(textpts) <- corine.crs
  
#   plottext <- paste(levels(newdat$name), " - #", MCP$id[i], " - ", newdat@data[which(newdat@data$mgroup==MCP$id[i])[1], "country"], " ", newdat@data[which(newdat@data$mgroup==MCP$id[i])[1], "year"], sep="")
#   text(textpts, plottext, cex=0.8, pos=4)
#   
#   plot circles showing radiuses within which new random MCP is generated
#   
#   randomradius <- c(50000,100000,200000,500000)
#   circleMCP <- list()
#   
#   for(j in 1:length(randomradius)) {
#   
#     circleMCP[[j]] <- getellipse(rx=randomradius[j], ry=randomradius[j], mid=c(4521624,3223129))
#   
#   }
#   
#   bufferpolygon <- SpatialPolygons(list(Polygons(list(Polygon(circleMCP[[1]]), Polygon(circleMCP[[2]]), Polygon(circleMCP[[3]]), Polygon(circleMCP[[4]])), ID="a")))
#   plot(bufferpolygon, add=T, lwd=1.5)
  
  
  # plot original data points (darkmagenta) from transmission cycle over resampled points (black) from transmission cycle for each stopover
  #plot(newdat.orig[which(newdat.orig$mgroup==mcp$id[i]),], pch=16, cex=0.4, col="darkmagenta", add=T)
  
  dev.off()
  
}

DEU.Chance29 <- PAoverlap.info[[a]][PAoverlap.info[[a]]$mgroup==29, "wdpaid"]

#PA.LochLomond <- subset(PA.subset, wdpaid==183407)
PA.DEU.Chance <- PA.subset[which(PA.subset@data$wdpaid %in% DEU.Chance29),]



# newdat.orig@data[which(newdat.orig@data$mgroup=="26"),] # non-resampled original point id=53488
# 
# newdat.orig@coords[209,]
# 
# newdat@data[newdat@data$id==53488,]


