##########################################################
#
#	Cuckoo sample stopover maps
#
# Sample stopover maps for cuckoo presentation BTO mini research day
#
#	Samantha Franks
# 16 April 2014
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
library(rworldmap) # getMap function
library(rworldxtra)

###-----------------------------------------------------------###
#         SET WORKING DIRECTORIES
###-----------------------------------------------------------###

parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

resampwd <- ("/resampled data with uncertainty 100 bootstraps")
origwd <- ("/original data + distance, movement groups, etc")

origwd.extravar <- ("/original data + extra grouping variables")
resampwd.extravar<- ("/resampled data + extra grouping variables")

corinewd <- paste(parentwd, "/data/corine raster 100x100", sep="")

corine.crs <- CRS("+init=epsg:3035")

setwd(parentwd)

###-----------------------------------------------------------###
#         LOAD BASE LAYERS - CORINE & PROTECTED AREA RASTERS
###-----------------------------------------------------------###

setwd(corinewd)

corine.crs <- CRS("+init=epsg:3035")

### ------- CORINE land cover raster file ------- ###
r <- raster("g100_06.tif")

Europebox <- extent(2483500,5890375,1276625,4286750) # bounding box extent for relevant parts of Europe (no Scandanavia or Turkey)
#Europeraster <- crop(r,Europebox)

### ------- Europe/Africa national boundaries layer ------- ###

dat.world <- getMap(resolution="high")

Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]

Eur.Afr <- spTransform(Eur.Afr,CRS("+init=epsg:4326")) # transform map to exactly the same projection as cuckoo newlongs/newlats data

Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]

###------------------------------------------------------------###
#   CREATE DATASETS FOR ANALYSIS (include extra grouping variables) USING SOURCE CODE
###------------------------------------------------------------###

# create list to hold centroids of stopover points for all cuckoos
centroid.stopovers <- list()

a <- 2

setwd(paste(datawd,origwd.extravar,sep=""))
originaldata <- read.csv(list.files()[a], header=TRUE)
coordinates(originaldata) <- c("long","lat")
proj4string(originaldata) <- CRS("+init=epsg:4326")
newdat.orig <- spTransform(originaldata, CRS=corine.crs)


setwd(paste(datawd,resampwd.extravar,sep=""))
resampleddata <- read.csv(list.files()[a], header=TRUE)
coordinates(resampleddata) <- c("newlongs","newlats")
proj4string(resampleddata) <- CRS("+init=epsg:4326")
newdat <- spTransform(resampleddata, CRS=corine.crs)


  
#################################################
#         MAPS
#################################################
  
# create new directory for specific bird to put output maps in, and setwd() as this new directory

mapoutputdir <- paste("C:/Users/samf/Documents/Git/cuckoos/output/sample maps", sep="")

setwd(mapoutputdir)

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


