##########################################################
#
#	Cuckoo stopover locations, land cover and protected areas
#
#	Samantha Franks
#	7 Nov 2013
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
library(rworldmap)
library(rworldxtra)
library(NCStats)

### -------- Working directories --------- ###

# GISwd <- c("C:/Users/samf/Documents/R/projects/cuckoos/ArcGIS/all original layers WGS84 export")
corinewd <- c("C:/Users/samf/Documents/R/projects/cuckoos/data/corine raster 100x100")
cuckoowd <- c("C:/Users/samf/Documents/R/projects/cuckoos/data")
cuckoowd2 <- c("C:/Users/samf/Documents/R/projects/cuckoos/data/resampled data with uncertainty 100 bootstraps + distance, movement groups, etc")


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

### --- PLOT CUCKOO AUTUMN MIGRATION POINTS FROM EUROPE TO AFRICA --- ###

for (a in 1:31){
  setwd("C:/Users/samf/Documents/R/projects/cuckoos/data/resampled data with uncertainty 100 bootstraps + distance, movement groups, etc")
  cuckoofiles <- list.files()
  
  cuckoo <- read.csv(cuckoofiles[a], header=T)
  autumnonly <- Subset(cuckoo, julian >= 131)
  
  # long/lat are in CRS("+proj=longlat +datum=WGS84")
  # newlongs/newlats are in CRS("+init=epsg:3395")
  # centroidlong/lat are in CRS("+init=epsg:3395")
  
  coordinates(autumnonly) <- c("long","lat")
  proj4string(autumnonly) <- CRS("+proj=longlat +datum=WGS84") 
  
  # clipEurAfr <- drawExtent()
  #EurAfrmap <- crop(Eur.Afr,clipEurAfr)
  
  
  if (length(levels(as.factor(autumnonly@data$year)))==1) {
    colpoints <- factor(autumnonly@data$year, labels=c("black"))  
  } else if (length(levels(as.factor(autumnonly@data$year)))==2) {
    colpoints <- factor(autumnonly@data$year, labels=c("darkmagenta","black"))  
  } else {
    colpoints <- factor(autumnonly@data$year, labels=c("royalblue","darkmagenta","black"))  
  }
  
  # view map to see successful crossings (by year, different colours) and whether route was SE or SW
  
  setwd("C:/Users/samf/Documents/R/projects/cuckoos/output/Europe-Africa autumn maps")

  tiff(paste(levels(autumnonly$name)," - Europe-Africa.tiff"),res=300,height=2000,width=2000,units="px")
  par(mar=c(0,0,0,0))
  
  plot(EurAfrmap)
  plot(autumnonly, pch=16, col=as.character(colpoints), cex=1.2, add=T)
  
  dev.off()
}


### -------- Protected area network -------- ###

###------------------------------------------------------------###
#         Clip cuckoo stopovers to Europe map
###------------------------------------------------------------###

# create list to hold centroids of stopover points for all cuckoos

centroid.stopovers <- list()

# source file creates an object called "original", which are the tidied original datapoints for all individuals
setwd("C:/Users/samf/Documents/R/projects/cuckoos/scripts/")
source("cuckoo original data tidying source.R")

# for loop for each cuckoo file starts here
for (a in 1:31) {
  
  ### --- RESAMPLED cuckoo data: Load and add geoinfo --- ###
  
  setwd("C:/Users/samf/Documents/R/projects/cuckoos/data/resampled data with uncertainty 100 bootstraps + distance, movement groups, etc")
  cuckoofiles <- list.files()
  
  cuckoo <- read.csv(cuckoofiles[a], header=T)
  
  # long/lat are in CRS("+proj=longlat +datum=WGS84")
  # newlongs/newlats are in CRS("+init=epsg:3395")
  # centroidlong/lat are in CRS("+init=epsg:3395")
  
  coordinates(cuckoo) <- c("newlongs","newlats")
  proj4string(cuckoo) <- CRS("+init=epsg:3395")
  cuckoo <- spTransform(cuckoo,CRS("+proj=longlat +datum=WGS84")) # transform newlongs/newlats from World Mercator/WGS84 to longlat WGS84
  
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
  coordinates(fullcuckoo) <- c("newlongs","newlats")
  proj4string(fullcuckoo) <- CRS("+proj=longlat +datum=WGS84")
  
  ### create Europe only (+ at sea points) and Africa only data sets
  Eurcuckoo <- fullcuckoo[which(fullcuckoo$continent=="Europe" | fullcuckoo$continent=="at sea"),]
  Afrcuckoo <- fullcuckoo[which(fullcuckoo$continent=="Africa"),]
  
  
  
  
  
  ### --- ORIGINAL cuckoo data: Load and add mgroup/mtype variables, geoinfo to  --- ###
  
  originalcuckoo <- original[[a]]
  
  # for each level of tcycle, look up the corresponding mgroup and mtype values
  firsts <- rep(NA,length(levels(as.factor(cuckoo$tcycle)))) # blank variable to hold row numbers of first row of data for each individual
  
  for (j in 1:length(levels(as.factor(cuckoo$tcycle)))){ # loop counts through tcycle levels
    firsts[j] <- min(which(cuckoo$tcycle==levels(as.factor(cuckoo$tcycle))[j])) # row number of first observation of cuckoo j
  }
  
  orig.mgroup <- cuckoo$mgroup[firsts]
  orig.mtype <- cuckoo$mtype[firsts]
  
  orig.mvmts <- data.frame(mgroup=orig.mgroup,mtype=orig.mtype)
  
  originalcuckoo$mgroup <- factor(originalcuckoo$tcycle, labels=orig.mgroup)
  originalcuckoo$mtype <- factor(originalcuckoo$tcycle, labels=orig.mtype)
  
  # to get right number of levels for original mgroup/mtype, convert from factor, to character, to integer, back to factor!
  originalcuckoo$mgroup <- as.factor(as.integer(as.character(originalcuckoo$mgroup)))
  originalcuckoo$mtype <- as.factor(as.character(originalcuckoo$mtype))
  
  # long/lat are in CRS("+proj=longlat +datum=WGS84")
  coordinates(originalcuckoo) <- c("long","lat")
  proj4string(originalcuckoo) <- CRS("+proj=longlat +datum=WGS84")
  
  ### over() retrieves values of Eur.Afr at points described by cuckoo
  # extracts continent and administrative country names from Eur.Afr
  geoinfo.orig <- over(originalcuckoo, Eur.Afr)[,c("REGION","ADMIN")]
  colnames(geoinfo.orig) <- c("continent", "country")
  geoinfo.orig <- droplevels(geoinfo.orig)
  
  ### add a new factor for continent/country NAs, which should be points at sea
  atsea.pts.orig <- which(is.na(geoinfo.orig$continent))
  
  continent.factors.orig <- c(levels(geoinfo.orig$continent), "at sea")
  levels(geoinfo.orig$continent) <- continent.factors.orig
  geoinfo.orig$continent[atsea.pts.orig] <- "at sea"
  
  country.factors.orig <- c(levels(geoinfo.orig$country),"at sea")
  levels(geoinfo.orig$country) <- country.factors.orig
  geoinfo.orig$country[atsea.pts.orig] <- "at sea"
  
  ### add geographic info to cuckoo dataset and redefine projection of newlongs/newlats
  fulloriginal <- data.frame(originalcuckoo, geoinfo.orig)
  coordinates(fulloriginal) <- c("long","lat")
  proj4string(fulloriginal) <- CRS("+proj=longlat +datum=WGS84")
  
  ### create Europe only (+ at sea points) and Africa only data sets
  Euroriginal <- fulloriginal[which(fulloriginal$continent=="Europe" | fulloriginal$continent=="at sea"),]
  Afroriginal <- fulloriginal[which(fulloriginal$continent=="Africa"),]
  
  ###------------------------------------------------------------###
  #         Trim European data - autumn stopovers only
  ###------------------------------------------------------------###
  
  # subset out only stopovers in Europe
  # do not include the first "stopover", which we assume is the capture and therefore breeding location
  # ???????????? other early stopovers (e.g. #2, or #3) may also need to be removed if they constitute movements around the breeding area (unlikely though, since distance threshold for a "stopover" is < 20km, and anything > 20km would likely consitute a migration movement away from the breeding area)
  # include only southward migration stopovers (month >= June 1)
  
  ### RESAMPLED CUCKOO DATA
  
  Eurcuckoo$mgroup <- as.factor(Eurcuckoo$mgroup)
  stopovers <- Eurcuckoo[which(Eurcuckoo$mtype=="C" | Eurcuckoo$mtype=="S" | is.na(Eurcuckoo$mtype)),] # includes capture (breeding) location, stopovers, as well as last known location (if in Europe)
  #stopovers <- stopovers[which(stopovers$mgroup!=1),] # this line removes the capture occasion
  autumn <- stopovers[which(stopovers$julian >= 131),] # all julian dates after May 10 (removes spring stopovers)
  autumn@data <- droplevels(autumn@data)
  
  autumn$mgroup <- as.factor(autumn$mgroup)
  
  # can add other subsetting variables here, and change whatever "dat" is called
  
  dat <- autumn
  
  # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
  
  newdat <- spTransform(dat, CRS=corine.crs)
  
  ### ORIGINAL CUCKOO DATA
  
  Euroriginal$mgroup <- as.factor(Euroriginal$mgroup)
  stopovers.orig <- Euroriginal[which(Euroriginal$mtype=="C" | Euroriginal$mtype=="S" | is.na(Euroriginal$mtype)),] # includes capture (breeding) location and stopovers
  # stopovers.orig <- stopovers.orig[which(stopovers.orig$mgroup!=1),] # this line removes the capture occasion
  autumn.orig <- stopovers.orig[which(stopovers.orig$julian >= 131),]  # all julian dates after May 10 (removes spring stopovers)
  autumn.orig@data <- droplevels(autumn.orig@data)
  
  autumn.orig$mgroup <- as.factor(autumn.orig$mgroup)
  
  # can add other subsetting variables here, and change whatever "dat" is called
  
  dat.orig <- autumn.orig
  
  # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
  
  newdat.orig <- spTransform(dat.orig, CRS=corine.crs)
  
  
  ###------------------------------------------------------------###
  #         Plot individual cuckoo stopovers on Europe raster map 
  ###------------------------------------------------------------###
  
  
  # create new directory for specific bird to put output maps in, and setwd() as this new directory
  
  dir.create(paste("C:/Users/samf/Documents/R/projects/cuckoos/output/Europe maps/", levels(newdat$name), sep=""))
  outputdir <- paste("C:/Users/samf/Documents/R/projects/cuckoos/output/Europe maps/", levels(newdat$name), sep="")
  setwd(outputdir)
  
  # export Europe autumn stopovers to a tiff file, with different colours for different years
  
  if (length(levels(as.factor(newdat@data$year)))==1) {
    colpoints <- factor(newdat@data$year, labels=c("black"))  
  } else if (length(levels(as.factor(newdat@data$year)))==2) {
    colpoints <- factor(newdat@data$year, labels=c("darkmagenta","black"))  
  } else {
    colpoints <- factor(newdat@data$year, labels=c("royalblue","darkmagenta","black"))  
  }
  
  tiff(paste(levels(newdat$name)," - Europe.tiff"),res=300,height=2000,width=2400,units="px")
  
  plot(Europeraster)
  plot(newdat, pch=16, cex=0.8, col=as.character(colpoints), add=T)
  coords <- cbind(2700000,2700000)
  textpts <- SpatialPoints(coords)
  proj4string(textpts) <- corine.crs
  text(textpts, levels(newdat$name))
  
  
  dev.off()
  
  ###-----------------------------------------------###
  #         Draw MCP around stopover points
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
  
  # ??? can convert below loop to a function
  
  for (i in 1:length(MCP)){
    polycoords <- MCP@polygons[[i]]@Polygons[[1]]@coords
    MCPbox[[i]] <- c(min(polycoords[,1])-10000, max(polycoords[,1])+10000, min(polycoords[,2])-10000, max(polycoords[,2])+10000)
    MCPbox
  }
  
  clipMCP <- function(MCPminmax){
    MCPbounds <- extent(MCPminmax)
    return(MCPbounds)
  }
  
  mcp.bbox <- lapply(MCPbox,clipMCP)
  
  ### --- CENTROIDS of stopover points & identifying LAST stopover --- ###
  
  extractcentroid <- function(spdf){ #function to extract centroid values from the polygons (list components) of a SpatialPolygonsDataFrame
    centroid <- spdf@labpt
    return(centroid)
  }
  
  extractid <- function(spdf){ #function to extract id name (mgroup)
    id <- spdf@ID
    return(id)
  }
  
  centroids <- t(sapply(MCP@polygons,extractcentroid)) # represents centroids of MCPs describing the stopovers
  stopovers <- t(sapply(MCP@polygons,extractid)) # represents centroids of MCPs describing the stopovers
  
  centroids <- as.data.frame(centroids)
  stopovers <- as.factor(stopovers)
  
  cent.stop <- data.frame(stopovers,centroids)
  colnames(cent.stop) <- c("mgroup","centroidlong","centroidlat")
  
  bird.name <- rep(levels(newdat$name), nrow(cent.stop))
  
  bird.year <- rep(NA,length(mcp.bbox))
  for (i in 1:length(mcp.bbox)){
    bird.year[i] <- newdat@data[which(newdat@data$mgroup==MCP$id[i])[1],"year"]
  }
  
  centroid.stopovers[[a]] <- data.frame(name=bird.name, year=bird.year, cent.stop)
  
  # plot the new cropped rasters for each stopover movement, with the relevant stopover points and MCP on it
  # export plots to a TIFF
  
  for (i in 1:length(mcp.bbox)) {
    croppedraster <- crop(Europeraster,mcp.bbox[[i]])
    
    tiff(paste(levels(newdat$name), " - movement ", MCP$id[i], ".tiff", sep=""),res=150,height=800,width=800,units="px")
    
    plot(croppedraster)
    plot(newdat[which(newdat$mgroup==MCP$id[i]),], pch=16, cex=0.8, col="black", add=T)
    
    # polygontoplot <- SpatialPolygons(list(MCP@polygons[[i]])) # MCP around stopover points
    # plot(polygontoplot, lwd = 2, add=T)
    coords <- cbind(MCPbox[[i]][1]+500,MCPbox[[i]][4]-3000)
    textpts <- SpatialPoints(coords)
    proj4string(textpts) <- corine.crs
    
    plottext <- paste(levels(newdat$name), " - #", MCP$id[i], " - ", newdat@data[which(newdat@data$mgroup==MCP$id[i])[1], "country"], " ", newdat@data[which(newdat@data$mgroup==MCP$id[i])[1], "year"], sep="")
    text(textpts, plottext, cex=0.8, pos=4)
    
    # plot original data points (darkmagenta) from transmission cycle over resampled points (black) from transmission cycle for each stopover
    plot(newdat.orig[which(newdat.orig$mgroup==MCP$id[i]),], pch=16, cex=0.8, col="darkmagenta", add=T)
    
    dev.off()
  }
}


##################################################
##################################################

###-------------------------------------------------------------###
#    Plot ALL cuckoo centroid stopover points on Europe map
###-------------------------------------------------------------###

setwd("C:/Users/samf/Documents/R/projects/cuckoos/output/Europe maps/")

tiff("all stopover centroids - Europe.tiff",res=300,height=2000,width=2400,units="px")

plot(Europeraster)

# for loop through list(centroid.stopovers), add individual's stopovers to the plot(Europeraster)
for (a in 1:31){
  
  usedat <- centroid.stopovers[[a]]
  
  # identify last stopovers before Sahara crossing
  # for each level of year, identify max(mgroup)
  
  yearlevels <- levels(as.factor(usedat$year))
  maxmgroup.byyear <- rep(NA, length(yearlevels))
  
  for (i in 1:length(yearlevels)){
    maxmgroup.byyear[i] <- max(as.integer(as.character(usedat[usedat$year==yearlevels[i], "mgroup"]))) # extract max mgroup number for each level of year
  }
  

  laststop <- factor(levels=c("N","Y"))
  laststop[which(as.integer(as.character(usedat$mgroup)) %in% maxmgroup.byyear)] <- "Y"
  laststop[which(!as.integer(as.character(usedat$mgroup)) %in% maxmgroup.byyear)] <- "N"
  usedat <- data.frame(usedat,laststop)
  
  if (length(usedat$laststop) == 1) {
    colpoints <- factor(usedat$laststop, labels=c("darkmagenta")) 
  } else {
    colpoints <- factor(usedat$laststop, labels=c("black","darkmagenta")) 
  }
  
  
  coordinates(usedat) <- c("centroidlong","centroidlat")
  proj4string(usedat) <- corine.crs
  
  plot(usedat, pch=16, cex=0.8, col=as.character(colpoints), add=T)
  # coords <- cbind(2700000,2700000)
  # textpts <- SpatialPoints(coords)
  # proj4string(textpts) <- corine.crs
  # text(textpts, levels(newdat$name))
  
}

dev.off()

##################################################
##################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################



###-----------------------------------------------###
#         Draw MCP around stopover points
###-----------------------------------------------###

# for each movement group of each cuckoo, draw a minimum convex polygon around the mgroup's points
# clip the raster to a bounding box with an extent that buffers the MCP by ~ 5km

# use mcp() to estimate home range
newdat@data <- droplevels(newdat@data)
stoppoints <- newdat[,c("mgroup")]
MCP.95 <- mcp(stoppoints, percent=95)

plot(MCP.95, colborder = "blue", lwd = 2, add=T)


plot(r)
boxScotland <- drawExtent()
boxEngland <- drawExtent()

system.time({newraster <- crop(r,box)})

plot(newraster)
boxScotland <- drawExtent()
boxEngland <- drawExtent()

rasterEng <- crop(r, boxEngland)
rasterScot <- crop(r, boxScotland)

plot(rasterEng)
plot(newdat, pch=16, cex=1, col = "blue" , add=T)
plot(MCP.95, colborder = "blue", lwd = 2, add=T)

plot(rasterScot)
plot(newdat, pch=16, cex=0.5, col = "black" , add=T)
plot(MCP.95, colborder = "blue", lwd = 2, add=T)


par(mfrow=c(1,1))
layout.show(layout(rbind(c(1,1))))
plot(rasterEng)
plot(newdat, pch=16, cex=0.3, col = "black" , add=T)
plot(BBMCP.90, colborder = "black", lwd = 2, add=T)

