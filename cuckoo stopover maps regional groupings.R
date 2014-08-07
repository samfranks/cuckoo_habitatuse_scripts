##############################################################################
#
#   PLOT CUCKOO STOPOVER LOCATIONS
#
#   Samantha Franks
#   23 Jun 2014
#
##############################################################################


library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
library(rgeos)
library(raster)
library(RgoogleMaps)
library(plotKML)

parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

resampwd <- ("/resampled data + extra grouping variables")
origwd <- ("/original data + extra grouping variables")

corinewd <- paste(parentwd, "/data/corine raster 100x100", sep="")
elevationwd <- paste(parentwd, "/data/elevation data", sep="")

corine.crs <- CRS("+init=epsg:3035")
epsg4326 <- CRS("+init=epsg:4326")

  setwd(parentwd)
  
  #### for loop for each cuckoo file starts here ####
  
  stopSPDFs4326 <- list()
  
  for (a in 1:31) { # START LOOP 
    
    ####--- LOAD DATA, ADD PROJECTION INFO, TRANSFORM ---####
    
    setwd(paste(datawd,resampwd,sep=""))
    
    dataset <- read.csv(list.files()[a], header=T)
    
    ### check that dataset is not an individual to exclude (Idemili & Karma); break out of current loop run and continue to next loop level
    if (dataset$name[1] == "Karma" | dataset$name[1] =="Idemili") {
      next
    }
    
    # convert mgroup to factor
    dataset$mgroup <- as.factor(dataset$mgroup)
    
    # create new concatenated variable, name.mgroup
    name.mgroup <- as.factor(paste(dataset$name,dataset$mgroup,sep=""))
    dataset2 <- data.frame(name.mgroup,dataset)
    
    # subset dataset by stopover sites (no breeding sites)
    dataset3 <- subset(dataset2, stopoversite=="Y")
    dataset3 <- droplevels(dataset3)
    
    # check if subsetted dataset has no stopovers; if no stopovers, break out of current count of "a" for loop and skip to next "a"
    if (nrow(dataset3) == 0) {
      next
    }
      
    ### convert to spatial points dataframe
    # long/lat (original) are in CRS("+init=epsg:4326")
    # long/lat (resampled) are in CRS("+init=epsg:3395")  
    # newlongs/newlats (resampled) are in CRS("+init=epsg:4326")
    # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
    
    coordinates(dataset3) <- c("newlongs","newlats")
    proj4string(dataset3) <- CRS("+init=epsg:4326")
    newdat <- spTransform(dataset3, CRS=corine.crs)
    
    ####--- CREATE STOPOVER POLYGONS ---####
      
    mgroupSPs <- list()
    
    for (i in 1:length(levels(newdat$mgroup))) {
      datsub <- subset(newdat, mgroup==levels(newdat@data$mgroup)[i])
      datsub@data <- droplevels(datsub@data)
      #mgroupSPs[[i]] <- gBuffer(datsub, width=500, id=levels(newdat@data$mgroup)[i])
      mgroupSPs[[i]] <- gBuffer(datsub, width=500, id=levels(datsub@data$name.mgroup))
    }
    
    allstops <- do.call(rbind, mgroupSPs)
    mgroupid <- sapply(slot(allstops, 'polygons'), function(i) slot(i, 'ID'))
    allstopSPDFs <- SpatialPolygonsDataFrame(allstops, data.frame(mgroup=mgroupid, row.names=mgroupid))
    
    ####--- CREATE .KML FILES FOR INDIVIDUAL ---####
    
    colours <- rainbow(28)
    
    ### KML POINTS
    
    # subset to 1 point per stopover
    x <- which(duplicated(newdat@data[,c("name.mgroup")])==FALSE)
    datreduced <- newdat[x,]
    
    # convert projection to epsg4326
    dat.kml <- spTransform(datreduced, CRS=epsg4326)
    
    setwd(paste(outputwd,"/stopover points and polygons/points",sep=""))  
    plotKML(dat.kml, paste(levels(dat.kml$name), "points.kml",sep=""), colour_scale=colours[a])
    #writeOGR(dat, paste(levels(dat$name), "points.kml",sep=""), paste(levels(dat$name), "points", sep=""), driver="KML")
    
    ### KML POLYGONS
    setwd(paste(outputwd,"/stopover points and polygons/polygons",sep=""))
    stopSPDFs4326[[a]] <- spTransform(allstopSPDFs, CRS=epsg4326)
    plotKML(stopSPDFs4326[[a]], paste(levels(dat.kml$name), "polygons.kml",sep=""), colour_scale=colours[a], alpha=0.75)
    #writeOGR(stopSPDFs4326[[a]], dsn=".", layer=paste(levels(dat$name), "polygons", sep=""), driver="ESRI Shapefile")
               
  }

####--- OUTPUT ALL STOPOVER POLYGONS ALL BIRDS ---####

allpolygons <- do.call(rbind, stopSPDFs4326)
#allpolygons <- rbind(stopSPDFs4326[[1]], stopSPDFs4326[[2]])

setwd(paste(outputwd,"/stopover points and polygons/polygons",sep=""))
writeOGR(allpolygons, dsn=".", layer="allpolygons", sep=""), driver="ESRI Shapefile")
kml(stopSPDFs4326, "allpolygons.kml", colour_scale=colours[a], fill=colours[a])



