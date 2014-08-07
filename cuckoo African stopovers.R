##########################################################
#
#	Cuckoo African stopovers
#
# Code extracts country & continent geographic info for cuckoo resampled location points, clips dataset to only points occurring in Africa, creates a MCP around each stopover
#
# Combines all stopover polygons for cuckoos with Africa locations and exports to a shapefile (also can export separate shapefiles for each individual cuckoo)
#
#	Samantha Franks
# 10 Dec 2013
# 23 Dec 2013 - rerun with re-edited resampled dat
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

###-----------------------------------------------------------###
#         Load base layers
###-----------------------------------------------------------###

### ------- Europe/Africa national boundaries layer ------- ###

dat.world <- getMap(resolution="high")

Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]

Eur.Afr <- spTransform(Eur.Afr,CRS("+proj=longlat +datum=WGS84")) # transform map to exactly the same projection as cuckoo newlongs/newlats data

Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]

# clipAfr <- drawExtent()
# Afrmap <- crop(Afr.only,clipAfr)

#writeOGR(Afrmap, dsn = ".", layer = "Africa map", driver = "ESRI Shapefile")

###------------------------------------------------------------###
#         Clip cuckoo stopovers to Europe map
###------------------------------------------------------------###

# create list to hold centroids of stopover points for all cuckoos

AfrMCPstopovers <- list()
Afrdata <- list()

# for loop for each cuckoo file starts here
for (a in 1:31) {
  
  ### --- RESAMPLED cuckoo data: Load and add geoinfo --- ###
  
  setwd("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/data/resampled data with uncertainty 100 bootstraps")
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
  
  if (nrow(Afrcuckoo@data) != 0) {
    
  ###------------------------------------------------------------###
  #         Trim Africa data - stopovers only
  ###------------------------------------------------------------###
  
  # subset out only stopovers in Africa (includes winter "stopovers")
   
  ### RESAMPLED CUCKOO DATA
  
  Afrcuckoo$mgroup <- as.factor(Afrcuckoo$mgroup)
  stopovers <- Afrcuckoo[which(Afrcuckoo$mtype=="S" | is.na(Afrcuckoo$mtype)),]
  
  stopovers@data <- droplevels(stopovers@data)
  
  # can add other subsetting variables here, and change whatever "dat" is called
  
  # remove stopovers (mgroups) with < 2 transmission cycles
  
  # for each level of mgroup, check the number of levels of tcycle - if < 2, then Subset data to remove that mgroup
  
  newstopovers <- stopovers
  
  for (i in 1:length(levels(stopovers$mgroup))) {
    temp <- newstopovers[newstopovers$mgroup==levels(stopovers$mgroup)[i],]
    temp@data <- droplevels(temp@data)
    if (length(levels(as.factor(temp$tcycle))) < 2) {
      newstopovers <- newstopovers[which(newstopovers$mgroup!=levels(stopovers$mgroup)[i]),]
      newstopovers@data <- droplevels(newstopovers@data)
    } else {
      newstopovers <- newstopovers
    }
  }
  
  newdat <- spTransform(newstopovers, CRS("+init=epsg:3395"))

  # write African data to list, to write to csv later
  
  Afrdata[[a]] <- newdat
  
  ###------------------------------------------------------------###
  #         Plot individual cuckoo stopovers on Africa map 
  ###------------------------------------------------------------###
  
  

  
  ###### set directory for export of all to one folder

   outputdir <- c("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/output/Africa stopovers/")
   setwd(outputdir)
  
  # export Africa stopovers to a tiff file, with different colours for different years
  
  if (length(levels(as.factor(newdat@data$year)))==1) {
    colpoints <- factor(newdat@data$year, labels=c("black"))  
  } else if (length(levels(as.factor(newdat@data$year)))==2) {
    colpoints <- factor(newdat@data$year, labels=c("darkmagenta","black"))  
  } else {
    colpoints <- factor(newdat@data$year, labels=c("royalblue","darkmagenta","black"))  
  }
  
#   tiff(paste(levels(newdat$name)," - Africa.tiff"),res=300,height=1000,width=1200,units="px")
#   
#   par(mar=c(1,1,1,1))
#   plot(Afrmap)
#   plot(newdat, pch=16, cex=0.8, col=as.character(colpoints), add=T)
#   coords <- cbind(-10,0)
#   textpts <- SpatialPoints(coords)
#   proj4string(textpts) <- CRS("+proj=longlat +datum=WGS84")
#   text(textpts, levels(newdat$name))
#   
#   
#   dev.off()
  
  ###-----------------------------------------------###
  #         Draw MCP around stopover points
  ###-----------------------------------------------###
  
  # for each movement group of each cuckoo, draw a MCP around the mgroup's points using mcp() to estimate home range
  # choose points for MCP such that all points are included (percent = 100)
  newdat@data <- droplevels(newdat@data)
  stoppoints <- newdat[,c("mgroup")]
  MCP <- mcp(stoppoints, percent=100)
  
  MCP@data <- data.frame(name=rep(levels(newdat$name),nrow(MCP@data)),MCP@data)
  colnames(MCP@data) <- c("name","mgroupid","MCParea")
  
  AfrMCPstopovers[[a]] <- MCP
  
#   writeOGR(MCP, dsn = ".", layer = paste(levels(newdat$name)," - Africa stopover MCPs", sep=""), driver = "ESRI Shapefile")
#   list.files(pattern = "stopover MCPs")
  } else {}
  
}
  
AfrMCPstopovers2 <- AfrMCPstopovers[!unlist(lapply(AfrMCPstopovers, is.null))]

Afrdata.complete <- Afrdata[!unlist(lapply(Afrdata, is.null))]

Afrdata.full <- do.call(rbind, Afrdata.complete)

write.csv(Afrdata.full, "Africa resampled location points.csv", row.names=F)

for (i in 1:length(AfrMCPstopovers2)){
AfrMCPstopovers2[[i]] <- spChFIDs(AfrMCPstopovers2[[i]], paste(AfrMCPstopovers2[[i]]$name, AfrMCPstopovers2[[i]]$mgroupid))
}

allAfrstopovers <- do.call(rbind, AfrMCPstopovers2)

allAfrstopovers2 <- spTransform(allAfrstopovers, CRS("+init=epsg:4326"))

writeOGR(allAfrstopovers2, dsn = ".", layer = "all birds - Africa stopover MCPs - 20131223 - m2", driver = "ESRI Shapefile")

#############################
#############################
#############################


