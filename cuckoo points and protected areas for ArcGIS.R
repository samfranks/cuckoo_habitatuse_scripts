##############################################################################
#
#   PLOT CUCKOO STOPOVER LOCATIONS
#
#   Samantha Franks
#   23 Aug 2014
#
##############################################################################


library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
library(rgeos)
library(raster)
library(RgoogleMaps)
library(plotKML)
library(reshape)


cluster <- FALSE

Mac <- FALSE

randomradius <- 500

if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  GISwd <- c("C:/Users/samf/Documents/GIS/cuckoos/protected areas")
  
}


if (cluster) parentwd <- c("/users1/samf/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")


####==== IMPORT POINT DATA ====####


setwd(paste(datawd, "/data for analysis/", sep=""))
present <- read.csv("presence data all variables points.csv", header=T)
absent <- read.csv(paste("absence data all variables points ", randomradius, " km.csv", sep=""), header=T)

present <- rename(present, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))

absent <- rename(absent, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))

newpresent <- rename(present, c("newlongs.epsg3035"="long.epsg3035", "newlats.epsg3035"="lat.epsg3035", "newlongs.epsg4326"="long.epsg4326", "newlats.epsg4326"="lat.epsg4326"))
newabsent <- rename(absent, c("randomlong.epsg3035"="long.epsg3035", "randomlat.epsg3035"="lat.epsg3035", "randomlong.epsg4326"="long.epsg4326", "randomlat.epsg4326"="lat.epsg4326"))

############# subset out only one round of absences
#newabsent <- subset(newabsent, nabsence==1)
newabsent <- newabsent[,-which(names(newabsent) %in% c("nabsence","random.scale"))]

alldata <- rbind(newpresent, newabsent)


####==== IMPORT PROTECTED AREA DATA ====####

PA.epsg4326 <- readOGR(GISwd, "terrestrial_PAs_Europe_GCS_WGS_1984_final")
natPA <- subset(PA.epsg4326, !grepl("Ramsar | Directive", desig_eng))
internatPA <- subset(PA.epsg4326, grepl("Ramsar | Directive", desig_eng))

####==== WRITE DATA FOR USE IN ARCGIS ====####

setwd(paste(outputwd,"/stopover points and polygons/", sep="")) 
writeOGR(internatPA, dsn=".", layer="international_PAs_epsg4326", driver="ESRI Shapefile")
writeOGR(natPA, dsn=".", layer="national_PAs_epsg4326", driver="ESRI Shapefile")

coordinates(newpresent) <- c("long.epsg4326", "lat.epsg4326")
proj4string(newpresent) <- CRS("+init=epsg:4326")

coordinates(newabsent) <- c("long.epsg4326", "lat.epsg4326")
proj4string(newabsent) <- CRS("+init=epsg:4326")

writeOGR(newpresent, dsn=".", layer="data_for_analysis_presences_epsg4326", driver="ESRI Shapefile")
writeOGR(newabsent, dsn=".", layer="data_for_analysis_absences_epsg4326", driver="ESRI Shapefile")