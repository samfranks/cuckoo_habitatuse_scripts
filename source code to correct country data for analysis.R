##############################################################################
#
#   ADD CORRECT COUNTRY DATA FOR ABSENCES AND CORRECT AT SEA POINTS COUNTRY DATA
#   Samantha Franks
#   9 Sep 2014
#
##############################################################################

####==== NOTES ====####

# correction requires a combination of R and ArcGIS
#--- in R ---#
# 1) load ESRI shapefile layer of countries that is used in ArcGIS project (from G:/Common Themes/ESRIData)
# 2) extract correct country names from this shapefile for ABSENT DATA only (absent data country variable is currently a result of the presence point that an absence was originally generated from): use over() in rgdal library
# 3) add new country data to absent data, and export as shapefile into ArcGIS

#--- in ArcGIS ---#
# 4) add shapefile layer of absences with new country data into ArcGIS

#=== for both absences and presences ===###
# 5) select only those observations where the new country field is a blank (NA) (or for presences, country=at sea), and export to new shapefile
# 6) with country blank-only data, use the Spatial Join tool in Overlay Toolbox and the Intersect match option, with a search radius of 10km
# 7) output the new shapefile as a "buffered" shapefile

library(sp)
library(rgeos)
library(rgdal)

cluster <- FALSE

Mac <- FALSE

present <- FALSE

randomradius <- 500

if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  
}


if (cluster) parentwd <- c("/users1/samf/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

####==== ADD CORRECT COUNTRY DATA ====####

### NOTES

### load ESRI shapefile layer of countries that is used for ArcGIS layer in map of presences/absences (from G:/Common Themes/ESRIData)
epsg4326 <- CRS("+init=epsg:4326")
country <- readOGR(dsn="C:/Users/samf/Documents/GIS/cuckoos", layer="ESRI_country_shapefile")
country <- spTransform(country, epsg4326)

### convert points to spatial points

### SPATIAL JOIN WITH BUFFER CAN BE IMPLEMENTED FOR PRESENCES IN ARCGIS
### COUNTRY EXTRACTION FOR ABSENCES IS FASTER IN R, SO EXTRACT NEW COUNTRIES USING ESRI COUNTRY SHAPEFILE IN R, THEN EXPORT TO ARCGIS

### Import present and absent datasets
setwd(paste(datawd, "/data for analysis/", sep=""))


if (present) {
  
  present <- read.csv("presence data all variables points.csv", header=T)
  
  present <- rename(present, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))
  
  present.sp <- present
  coordinates(present.sp) <- c("newlongs.epsg4326","newlats.epsg4326")
  proj4string(present.sp) <- epsg4326
  
  geoinfo.present <- over(present.sp, country)[,c("CNTRY_NAME", "SQKM_CNTRY")]
  geoinfo.present <- droplevels(geoinfo.present)
  geoinfo.present <- rename(geoinfo.present, c("CNTRY_NAME"="country2", "SQKM_CNTRY"="area.sqkm"))
  
  present.sp2 <- present.sp
  present.sp2@data <- data.frame(present.sp@data, geoinfo.present)
  presentDF <- as.data.frame(present.sp2)
}

if (!present) {
  absent <- read.csv(paste("absence data all variables points ", randomradius, " km.csv", sep=""), header=T)
  
  absent <- rename(absent, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))
  
  
  absent.sp <- absent
  coordinates(absent.sp) <- c("randomlong.epsg4326", "randomlat.epsg4326")
  proj4string(absent.sp) <- epsg4326
  
  geoinfo.absent <- over(absent.sp, country)[,c("CNTRY_NAME", "SQKM_CNTRY")]
  geoinfo.absent <- droplevels(geoinfo.absent)
  geoinfo.absent <- rename(geoinfo.absent, c("CNTRY_NAME"="country2", "SQKM_CNTRY"="area.sqkm"))
  
  absent.sp2 <- absent.sp
  absent.sp2@data <- data.frame(absent.sp@data, geoinfo.absent)
  absentDF <- as.data.frame(absent.sp2)
  
}

setwd(paste(outputwd, "/stopover points and polygons", sep=""))
writeOGR(absent.sp2, dsn=".", layer=paste("absences_", randomradius, "km_new_country_extraction", sep=""), driver="ESRI Shapefile")
# writeOGR(present.sp2, dsn=".", layer="presences_new_country_extraction", driver="ESRI Shapefile")

# ####==== ADD CORRECT COUNTRY DATA ====####
# 
# # absences country data is generated from presence country data, rather than actual country an absence falls in
# # code + work in ArcGIS also adds country info for at sea points by a buffered Spatial Join (not possible in R)
# 
# setwd(paste(parentwd, "/scripts", sep=""))
# source("source code to correct country data for analysis.R")


########################################################

############    do work in ArcGIS here   #################

########################################################

### re-import presence and absence data shapefiles from ArcGIS

if (present) {
  
  # 1) presence data, with at sea points only (now with country name)
  present.atsea <- readOGR(paste(outputwd, "/stopover points and polygons", sep=""), layer="data_for_analysis_presences_epsg4326_atsea_buffered")
  present.atsea2 <- as.data.frame(present.atsea[,3:length(present.atsea@data)])
  colnames(present.atsea2)[1:48] <- colnames(presentDF)[1:48]
  present.atsea2 <- rename(present.atsea2, c("presenc"="presence", "CNTRY_NAME"="country2", "SQKM_CNTRY"="area.sqkm", "coords.x1"="newlongs.epsg4326", "coords.x2"="newlats.epsg4326"))
  # remove at sea points from present.newcountry dataset and rbind present.atsea points
  present.atsea3 <- present.atsea2[,names(presentDF)]
  
  # from at sea points, want to use country2 data
  present.atsea4 <- data.frame(present.atsea3, usecountry=present.atsea3$country2)
  
  # from all other points, want to use country data
  presentDF2 <- data.frame(presentDF, usecountry=presentDF$country)
  
  present.all <- rbind(present.atsea4, subset(presentDF2, usecountry!="at sea"))
  present.all <- droplevels(present.all)
  
  setwd(paste(datawd, "/data for analysis/", sep=""))
  write.csv(present.all, "presence data all variables points all with country data.csv", row.names=FALSE)
  
}

# 2) absence data, at sea points only

if (!present) {
  
  absent.atsea <- readOGR(paste(outputwd, "/stopover points and polygons", sep=""), layer=paste("data_for_analysis_absences_", randomradius, "km_epsg4326_atsea_buffered", sep=""))
  absent.atsea2 <- as.data.frame(absent.atsea)
  absent.atsea2 <- absent.atsea2[,c(3:59)]
  absent.atsea3 <- absent.atsea2[,c(1:49, 56:57, 50:51, 54:55)]
  colnames(absent.atsea3) <- names(absentDF)
  
  # from at sea points, want to use country2 data
  absent.atsea4 <- data.frame(absent.atsea3, usecountry=absent.atsea3$country2)
  
  # from all other absence points, want to use country2 data
  absentDF2 <- data.frame(absentDF, usecountry=absentDF$country2)
  
  absent.all <- rbind(absent.atsea4, subset(absentDF2, usecountry!="NA"))
  absent.all <- droplevels(absent.all)
  
  setwd(paste(datawd, "/data for analysis/", sep=""))
  write.csv(absent.all, paste("absence data all variables points all with country data ", randomradius, " km.csv", sep=""), row.names=FALSE)
  
}