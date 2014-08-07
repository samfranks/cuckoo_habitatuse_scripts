##########################################################
#
#  Cuckoo land cover and protected area extractions
#
#	Samantha Franks
#	28 Nov 2013
# 17 Dec 2013 - a) run with corrected mgroups + mtypes from rerun resampling data; 2) now also includes code for to run originl data
# 24 Dec 2013 - now removes mgroups with points from only a single tcycle
# 31 Mar 2014
#
##########################################################

###--- TWO APPROACHES FOR PROTECTED AREA ANALYSIS

### APPROACH 1: do cuckoo locations fall within protected areas?
# Extract protected area information from actual cuckoo location points (either original or resampled)

### APPROACH 2: do cuckoo stopover polygons overlap with protected areas, and if so, what proportion of their area overlaps?
# Extract protected area overlap information with intersection of stopover polygons and protected area polygons

# Approach 2 address the question of "what type of habitat is generally comprising the landscape in which cuckoos are stopping?" rather than specifically what habitat THEY ARE USING.  This approach captures the amount of protected area in the landscape of the stopover polygon, but not necessarily its actual use

# Approach 1 will better capture the ACTUAL USE of protected areas

# BUT, Approach 2 will better capture use of edge habitat around protected areas, which will be missed by Approach 1

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
library(PBSmapping)

# for (a in 1:31) {
# source("C:/Users/samf/Documents/Git/cuckoos/scripts/source code to extract time and location-specific cuckoo points.R")
# }

if (!("Europeraster" %in% ls())) { ### run prep code to load base layers if not in workspace already
  
  ###--------------------------------------------------------###
  #        LOAD BASE LAYERS
  ###--------------------------------------------------------###
  
  ### CORINE layer ###
  
  # projection and datum information for all maps
  # epsg: 3395 - +proj=merc +ellps=wGS84 +datum=WGS84
  # epsg: 4326 - +proj=latlong +ellps=wGS84 +datum=WGS84
  # epsg: 3035 - LAEA CRS, used for the Corine Land Cover 2006 data
  
  corine.crs <- CRS("+init=epsg:3035")
  
  ### ------- CORINE land cover raster file ------- ###
  corinewd <- c("C:/Users/samf/Documents/Git/cuckoos/data/corine raster 100x100")
  
  setwd(corinewd)
  r <- raster("g100_06.tif")
  
  Europebox <- extent(2483500,5890375,1276625,4286750) # bounding box extent for relevant parts of Europe (no Scandanavia or Turkey)
  Europeraster <- crop(r,Europebox)
  
}

### ------- Europe/Africa national boundaries layer ------- ###

if (!("Eur.Afr" %in% ls())) { ### run prep code to load base layers if not in workspace already
  
  dat.world <- getMap(resolution="high")
  
  Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
  Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
  Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]
  
  Eur.Afr <- spTransform(Eur.Afr,CRS("+proj=longlat +datum=WGS84")) # transform map to exactly the same projection as cuckoo newlongs/newlats data
  
  Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
  Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]
  
}
  
### ------- Protected area GIS layer (raster or shapefile) ------- ###

if (!("PAs" %in% ls())) { ### run prep code to load base layers if not in workspace already
  
  #   GISwd <- c("D:/Sam Franks/GIS/cuckoos")
  #   setwd(GISwd)
  #   
  #   PA <- raster("Europe PA raster 50m x 50m.tif")
  #   
  #   colors <- c("white",rep("blue",89662))
  #   PA@legend@colortable <- colors
  #   # default of rpoly colortable is logical(0)
  
  #if (!cluster)
  GISwd <- c("D:/Sam Franks/GIS/cuckoos")
  
  #if (cluster)
  #  GISwd <- c("/users1/samf/cuckoos")
  
  PAs <- readOGR(GISwd, "terrestrial PAs mainland W Europe corine countries only EPSG 3035")
  PAs <- spTransform(PAs,CRS("+init=epsg:3035"))
  
  PA.subset <- subset(PAs, desig_eng!="UNESCO-MAB Biosphere Reserve")
  PA.subset@data <- droplevels(PA.subset@data)
  
}

###-----------------------------------------------------------###
#         SET WORKING DIRECTORIES
###-----------------------------------------------------------###

setwd("C:/Users/samf/Documents/Git/cuckoos")
parentwd <- getwd()

datawd <- ("C:/Users/samf/Documents/Git/cuckoos/data")
outputwd <- ("C:/Users/samf/Documents/Git/cuckoos/output")

resampwd <- ("/resampled data + extra grouping variables")
origwd <- ("/original data + extra grouping variables")


##############################################################
##############################################################
#
# EXTRACT CORINE & PA VALUES FOR LOCATION POINTS
#
##############################################################
##############################################################

### TO RUN ORIGINAL or RESAMPLED data, change below line

original <- TRUE # or FALSE if resampled

# create list to hold datasets with extracted corine and protected area values for each cuckoo

corine.values <- list()
PA.extractedinfo <- list()
corine.pa.values <- list()
PAextractinfotime <- list()


# for loop for each cuckoo file starts here
for (a in 1:31) {
  
  ###--- LOAD DATA, ADD PROJECTION INFO, TRANSFORM ---###
  
  if (original) setwd(paste(datawd,origwd,sep=""))
  if (!original) setwd(paste(datawd,resampwd,sep=""))
    
  dataset <- read.csv(list.files()[a], header=T)
  
  # long/lat are in CRS("+proj=longlat +datum=WGS84")
  # newlongs/newlats are in CRS("+init=epsg:3395")
  # centroidlong/lat are in CRS("+init=epsg:3395")
  # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
  
  # ---------------------- RESAMPLED DATA -----------------------#
  if (!original) {
    coordinates(dataset) <- c("newlongs","newlats")
    proj4string(dataset) <- CRS("+init=epsg:3395") 
    newdataset <- spTransform(dataset, CRS=corine.crs) 
  }
  
  # ---------------------- ORIGINAL DATA -----------------------#
  if (original) {
    coordinates(dataset) <- c("long","lat")
    proj4string(dataset) <- CRS("+init=epsg:4326")
    newdataset <- spTransform(dataset, CRS=corine.crs)
  }
    
  ###--- LAND-COVER POINT VALUE EXTRACTIONS ---###
  
  corine.values[[a]] <- data.frame(newdataset, corine.values=extract(Europeraster, newdataset@coords))
  
 ###--- PROTECTED AREA POINT VALUE EXTRACTIONS ---###
  
  PAextractinfotime[[a]] <- system.time({
    
    PA.extractedinfo[[a]] <- over(newdataset, PA.subset)[,c("wdpaid","country","name","desig_eng","desig_type","iucn_cat","rep_area")]
    
    colnames(PA.extractedinfo[[a]]) <- c("wdpaid", "sitename","desig_eng","sitecountry", "sub_loc", "desig_type", "iucn_cat", "rep_area")
    
    corine.pa.values[[a]] <- data.frame(corine.values[[a]], PA.extractedinfo[[a]])
    
  })

  
} # END LOOP through each cuckoo file

# set directory of csv output dataset with corine values
if (!original) setwd(paste(datawd,"/resampled data corine and PA extracted values",sep=""))
if (original) setwd(paste(datawd,"/original data corine and PA extracted values",sep=""))


for (i in 1:31){
   write.csv(corine.pa.values[[i]], paste(corine.pa.values[[i]]$name[1],".csv", sep=""), row.names=FALSE)
}


