######################## HEADER ############################
#
#  Cuckoo land cover, protected area, elevation and climate extractions
#
#	Samantha Franks
#	28 Nov 2013
# 17 Dec 2013 - a) run with corrected mgroups + mtypes from rerun resampling data; 2) now also includes code for to run originl data
# 24 Dec 2013 - now removes mgroups with points from only a single tcycle
# 31 Mar 2014
# 2 May 2014: changed stopover polygons so that no longer MCPs, but are derived from 500m buffered points; added elevation raster data extraction
# 6-7 May 2014: added elevation data extraction code
# 13 May 2014: added climate data extraction code and is working to increase buffer at which climate data is extracted if NA cells are found
# 2 July 2014: changed climate data to extract SPEI data, bootstrap the number of absences, and output a ESRI shapefile of the bootstrapped absence polygons for each individual
#
###########################################################

################### TO DO #########################

# need to add data extractions to this script for:
#   - elevation data
#   - climate data

###--- TWO APPROACHES FOR PROTECTED AREA ANALYSIS

### APPROACH 1: do cuckoo locations fall within protected areas?
# Extract protected area information from actual cuckoo location points

### APPROACH 2: do cuckoo stopover polygons overlap with protected areas, and if so, what proportion of their area overlaps?
# Extract protected area overlap information with intersection of stopover polygons and protected area polygons

# Approach 2 address the question of "what type of habitat is generally comprising the landscape in which cuckoos are stopping?" rather than specifically what habitat THEY ARE USING.  This approach captures the amount of protected area in the landscape of the stopover polygon, but not necessarily its actual use

# Approach 1 will better capture the ACTUAL USE of protected areas

# BUT, Approach 2 will better capture use of edge habitat around protected areas, which will be missed by Approach 1

################### NOTES #########################

# CORINE land cover data 2006 (V.17 12/2013) is downloaded as a raster .tif from the European Environment Agency http://www.eea.europa.eu/data-and-maps/data/corine-land-cover-2006-raster-3

# Elevation raster data is downloaded from within R using getData("alt").  Data were aggregated from SRTM 90 m resolution data between -60 and 60 latitude, and is in EPSG 4326 (WGS84 lat/long) projection.  Relevant country rasters are downloaded in the Europe_elevation_raster.R file, merged, and then written to a .tif raster file (Europe elevation raster.tif)

# Climate data uses ENSEMBLES European climate data downloaded as NETCDF files from http://www.ecad.eu/download/ensembles/download.php
# uses the 15-year chunk datafiles on the 0.25 degree grid for daily mean temperature (TG) and daily precipitation sum (RR) for 1995-2013

# Protected area data uses vector shapefile downloaded for all of Europe from protectedplanet.net.  Shapefile was loaded into GIS first, filtered for relevant countries to make the file smaller (excluded Russia, Fennoscandia, and some of the Eastern European countries), converted to EPSG 3035 projection, then exported as a vector file and loaded into R, where it was pruned further to remove all UNESCO MAB biosphere reserves

#### LOAD WORKSPACE ####
#load("~/Git/cuckoos/PA corine extraction.RData")
#load("~/cuckoos/PA corine temp precip extraction.RData") # load workspace on cluster
load("~/Git/cuckoos/PA corine temp precip extraction.RData") # load workspace on PC


#### TO RUN ON CLUSTER, CHANGE BELOW LINE ####
cluster <- FALSE

#### TO RUN RANDOM MCPs, CHANGE BELOW LINE ####
generaterandom <- TRUE
#randomradius <- c(50000, 100000, 200000, 500000) # change to desired search radius for randomized pseudoabsence polygon, in metres
randomradius <- c(50000) # change to desired search radius for randomized pseudoabsence polygon, in metres

#### TO CHANGE NUMBER OF PSEUDO-ABSENCES GENERATED, CHANGE BELOW LINE ####
n.absences <- 10

#### TO CHANGE NUMBER OF INDIVIDUALS, CHANGE BELOW LINE ####
n.birds <- 31

library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
library(rgeos)
library(raster)
library(geosphere)
library(shape)
library(adehabitatHR)
library(SPEI)

allPAobjects <- which(ls() == "PA.epsg4326" | ls() == "PA.epsg3035" | ls() == "r" | ls() == "corine.crs" | ls() == "temp1980.2013" | ls() == "precip1980.2013")

#allPAobjects <- which(ls() == "PA.subset" | ls() == "r" | ls() == "Europeraster" | ls() == "corine.crs")

toremove <- ls()[-allPAobjects]

rm(list=toremove)

#################################################
#
####      SET WORKING DIRECTORIES     ####
#
#################################################

if (!cluster){
  parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
}

if (cluster){
  parentwd <- c("/users1/samf/cuckoos")
}


datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

resampwd <- ("/resampled data + extra grouping variables")
origwd <- ("/original data + extra grouping variables")

corinewd <- paste(parentwd, "/data/corine raster 100x100/v17", sep="")
elevationwd <- paste(parentwd, "/data/elevation data", sep="")

if (!cluster) {
  GISwd <- c("C:/Users/samf/Documents/GIS/cuckoos/protected areas")
}

if (cluster) {
  GISwd <- c("/users1/samf/cuckoos")
}

corine.crs <- CRS("+init=epsg:3035")
epsg4326 <- CRS("+init=epsg:4326")

setwd(parentwd)

#################################################
#
####      LOAD BASE LAYERS (if needed)     ####
#
#################################################


###--- projection and datum information for all data layers ---###

# CORINE => epsg: 3035, LAEA CRS
# Protected Areas (PA.subset) => epsg: 3035
# Climate => epsg: 4326, +proj=latlong +ellps=wGS84 +datum=WGS84
# Elevation => 
# epsg: 3395 - +proj=merc +ellps=wGS84 +datum=WGS84
# epsg: 4326 - +proj=latlong +ellps=wGS84 +datum=WGS84
# epsg: 3035 - LAEA CRS, used for the Corine Land Cover 2006 data

### -------- CLIMATE (runs source file) --------

# 25 Jun 2014: NEW CLIMATE DATA, for calculating SPEI
# 2 raster objects, 1 for total monthly precipitation, 1 for mean temperature for all cells and all months between 1980-2013 (each layer in the rasterstack is a month)
# precip1980.2013, temp1980.2013
# climate data is the same every time, so needs only be run once at beginning of script


if (!("precip1980.2013" %in% ls())) {
  climvarname <- "precipitation"
  setwd(paste(parentwd, "/scripts", sep=""))
  source("extract_netCDF_climate_data_spei.R")
  
  
  climvarname <- "temperature"
  setwd(paste(parentwd, "/scripts", sep=""))
  source("extract_netCDF_climate_data_spei.R")
}

### ------- CORINE land cover raster file -------

if (!("r" %in% ls())) { ### run prep code to load base layers if not in workspace already
  
  setwd(corinewd)
  r <- raster("g100_06.tif")
  
  #Europebox <- extent(2483500,5890375,1276625,4286750) # bounding box extent for relevant parts of Europe (no Scandanavia or Turkey)
  #Europeraster <- crop(r,Europebox)
  
}

### ------- ELEVATION -------

setwd(elevationwd)
elevation <- raster("Europe elevation raster.tif")

### ------- PROTECTED AREAS (vector shapefile) -------

if (!("PA.epsg3035" %in% ls())) { ### run prep code to load base layers if not in workspace already
  
  #   setwd(GISwd)
  #   
  #   PA <- raster("Europe PA raster 50m x 50m.tif")
  #   
  #   colors <- c("white",rep("blue",89662))
  #   PA@legend@colortable <- colors
  #   # default of rpoly colortable is logical(0)
  
  #PAs <- readOGR(GISwd, "terrestrial PAs mainland W Europe corine countries only EPSG 3035")
  PA.epsg4326 <- readOGR(GISwd, "terrestrial_PAs_Europe_GCS_WGS_1984_final")
  
#   PA.subset <- subset(PAs, desig_eng!="UNESCO-MAB Biosphere Reserve" & rep_area!=0)
#   #PA.subset <- subset(PAs, desig_eng!="UNESCO-MAB Biosphere Reserve")
#   PA.subset@data <- droplevels(PA.subset@data)
  
  PA.epsg3035 <- spTransform(PA.epsg4326, CRS=corine.crs)
  
}

##########################################################################################
######              EXTRACT CORINE, PA, ELEVATION, and CLIMATE DATA                 ######
##########################################################################################

### sets the number of bootstraps if random stopovers need generating

# counts the radius with at which absences are generated (50,100,200,500km)
if (generaterandom) {
  numberrepeats <- length(randomradius)
}

# if running presence polygons, then both loops only run a single time
if (!generaterandom) {
  numberrepeats <- 1
  n.absences <- 1
}

#### for loop which counts levels of randomradius ####

for (z in 1:numberrepeats) {  # counts through different levels of randomradius when generaterandom=TRUE
  
  mgroupyear.combn <- list()
  #climate.values <- list()
  #fullclimate <- list()
  #corine.values <- list()
  #elevation.values <- list()
  #fullelevation <- list()
  #MCPoverlap <- list()
  #PAoverlap.info <- list()
  randomstoplist <- list()
  allstoplist <- list()
  
  #=============================================================================================#
  #==== SET-UP OPEN FILE CONNECTIONS TO WRITE EXTRACTED DATA - ALL BIRD FILES  ====
  #=============================================================================================#
  
  # close file connection at end of z in 1:numberrepeats loop
  
  if (generaterandom) { # for random stopovers
    
    ### CLIMATE (all birds together)    
    setwd(paste(datawd,"/corine PA elevation spei extracted values/climate/",sep=""))
    climatefile <- file(paste("climate random ", randomradius[z]/1000, " km radius all birds.csv", sep=""), "w")
    
    ### ELEVATION (all birds together)
    setwd(paste(datawd,"/corine PA elevation spei extracted values/elevation",sep=""))
    elevationfile <- file(paste("elevation random ", randomradius[z]/1000, " km radius all birds.csv", sep=""), "w")    
    
    ### PROTECTED AREAS 
    setwd(paste(datawd,"/corine PA elevation spei extracted values/protected areas",sep=""))
    PAdetailsfile <- file(paste("PA site details random ", randomradius[z]/1000, " km radius all birds.csv", sep=""), "w")    
    PAoverlapfile <- file(paste("stopover overlap area random ", randomradius[z]/1000, " km radius all birds.csv", sep=""), "w")
    
  }
  
  if (!generaterandom) { # for actual stopovers
    
    ### CLIMATE (all birds together)    
    setwd(paste(datawd,"/corine PA elevation spei extracted values/climate/",sep=""))
    climatefile <- file(paste("climate all birds.csv", sep=""), "w")
    
    ### ELEVATION (all birds together)
    setwd(paste(datawd,"/corine PA elevation spei extracted values/elevation",sep=""))
    elevationfile <- file(paste("elevation all birds.csv", sep=""), "w")
    
    ### PROTECTED AREAS
    setwd(paste(datawd,"/corine PA elevation spei extracted values/protected areas",sep=""))
    PAdetailsfile <- file(paste("PA site details all birds.csv", sep=""), "w")    
    PAoverlapfile <- file(paste("stopover overlap area all birds.csv", sep=""), "w")
    
  }
  
  #######################################################
  #######################################################
  #######################################################
  
  #### for loop for each cuckoo file starts here ####  
  
  for (a in 1:n.birds) { # START LOOP
    
    ###=====================================================================###
    #====                  SET-UP DATA & STOPOVER POLYGONS                ====
    ###=====================================================================###
    
    ###------------------ LOAD DATA, SUBSET, ADD PROJECTION INFO, TRANSFORM ------------------
    
    ### LOAD DATA & SUBSET
    
    setwd(paste(datawd,resampwd,sep=""))    
    dataset <- read.csv(list.files()[a], header=T)
    
    ### check that dataset is not an individual to exclude (Idemili & Karma); break out of current loop run and continue to next loop level
    if (dataset$name[1] == "Karma" | dataset$name[1] =="Idemili") {
      next
    }
    
    # convert mgroup to factor
    dataset$mgroup <- as.factor(dataset$mgroup)
    
    # create new concatenated variable, name.mgroup
    #name.mgroup <- as.factor(paste(dataset$name,dataset$mgroup,sep=""))
    #dataset2 <- data.frame(name.mgroup,dataset)
    dataset2 <- data.frame(dataset)
    
    # subset dataset by stopover sites (no breeding sites)
    dataset3 <- subset(dataset2, stopoversite=="Y")
    dataset3 <- droplevels(dataset3)
    
    # check if subsetted dataset has no stopovers; if no stopovers, break out of current count of "a" for loop and skip to next "a"
    if (nrow(dataset3) == 0) {
      next
    }
    
    # summary table of additional grouping variables for combinations of year and mgroup
    mgroupyear.combn[[a]] <- unique(dataset3[c("name", "mgroup", "year", "laststop", "strategy", "Sahara.success")])
    
    ### PROJECT & TRANSFORM    
    # long/lat (original) are in CRS("+init=epsg:4326")
    # long/lat (resampled) are in CRS("+init=epsg:3395")  
    # newlongs/newlats (resampled) are in CRS("+init=epsg:4326")
    # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
    
    coordinates(dataset3) <- c("newlongs","newlats")
    proj4string(dataset3) <- CRS("+init=epsg:4326")
    newdat <- spTransform(dataset3, CRS=corine.crs)
    
    ###------------------- CALCULATE BUFFERED POINTS USED FOR RESAMPLED STOPOVER POLYGONS -------------------
    
    mgroupSPs <- list()
    
    for (i in 1:length(levels(newdat$mgroup))) {
      datsub <- subset(newdat, mgroup==levels(newdat@data$mgroup)[i])
      datsub@data <- droplevels(datsub@data)
      mgroupSPs[[i]] <- gBuffer(datsub, width=500, id=levels(newdat@data$mgroup)[i])
      # mgroupSPs[[i]] <- gBuffer(datsub, width=500, id=levels(datsub@data$name.mgroup))
    }
    
    allstops <- do.call(rbind, mgroupSPs)
    #name.mgroupid <- sapply(slot(allstops, 'polygons'), function(i) slot(i, 'ID'))
    mgroupid <- sapply(slot(allstops, 'polygons'), function(i) slot(i, 'ID'))
    #allstopSPDFs <- SpatialPolygonsDataFrame(allstops, data.frame(name=levels(newdat@data$name), mgroup=levels(newdat@data$mgroup), name.mgroup=name.mgroupid, row.names=name.mgroupid))
    #allstopSPDFs <- SpatialPolygonsDataFrame(allstops, data.frame(mgroupyear.combn[[a]], name.mgroup=name.mgroupid, row.names=name.mgroupid))
    allstopSPDFs <- SpatialPolygonsDataFrame(allstops, data.frame(mgroupyear.combn[[a]], row.names=mgroupid))
    
    
    ###=====================================================================###
    #====                        DATA EXTRACTIONS                          ====
    ###=====================================================================###
    
    #######################################################
    #######################################################
    
    ###============= FUNCTIONS ==================
    
    ###--- FUNCTION to extract CLIMATE DATA with specified arguments ---###    
    extractwithargs <- function(rastername,SPDF){
      x <- extract(rastername, coordinates(SPDF), small=TRUE, buffer=buffer.metres[countinner])
      return(x)
    }
    
    ###--- FUNCTION TO EXTRACT CORINE DATA ---###
    
    newextract <- function(rname, spdf) {
      # split SPDF object into separate SPs by mgroup
      SPlist <- list()
      for (i in 1:nrow(spdf@data)){
        SPlist[[i]] <- subset(spdf, spdf@data$mgroup==levels(spdf@data$mgroup)[i])
      }
      
      # create a list of extents for each mgroup's SP
      extentlist <- lapply(SPlist, extent)
      
      # crop raster according to each mgroup's extent
      rlist <- list()
      for (i in 1:length(extentlist)) {
        rlist[[i]] <- crop(r, extentlist[[i]])
      }
      
      # for each mgroup's SP, extract values from cropped raster        
      extracted.values <- list()
      extracted.values.unlisted <- list()
      for (i in 1:length(SPlist)) {
        extracted.values[[i]] <- extract(rlist[[i]], SPlist[[i]])
        extracted.values.unlisted[[i]] <- unlist(extracted.values[[i]])
      }
      
      return(extracted.values.unlisted)        
    }
    
    #######################################################
    #######################################################
    
    ###============= SET-UP OPEN FILE CONNECTIONS TO WRITE EXTRACTED CORINE DATA - BY INDIVIDUAL BIRD ==================
    
    # file connections are closed at end of (a in 1:31) loop (counts through birds)
    
    if (generaterandom) {
      ### CORINE (one bird per file)      
      setwd(paste(datawd,"/corine PA elevation spei extracted values/corine random",sep=""))
      setwd(paste(randomradius[z]/1000, " km radius", sep=""))
      corinefile <- file(paste(levels(newdat@data$name), " corine.csv", sep=""), "w")
    }
    
    if (!generaterandom) {
      ### CORINE (one bird per file)      
      setwd(paste(datawd,"/corine PA elevation spei extracted values/corine",sep=""))
      corinefile <- file(paste(levels(newdat@data$name), ".csv", sep=""), "w")
    }
    
    #######################################################
    #######################################################
    
    ###============================= LOOP THAT GENERATES n.absences NEW RANDOM STOPOVERS PER BIRD ===============================
    # if generaterandom = FALSE, then n.absences = 1
    
    for (b in 1:n.absences) { # START loop that bootstraps absence polygons
      
      ###---------------------------------------###
      #----   CLIMATE AND CORINE  EXTRACTION ----
      ###---------------------------------------###
      
      ### CLIMATE: Coarse resolution of climate raster means that many coastal areas aren't well covered, stopover polygons extract NANs at the 25km buffer level; for stopovers that extract NAs from the climate raster, increase the buffer radius by 25km steps until it captures cells with values for the problem stopovers.  i.e. extract the mean climate values from within a 25km buffer of the stopover centroid, and if no climate data within 25km of the stopover, then increase the buffer until climate data is available and extract it
      
      ### CORINE: Random polygon placement means that some random stopovers may be placed in areas which either have poor terrestrial coverage (> 50% water) or poor coverage by the Corine raster. For generaterandom=TRUE, want to check that random stopovers pick up corine (as well as climate) raster cells 
      
      ### SPEI
      # 6-month SPEI values are calculated based on the current month's plus 5 previous months data
      # this timescale should provide a good representation of vegetation growing conditions over time periods likely to be relevant for determining the quality of habitats that cuckoos use
      # SPEI-month values to use will be March and August
      # SPEI-Mar represents the conditions over the winter (Oct-Mar period), which is of greatest relevance in the Mediterranean where winter rainfall is most important for vegetation growth
      # SPEI-Aug represents the conditions over the spring and summer (Mar-Aug period), which is of greater relevance in the more northern parts of Europe where rainfall early in the growing season sets the conditions for vegetation growth
      # from Vicente-Serrano et al. (2012 PNAS) a 6 month SPEI may be the best trade-off in using a timescale that will be relatively well correlated with vegetation activity (as measured by NDVI and other indices)
      
      # doing extract() on the climate data rasterstack returns a matrix where columns correspond to month/years and rows correspond to a grid cell (buffer on coordinates will make some extractions come up with multiple grid cells), sometimes dozens if the stopover centroid has a very large buffer on it because of extracting NAs from the climate raster
      
      ## OPTIONS for SPEI ##
      # 1) for each row (grid cell) of extracted climate data, add to climate dataframe(year, month, precip, tmean) then run the SPEI calculation, then average the SPEIs from all grid cells in a given month/year
      # OR 2)****** average extracted data across all grid cells for a given month/year, then run the SPEI calculation on the averaged data (MMT = mean monthly totals) 
      
      
      ###----------------------- CLIMATE/CORINE EXTRACTION, IF GENERATERANDOM = TRUE --------------------------------
      
      # climate raster is coarser resolution than corine, so more likely to generate NAs for cell in which stopover polygon falls
      # run source file to generate random buffered polygons and extract climate data such that all stopovers have associated climate data at the smallest buffer level that generates non-NA values
      
      if (generaterandom) { # START generaterandom = TRUE loop
        
        repeat { # OUTER repeat loop, for corine extractions > 50% marine habitats (open water & wetlands)
          
          ####---- CLIMATE EXTRACTION (extracted inside source) ----####
          setwd(paste(parentwd, "/scripts", sep=""))
          source("cuckoo_generate_random_buffered_polygons_spei.R")
          #     randomtry[a] <- countouter
          
          ####---- SPEI calculation ----####
          MMTprecip <- lapply(lapply(useprecip, function(x) {apply(x,2,mean,na.rm=TRUE)}), function(x) {data.frame(yearmonth=names(x),precip=x)})
          MMTtemp <- lapply(lapply(usetemp, function(x) {apply(x,2,mean,na.rm=TRUE)}), function(x) {data.frame(yearmonth=names(x),temp=x)})
          
          climdat <- list()
          climdat.SPEI <- list()
          
          #spei.data <- randomstop4326@data
          
          spei.Mar <- rep(NA, nrow(randomstop4326@data))
          spei.Aug <- rep(NA, nrow(randomstop4326@data))
          
          spei.data <- data.frame(randomstop4326@data, spei.Mar, spei.Aug)
          
          for (i in 1:length(MMTprecip)) {
            
            precip <- MMTprecip[[i]]$precip
            temp <- MMTtemp[[i]]$temp
            climdat[[i]] <- data.frame(year=as.integer(substring(MMTprecip[[i]]$yearmonth,2,5)), month=as.factor(substring(MMTprecip[[i]]$yearmonth,6)),precip, temp)
            climdat[[i]]$PET <- thornthwaite(climdat[[i]]$temp, coordinates(subset(randomstop4326, mgroup %in% names(MMTprecip)[i]))[2])
            spei6 <- spei(climdat[[i]]$precip-climdat[[i]]$PET, 6, na.rm=TRUE)
            climdat.SPEI[[i]] <- data.frame(climdat[[i]], spei6=as.numeric(spei6$fitted))
            
            spei.data$spei.Mar[which(spei.data$mgroup %in% names(MMTprecip)[i])] <- climdat.SPEI[[i]][climdat.SPEI[[i]]$year==spei.data[spei.data$mgroup==names(MMTprecip)[i],"year"] & climdat.SPEI[[i]]$month=="Mar","spei6"]
            spei.data$spei.Aug[which(spei.data$mgroup %in% names(MMTprecip)[i])] <- climdat.SPEI[[i]][climdat.SPEI[[i]]$year==spei.data[spei.data$mgroup==names(MMTprecip)[i],"year"] & climdat.SPEI[[i]]$month=="Aug","spei6"]            
            #spei.data$spei.Mar[which(spei.data$mgroup %in% names(MMTprecip)[i])] <- spei.Mar
            #spei.data$spei.Aug[which(spei.data$mgroup %in% names(MMTprecip)[i])] <- spei.Aug      
            
          }          
          
          #######################################################
          
          #--- CORINE ---#
          
          corine.values <- newextract(r, randomstopSPDFs)  # formerly corine.values[[a]]
          
          check <- logical()
          
          # corine categories to check for in > 0.5 of extracted values are any marine categories (37 = Salt marshes, 38 = Salines, 39 =  Intertidal flats, 42 = Coastal lagoons, 43 = Estuaries, 44 = Sea and ocean, as well as NAs)
          for (i in 1:length(corine.values)) {
            if ( (length(which(corine.values[[i]]==37 | corine.values[[i]]==38 | corine.values[[i]]==39 | corine.values[[i]]==41 | corine.values[[i]]==42 | corine.values[[i]]==43 | corine.values[[i]]==44 | is.na(corine.values[[i]])))/length(corine.values[[i]]) < 0.5) ) {
              check[i] <- TRUE
            } else {
              check[i] <- FALSE
            }
          }
          
          # checks whether all random stopovers have < 50% marine waters (coastal lagoons, estuaries, sea & ocean) and wetlands
          # if ANY random stopovers come up with > 50% marine habitat, then ALL stopovers are randomly generated again (rather than just the problem stopovers) - too hard to regenerate just the problem one!
          if (all(check==TRUE)) break # breaks OUTER repeat loop, for corine extractions
          
        } # END outer repeat loop, for corine extractions
        
        corine.MCP <- list()
        
        for (i in 1:length(corine.values)) {
          corine.MCP[[i]] <- data.frame(
            name=levels(newdat@data$name), 
            mgroup=levels(newdat@data$mgroup)[i], 
            year=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "year"][1], 
            laststop=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "laststop"][1],
            strategy=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "strategy"][1],
            Sahara.success=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "Sahara.success"][1],
            corine.values=corine.values[[i]])
        }
        
        ### WRITE CLIMATE SPEI DATA TO FILE
        if (a==1 & b==1) {
          write.table(data.frame(spei.data, nabsence=b), file=climatefile, row.names=FALSE, sep=",")
        } else {
          write.table(data.frame(spei.data, nabsence=b), file=climatefile, row.names=FALSE, col.names=FALSE, sep=",")
        }
        
        ### WRITE CORINE DATA TO FILE
        fullcorine <- do.call(rbind,corine.MCP) # unlists random polygons
        
        if (b==1) {
          write.table(data.frame(fullcorine, nabsence=b), file=corinefile, row.names=FALSE, sep=",")
        } else {
          write.table(data.frame(fullcorine, nabsence=b), file=corinefile, row.names=FALSE, col.names=FALSE, sep=",")
        }
        
        
      } # END generaterandom = TRUE loop
      
      
      #######################################################
      #######################################################
      
      ###----------------------- CLIMATE/CORINE EXTRACTION, IF GENERATERANDOM = FALSE --------------------------------
      
      if (!generaterandom) { # START generaterandom = FALSE loop
        
        #==== CLIMATE ====
        
        allstop4326 <- spTransform(allstopSPDFs, CRS=epsg4326)
        
        #==== Extract precip & temp data for all year/month climate rasters for each stopover ====
        
        buffer.metres <- c(25000,50000,75000,100000,125000,150000) # different than the random stopovers code, which only buffers up to 50km.  Less leeway with the actual stopovers, so buffer up to 150km in the hopes of finding a non-NAN raster cell, if necessary
        
        countinner <- 1 # counts the buffer levels, starts at smallest buffer (25 km)
        
        extractprecip <- extractwithargs(precip1980.2013, allstop4326)
        names(extractprecip) <- allstop4326@data$mgroup
        extracttemp <- extractwithargs(temp1980.2013, allstop4326)
        names(extracttemp) <- allstop4326@data$mgroup
        
        problems <- do.call(rbind,lapply(extractprecip, function(x) {all(is.na(x))}))
        probmgroupids <- allstop4326@data[which(problems==TRUE), "mgroup"] # returns mgroupid of stopover which produces ALL NAs at current buffer level
        
        if (all(problems==FALSE)) { # if there are no problem stopovers that extract all NAs
          
          useprecip <- extractprecip
          usetemp <- extracttemp
          
        } else {
          
          okmgroups <- allstop4326@data[which(problems==FALSE), "mgroup"]
          useprecip <- extractprecip[which(names(extractprecip) %in% as.character(okmgroups))]
          usetemp <- extracttemp[which(names(extractprecip) %in% as.character(okmgroups))]
          
          repeat {        
            
            probstops <- subset(allstop4326, mgroup %in% probmgroupids)
            countinner <- countinner + 1
            
            prob.extractprecip <- extractwithargs(precip1980.2013, probstops)
            names(prob.extractprecip) <- probstops@data$mgroup
            prob.extracttemp <- extractwithargs(temp1980.2013, probstops)
            names(prob.extracttemp) <- probstops@data$mgroup
            
            problems <- do.call(rbind,lapply(prob.extractprecip, function(x) {all(is.na(x))}))      
            
            # if all these problem mgroups now have values (ie. all are non-NA), then break out of the repeat loop and append extracted data to the list of data from the ok mgroups
            if (all(problems==FALSE)) {
              useprecip <- c(useprecip,prob.extractprecip)
              usetemp <- c(usetemp,prob.extracttemp)
              break
              
            } else {
              
              # if there is still > 1 mgroup with a problem, then need to continue repeat loop; first, check if any mgroups that were a problem are now ok (ie. are non-NA) at the new buffer level. For ones that are ok, append data to the list of data from ok mgroups (created at the top of the "else" statement), and repeat the loop with the mgroups that are still a problem
              okmgroups <- probstops@data[which(problems==FALSE), "mgroup"]
              rerunmgroups.precip <- prob.extractprecip[which(names(prob.extractprecip) %in% as.character(okmgroups))]
              rerunmgroups.temp <- prob.extracttemp[which(names(prob.extractprecip) %in% as.character(okmgroups))]
              
              useprecip <- c(useprecip, rerunmgroups.precip)
              usetemp <- c(usetemp, rerunmgroups.temp)
              
              probmgroupids <- probstops@data[which(problems==TRUE), "mgroup"]
              
            } # end inner IF statement (if problem mgroups are all still NA)
          } # end REPEAT loop - cycles through buffer levels
          
          
        } # end outer IF statement (if some mgroups are NA)
        
        #==== SPEI calculation  ====        
        MMTprecip <- lapply(lapply(useprecip, function(x) {apply(x,2,mean,na.rm=TRUE)}), function(x) {data.frame(yearmonth=names(x),precip=x)})
        MMTtemp <- lapply(lapply(usetemp, function(x) {apply(x,2,mean,na.rm=TRUE)}), function(x) {data.frame(yearmonth=names(x),temp=x)})
        
        climdat <- list()
        climdat.SPEI <- list()
        
        spei.Mar <- rep(NA, nrow(allstop4326@data))
        spei.Aug <- rep(NA, nrow(allstop4326@data))
        
        spei.data <- data.frame(allstop4326@data, spei.Mar, spei.Aug)
        
        for (i in 1:length(MMTprecip)) {
          
          precip <- MMTprecip[[i]]$precip
          temp <- MMTtemp[[i]]$temp
          climdat[[i]] <- data.frame(year=as.integer(substring(MMTprecip[[i]]$yearmonth,2,5)), month=as.factor(substring(MMTprecip[[i]]$yearmonth,6)),precip, temp)
          climdat[[i]]$PET <- thornthwaite(climdat[[i]]$temp, coordinates(subset(allstop4326, mgroup %in% names(MMTprecip)[i]))[2])
          spei6 <- spei(climdat[[i]]$precip-climdat[[i]]$PET, 6, na.rm=TRUE)
          climdat.SPEI[[i]] <- data.frame(climdat[[i]], spei6=as.numeric(spei6$fitted))
          
          spei.data$spei.Mar[which(spei.data$mgroup %in% names(MMTprecip)[i])] <- climdat.SPEI[[i]][climdat.SPEI[[i]]$year==spei.data[spei.data$mgroup==names(MMTprecip)[i],"year"] & climdat.SPEI[[i]]$month=="Mar","spei6"]
          spei.data$spei.Aug[which(spei.data$mgroup %in% names(MMTprecip)[i])] <- climdat.SPEI[[i]][climdat.SPEI[[i]]$year==spei.data[spei.data$mgroup==names(MMTprecip)[i],"year"] & climdat.SPEI[[i]]$month=="Aug","spei6"]            
          #spei.data$spei.Mar[which(spei.data$mgroup %in% names(MMTprecip)[i])] <- spei.Mar
          #spei.data$spei.Aug[which(spei.data$mgroup %in% names(MMTprecip)[i])] <- spei.Aug      
          
        }
        
        #       # plot SPEI
        #       par(mfrow=c(2,2))
        #       for (i in 1:length(MMTprecip)) {
        #         
        #         precip <- MMTprecip[[i]]$precip
        #         temp <- MMTtemp[[i]]$temp
        #         climdat[[i]] <- data.frame(year=as.integer(substring(MMTprecip[[i]]$yearmonth,2,5)), month=as.factor(substring(MMTprecip[[i]]$yearmonth,6)),precip, temp)
        #         climdat[[i]]$PET <- thornthwaite(climdat[[i]]$temp, coordinates(subset(allstop4326, name.mgroup %in% names(MMTprecip)[i]))[2])
        #         #spei6 <- spei(ts(climdat[[i]]$precip-climdat[[i]]$PET, 6))
        #         #plot(spei6)
        #         
        #         plot(spei(ts(climdat[[i]]$precip-climdat[[i]]$PET,freq=12,start=c(1980,1)),12))
        #         mtext(names(MMTprecip)[i], side=4)
        #         
        #       }
        
        ### CLIMATE SUMMARY, ACTUAL STOPOVERS ###
        # 18 birds captured non-NAN cells at 25km, 5 at 50km, 4 at 75km, 1 at 100km, 1 at 150km
        # See per bird information summary file at C:\Users\samf\Documents\Git\cuckoos\output\all cuckoo stopovers and non-NAN climate raster buffer values
        
        #######################################################
        
        #--- CORINE ---#
        
        corine.values <- newextract(r, allstopSPDFs)
        
        corine.MCP <- list()
        
        for (i in 1:length(corine.values)) {
          corine.MCP[[i]] <- data.frame(
            name=levels(newdat@data$name), 
            mgroup=levels(newdat@data$mgroup)[i], 
            year=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "year"][1], 
            laststop=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "laststop"][1],
            strategy=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "strategy"][1],
            Sahara.success=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "Sahara.success"][1],
            corine.values=corine.values[[i]])
        }
        
        ### WRITE CLIMATE SPEI DATA TO FILE
        if (a==1 & b==1) {
          write.table(spei.data, file=climatefile, row.names=FALSE, sep=",")
        } else {
          write.table(spei.data, file=climatefile, row.names=FALSE, col.names=FALSE, sep=",")
        }
        
        ### WRITE CORINE DATA TO FILE
        fullcorine <- do.call(rbind,corine.MCP) # unlists random polygons
        
        if (b==1) {
          write.table(fullcorine, file=corinefile, row.names=FALSE, sep=",")
        } else {
          write.table(fullcorine, file=corinefile, row.names=FALSE, col.names=FALSE, sep=",")
        }
        
        
      } # end IF !generaterandom
      
      
      #######################################################
      #######################################################
      #######################################################
      
      ###-------- REAL OR RANDOM STOPS --------
      
      if (generaterandom) {
        mcp <- randomstopSPDFs
        
        ####============ CREATE NEW SPDF with unique row.ids for bootstrapped absence polygon
        randomstopSPDFs@data <- data.frame(randomstopSPDFs@data, scale=randomradius[z]/1000, nabsence=b, name.mgroup.nbootstrap=paste(randomstopSPDFs@data$name,randomstopSPDFs@data$mgroup,"_",b, sep=""))
        
        # from http://gis.stackexchange.com/questions/37503/rename-a-spatialpolygon-class-object-in-r
        # DISSOLVE BASED ON "REGION" COLUMN USING rgeos::gUnionCascaded 
        randompolygons <- gUnionCascaded(randomstopSPDFs, id=randomstopSPDFs@data$name.mgroup.nbootstrap)
        
        # CREATE A DATAFRAME OF VALUES RETAINED FROM "REGION" COLUMN
        sdf <- data.frame(ID=row.names(randompolygons))
        
        # ASSIGN "REGION" VALUES TO rownames SO DATA MATCHES slots CORRECTLY
        row.names(sdf) <- row.names(randompolygons)
        
        # CREATE A SpatialPolygonsDataFrame OBJECT WITH A data slot HOLDING REGION IDS.  
        randompolygons <-  SpatialPolygonsDataFrame(randompolygons, sdf)
        randomstoplist[[b]] <- randompolygons
        
        # ADD rest of original dataframe
        randomstoplist[[b]]@data <- data.frame(randompolygons@data, randomstopSPDFs@data[,c("name","mgroup","year","laststop","strategy","Sahara.success","scale","nabsence")])
        names(randomstoplist[[b]]@data) <- c("ID","name","mgroup","year","last","strat","Sahsuc","scale","nboot")
      }
      
      if (!generaterandom) {
        mcp <- allstopSPDFs
        allstopSPDFs@data <- data.frame(allstopSPDFs@data, name.mgroup=paste(allstopSPDFs@data$name, allstopSPDFs@data$mgroup, sep=""))
        
        # from http://gis.stackexchange.com/questions/37503/rename-a-spatialpolygon-class-object-in-r
        # DISSOLVE BASED ON "REGION" COLUMN USING rgeos::gUnionCascaded 
        truepolygons <- gUnionCascaded(allstopSPDFs, id=allstopSPDFs@data$name.mgroup)
        
        # CREATE A DATAFRAME OF VALUES RETAINED FROM "REGION" COLUMN
        sdf <- data.frame(ID=row.names(truepolygons))
        
        # ASSIGN "REGION" VALUES TO rownames SO DATA MATCHES slots CORRECTLY
        row.names(sdf) <- row.names(truepolygons)
        
        # CREATE A SpatialPolygonsDataFrame OBJECT WITH A data slot HOLDING REGION IDS.  
        truepolygons <-  SpatialPolygonsDataFrame(truepolygons, sdf)
        allstoplist[[a]] <- truepolygons
        
        # ADD rest of original dataframe
        allstoplist[[a]]@data <- data.frame(truepolygons@data, allstopSPDFs@data[,c("name","mgroup","year","laststop","strategy","Sahara.success")])
        names(allstoplist[[a]]@data) <- c("ID","name","mgroup","year","last","strat","Sahsuc")
      }
      
      ###-------------------------------###
      
      #######################################################
      #######################################################
      #######################################################
      
      ###---------------------------###
      #----   ELEVATION  ----
      ###---------------------------###
      
      mcp.elevtransform <- spTransform(mcp, CRS=elevation@crs)
      
      elevation.values <- extract(elevation, mcp.elevtransform)    # formerly elevation.values[[a]]
      
      elevation.dataset <- list()
      
      for (i in 1:length(elevation.values)) {
        elevation.dataset[[i]] <- data.frame(
          name=levels(newdat@data$name), 
          mgroup=levels(newdat@data$mgroup)[i], 
          year=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "year"][1], 
          laststop=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "laststop"][1],
          strategy=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "strategy"][1],
          Sahara.success=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "Sahara.success"][1],
          elevation.m=elevation.values[[i]])
      }
      
      fullelevation <- do.call(rbind, elevation.dataset)
      
      ### WRITE ELEVATION DATA TO FILE
      
      if (generaterandom) {
        if (a==1 & b==1) {
          write.table(data.frame(fullelevation, nabsence=b), file=elevationfile, row.names=FALSE, sep=",")
        } else {
          write.table(data.frame(fullelevation, nabsence=b), file=elevationfile, row.names=FALSE, col.names=FALSE, sep=",")
        }
      }
      
      if (!generaterandom) {
        if (a==1 & b==1) {
          write.table(fullelevation, file=elevationfile, row.names=FALSE, sep=",")
        } else {
          write.table(fullelevation, file=elevationfile, row.names=FALSE, col.names=FALSE, sep=",")
        }
      }
      
      #######################################################
      #######################################################
      #######################################################
      
      ###---------------------------###
      #----   PROTECTED AREA  ----
      ###---------------------------###
      
      mcp.id <- sapply(slot(mcp, "polygons"), function(i) slot(i, "ID"))
      mcp.area <- gArea(mcp, byid=TRUE)/1000000
      
      # subset PA dataset to only those which intersect with MCPs to make gIntersection calculation faster
      #subsetintersect <- gIntersects(mcp, PA.subset, byid=TRUE, prepared=TRUE)
      subsetintersect <- gIntersects(mcp, PA.subset, byid=TRUE)
      newsubinter <- as.data.frame(subsetintersect)
      
      # for each column (MCP id) of "subsetintersect", extract row number (PA id) that is TRUE (ie. which PAs the MCPs intersect with)
      # compile into vector of id row numbers to use to subset the PA.subset shapefile to only the intersecting PAs
      ids <- numeric()
      for (i in 1:length(mcp)){
        add <- which(newsubinter[,i] == TRUE)
        ids <- c(ids, add)
      }
      
      ### check if any MCPs intersect with PAs, if not, then output 0's to MCPoverlap variable, otherwise, extract overlap info ###
      if (length(ids) == 0) {
        MCPoverlap <- data.frame(    # formerly MCPoverlap[[a]]
          name=levels(newdat@data$name), 
          MCPid=mcp.id,
          MCP.area.km2=mcp.area,
          PA.overlap.area.km2=rep(0, length(mcp)), 
          prop.MCP.overlap=rep(0, length(mcp))
        )
      } else {
        
        # subset the PA.subset SPDF to only the intersecting polygons
        # add unique(ids) because some polygons are extracted twice (presumably if they are more linear features such as rivers and intersect the MCP in two different and distinctly separate places?)
        PA.subsetintersect <- PA.subset[ids,]
        PA.subsetintersect@data <- droplevels(PA.subsetintersect@data)
        
        # derive polygon intersections between MCPs and PAs
        # duplicated polygon intersections sometimes found (not sure why, maybe because of duplicated polygons in the PA shapefile dataset?  in any case, remove duplicates)
        PAintersection <- gIntersection(mcp, PA.subsetintersect, byid=TRUE)
        
        if (any(duplicated(names(PAintersection)))){
          PAintersection2 <- PAintersection[-which(duplicated(names(PAintersection)))]} else {
            PAintersection2 <- PAintersection
          }
        
        # extract parent polygon ids using substrings - names(PAintersection) using unlist(strsplit) by " "
        rowids <- names(PAintersection2)
        splitids <- strsplit(rowids, " ")
        parentMCP <- numeric()
        parentPA <- numeric()
        for (i in 1:length(splitids)){
          parentMCP[i] <- as.numeric(unlist(splitids[[i]])[1])
          parentPA[i] <- as.numeric(unlist(splitids[[i]])[2])
        }
        
        # merge overlapping polygons, convert from m^2 to km^2
        # mergePA polygons will no longer be identified with associated attributes of underlying PAs - to go back and extract specific attributes for underlying PA polygons, use the split string vector parentPA and look for row of interest in PA.subset using as.numeric(row.names(PA.subsetintersect)))
        mergePA <- gUnaryUnion(PAintersection2, id=parentMCP)
        PA.overlap.area <- gArea(mergePA, byid=TRUE)/1000000
        
        # 
        PA.overlap.area2 <- numeric()
        
        PA.overlap.area2[which(mcp.id %in% names(PA.overlap.area))] <- PA.overlap.area
        PA.overlap.area2[which(!(mcp.id %in% names(PA.overlap.area)))] <- 0
        
        # max(gArea(PAintersection, byid=TRUE))/1000000
        
        # proportion area of MCP overlapping with PA (for MCPs which intersect with a PA)
        #MCP.area.km2 <- mcp@data[which(names(PA.overlap.area) %in% mcp$id),c("area")]
        MCP.area.km2 <- mcp.area
        prop.MCP.overlap <- PA.overlap.area2/MCP.area.km2        
        
        MCPoverlap <- data.frame(     # formerly MCPoverlap[[a]]
          name=levels(newdat@data$name), 
          MCPid=mcp.id, 
          MCP.area.km2, 
          PA.overlap.area.km2=PA.overlap.area2, 
          prop.MCP.overlap
        )
        
        ### WRITE PROTECTED AREA OVERLAP DATA TO FILE
        if (generaterandom) {
          if (a==1 & b==1) {
            write.table(data.frame(MCPoverlap, nabsence=b), file=PAoverlapfile, row.names=FALSE, sep=",")
          } else {
            write.table(data.frame(MCPoverlap, nabsence=b), file=PAoverlapfile, row.names=FALSE, col.names=FALSE, sep=",")
          }
        }
        
        if (!generaterandom){
          if (a==1 & b==1) {
            write.table(MCPoverlap, file=PAoverlapfile, row.names=FALSE, sep=",")
          } else {
            write.table(MCPoverlap, file=PAoverlapfile, row.names=FALSE, col.names=FALSE, sep=",")
          }
        }
        
        # attribute info about overlapping PAs
        indPA.overlap <- gArea(PAintersection2, byid=TRUE)/1000000 
        #   names(indPA.overlap) <- names(PAintersection)
        
        
        indPA.overlap2 <- indPA.overlap
        
        PAoverlap.info <- data.frame(     # formerly PAoverlap.info[[a]]
          birdname=levels(newdat@data$name),
          mgroup=parentMCP,
          PA.subsetintersect@data[,c("wdpaid","country","name","desig_eng","desig_type","iucn_cat","rep_area")],
          MCP.overlap.area.km2=indPA.overlap2
        )
        
        ### WRITE PROTECTED AREA SITE ATTRIBUTES DATA TO FILE
        if (generaterandom) {
          if (a==1 & b==1) {
            write.table(data.frame(PAoverlap.info, nabsence=b), file=PAdetailsfile, row.names=FALSE, sep=",")
          } else {
            write.table(data.frame(PAoverlap.info, nabsence=b), file=PAdetailsfile, row.names=FALSE, col.names=FALSE, sep=",")
          }
        }
        
        if (!generaterandom) {
          if (a==1 & b==1) {
            write.table(PAoverlap.info, file=PAdetailsfile, row.names=FALSE, sep=",")
          } else {
            write.table(PAoverlap.info, file=PAdetailsfile, row.names=FALSE, col.names=FALSE, sep=",")
          }
        }
        
      } # end IF statement checking for any PA overlap
      
    } # end FOR loop counting through the number of times a random polygon (pseudo-absence) is generated
    
    ### write an individual's bootstrapped absence polygons to ESRI shapefile (after counting through 1:n.absences)
    if (generaterandom) {
      allrandombootstrappedpolygons <- do.call(rbind, randomstoplist)
      setwd(paste(outputwd,"/bootstrapped absence polygons",sep=""))
      writeOGR(allrandombootstrappedpolygons, dsn=".", layer=paste(levels(newdat@data$name), randomradius[z]/1000, "km random bootstrapped polygons", sep=" "), driver="ESRI Shapefile")
    }   
    
    close(corinefile)
    
  }  # end FOR loop counting through individuals 1:31
  
  close(climatefile)
  close(elevationfile)
  close(PAdetailsfile)
  close(PAoverlapfile)
  
  # 
  
  ### write ALL individuals' presence polygons to ESRI shapefile
  if (!generaterandom) {
    non.null.list <- allstoplist[!sapply(allstoplist, is.null)]
    alltruepolygons <- do.call(rbind, non.null.list)
    setwd(paste(outputwd,"/presence polygons",sep=""))
    writeOGR(alltruepolygons, dsn=".", layer="all birds presence polygons", driver="ESRI Shapefile")   
  }
  
} # end FOR loop counting through numberrepeats (random radius levels if absences are being generated)

# x <- SpatialPoints(coordinates(allstopSPDFs), proj4string=corine.crs)
# poly <- subset(allstopSPDFs, mgroup==5)
# x <- extract(r, poly)
# 
# y <- SpatialPoints(coordinates(poly), proj4string=corine.crs)
