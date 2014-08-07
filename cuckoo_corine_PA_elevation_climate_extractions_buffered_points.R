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
#
###########################################################

################### TO DO #########################

# need to add data extractions to this script for:
#   - elevation data
#   - climate data

###--- TWO APPROACHES FOR PROTECTED AREA ANALYSIS

### APPROACH 1: do cuckoo locations fall within protected areas?
# Extract protected area information from actual cuckoo location points (either original or resampled)

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
load("~/Git/cuckoos/PA corine temp precip extraction.RData")


#### TO RUN ON CLUSTER, CHANGE BELOW LINE ####
cluster <- FALSE

#### TO RUN RANDOM MCPs, CHANGE BELOW LINE ####
generaterandom <- FALSE
randomradius <- c(50000, 100000, 200000, 500000) # change to desired search radius for randomized pseudoabsence polygon, in metres

library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
library(rgeos)
library(raster)
library(geosphere)
library(shape)
library(adehabitatHR)
library(SPEI)

# allPAobjects <- which(ls() == "PA.subset" | ls() == "r" | ls() == "corine.crs" | ls() == "temp1980.2013" | ls() == "precip1980.2013")
# 
# #allPAobjects <- which(ls() == "PA.subset" | ls() == "r" | ls() == "Europeraster" | ls() == "corine.crs")
# 
# toremove <- ls()[-allPAobjects]
# 
# rm(list=toremove)

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
  GISwd <- c("C:/Users/samf/Documents/GIS/cuckoos")
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

# runs source file with output 2 rasterbricks for each climate variable (meanwintertemp, meanwinterprecip, meanspringtemp, meanspringprecip), where winter is Nov-Feb and spring is Mar-May
# run source file for each of precipitation and temperature
# OUTPUT: a list of RasterBricks called climate.data
# climate data is the same every time, so needs only be run once at beginning of script

# 25 Jun 2014: NEW CLIMATE DATA
# 2 raster objects, 1 for total monthly precipitation, 1 for mean temperature for all cells and all months between 1980-2013
# precip1980.2013, temp1980.2013

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

if (!("PA.subset" %in% ls())) { ### run prep code to load base layers if not in workspace already
  
  #   setwd(GISwd)
  #   
  #   PA <- raster("Europe PA raster 50m x 50m.tif")
  #   
  #   colors <- c("white",rep("blue",89662))
  #   PA@legend@colortable <- colors
  #   # default of rpoly colortable is logical(0)
  
  PAs <- readOGR(GISwd, "terrestrial PAs mainland W Europe corine countries only EPSG 3035")
  
  PA.subset <- subset(PAs, desig_eng!="UNESCO-MAB Biosphere Reserve" & rep_area!=0)
  #PA.subset <- subset(PAs, desig_eng!="UNESCO-MAB Biosphere Reserve")
  PA.subset@data <- droplevels(PA.subset@data)
  
  PA.subset <- spTransform(PA.subset, CRS=corine.crs)
  
}

#################################################
#
# EXTRACT CORINE, PA, ELEVATION, and CLIMATE DATA (resampled data only) ####
#
#################################################

### TO RUN ORIGINAL or RESAMPLED data, change below line

original <- FALSE # or FALSE if resampled

# create list to hold datasets with extracted corine and protected area values for each cuckoo

if (generaterandom) {
  numberrepeats <- length(randomradius)
}

if (!generaterandom) {
  numberrepeats <- 1
}

#### for loop which counts levels of randomradius ####

for (z in 1:numberrepeats) {  # counts through different levels of randomradius when generaterandom=TRUE
  
  mgroupyear.combn <- list()
  climate.values <- list()
  fullclimate <- list()
  corine.values <- list()
  fullcorine <- list()
  elevation.values <- list()
  fullelevation <- list()
  MCPoverlap <- list()
  PAoverlap.info <- list()
  
  #######################################################
  #######################################################
  #######################################################
  
  #### for loop for each cuckoo file starts here ####
  
  for (a in 1:31) { # START LOOP 
    
    ###--- LOAD DATA, ADD PROJECTION INFO, TRANSFORM ---###
    
    if (original) setwd(paste(datawd,origwd,sep=""))
    if (!original) setwd(paste(datawd,resampwd,sep=""))
    
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
    
    # long/lat (original) are in CRS("+init=epsg:4326")
    # long/lat (resampled) are in CRS("+init=epsg:3395")  
    # newlongs/newlats (resampled) are in CRS("+init=epsg:4326")
    # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
    
    # -------------------- RESAMPLED DATA ---------------------#
    if (!original) {
      coordinates(dataset3) <- c("newlongs","newlats")
      proj4string(dataset3) <- CRS("+init=epsg:4326")
      newdat <- spTransform(dataset3, CRS=corine.crs)
    }
    
    # -------------------- ORIGINAL DATA ----------------------#
    if (original) {
      coordinates(dataset3) <- c("newlongs","newlats")
      proj4string(dataset3) <- CRS("+init=epsg:4326")
      newdat <- spTransform(dataset3, CRS=corine.crs)
      
    }
    
    # summary table of additional grouping variables for combinations of year and mgroup
    mgroupyear.combn[[a]] <- unique(newdat@data[c("name", "mgroup", "year", "laststop", "strategy", "Sahara.success")])
    
    ###-----------------------------------------------###
    #---  CALCULATE BUFFERED POINTS USED FOR RESAMPLED STOPOVER LOCATIONS ----
    ###-----------------------------------------------###
    
      
    ######## NOT USING THIS KERNEL METHOD ############
    # # for each movement group of each cuckoo, calculate the UD using function kernelUD() in library(adehabitatHR)
    # newdat@data$mgroup <- as.factor(newdat@data$mgroup)
    # stoppoints <- newdat[,c("mgroup")]
    # kud <- kernelUD(stoppoints, h="href")
    # # kud1 <- kernelUD(stoppoints, h="LSCV") # this doesn't work - doesn't minimize
    # homerange <- getverticeshr(kud)
    ##################################################
    
    mgroupSPs <- list()
    
    for (i in 1:length(levels(newdat$mgroup))) {
      datsub <- subset(newdat, mgroup==levels(newdat@data$mgroup)[i])
      datsub@data <- droplevels(datsub@data)
      #mgroupSPs[[i]] <- gBuffer(datsub, width=500, id=levels(newdat@data$mgroup)[i])
      mgroupSPs[[i]] <- gBuffer(datsub, width=500, id=levels(datsub@data$name.mgroup))
    }
    
    allstops <- do.call(rbind, mgroupSPs)
    name.mgroupid <- sapply(slot(allstops, 'polygons'), function(i) slot(i, 'ID'))
    allstopSPDFs <- SpatialPolygonsDataFrame(allstops, data.frame(name=levels(newdat@data$name), mgroup=levels(newdat@data$mgroup), name.mgroup=mgroupid, row.names=mgroupid))
    
#     mgroupSPs <- list()
#     
#     for (i in 1:length(levels(newdat$mgroup))) {
#       newdatsub <- subset(newdat, mgroup==levels(newdat$mgroup)[i])
#       newdatsub@data <- droplevels(newdatsub@data)
#       mgroupSPs[[i]] <- gBuffer(newdatsub, width=500, id=levels(newdat$mgroup)[i])
#     }
#     
#     allstops <- do.call(rbind, mgroupSPs)
#     mgroupid <- sapply(slot(allstops, 'polygons'), function(i) slot(i, 'ID'))
#     allstopSPDFs <- SpatialPolygonsDataFrame(allstops, data.frame(mgroup=mgroupid, row.names=mgroupid))
    
    #   # extract IDs of all stopover polygons
    #   sapply(slot(allstopSPDFs, 'polygons'), function(i) slot(i, 'ID'))
    
    # plot(croppedraster)
    # plot(allstopSPDFs, add=T, lwd=1.5)
    # plot(newdat, add=T, pch=16, cex=0.6)
    
    ###########################
    #
    ####  DATA EXTRACTIONS ####
    #
    ###########################
    
    ###---------------------------###
    #----   CLIMATE AND CORINE  ----
    ###---------------------------###
    
    ### CLIMATE: Coarse resolution of climate raster means that many coastal areas aren't well covered, stopover polygons extract NANs at the 25km buffer level; for stopovers that extract NAs from the climate raster, increase the buffer radius by 25km steps until it captures cells with values for the problem stopovers.  i.e. extract the mean climate values from within a 25km buffer of the stopover centroid, and if no climate data within 25km of the stopover, then increase the buffer until climate data is available and extract it
    ### CORINE: Random polygon placement means that some random stopovers may be placed in areas which either have poor terrestrial coverage (> 50% water) or poor coverage by the Corine raster. For generaterandom=TRUE, want to check that random stopovers pick up corine (as well as climate) raster cells
    
    ###--- FUNCTIONS TO EXTRACT CLIMATE DATA ---###
    
    # FUNCTION extracts the raster climate data for each stopover polygon for an individual, takes the mean of all the raster cells within the buffer zone    
    extract.climate <- function(rastername, polygonname){
      extracted.data <- extract(rastername, coordinates(polygonname), small=TRUE, buffer=buffer.metres[countinner], fun=mean)
      #extracted.data <- extract(rastername, coordinates(polygonname), small=TRUE, buffer=25000, fun=mean) 
      colnames(extracted.data) <- c("year1","year2","year3")
      rownames(extracted.data) <- 1:nrow(extracted.data)  
      mgroupid <- polygonname@data  
      extracted.data.withid <- data.frame(extracted.data, mgroup=mgroupid)
      return(extracted.data.withid)
    }
    
    # FUNCTION converts the list of mean climate data per stopover from the previous function to a data.frame showing mgroup, year, season, precip, and temp data
    convertclimate <- function(output){
      nameslist <- strsplit(names(output), ".", fixed=TRUE)
      namesmatrix <- do.call(rbind, nameslist)
      
      all <- do.call(rbind, output)
      season <- rep(namesmatrix[,2], sapply(output, nrow))
      climvar <- rep(namesmatrix[,3], sapply(output, nrow))
      
      climate.dataset <- data.frame(all, season, climvar)      
      climate.dataset2 <- melt(climate.dataset, id=c("mgroup","season","climvar"))      
      climate.dataset3 <- reshape(climate.dataset2, v.names="value", timevar="climvar", idvar=c("season","mgroup","variable"), direction="wide")
      names(climate.dataset3) <- c("mgroup","season","year","precip.mm","temp.C")      
      climate.dataset3$year <- factor(climate.dataset3$year, labels=c("2011","2012","2013"))
      
      return(climate.dataset3)
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
    
    ###--------- CLIMATE AND CORINE EXTRACTION -----------
    
    
    #--- IF GENERATERANDOM = TRUE ---#
    # run source file to generate random buffered polygons and extract climate data such that all stopovers have associated climate data at the smallest buffer level that generates non-NA values
    if (generaterandom) { # START generaterandom = TRUE loop
      
      repeat { # OUTER repeat loop, for corine extractions > 50% water
        
        #--- CLIMATE (extracted inside source) ---#
        setwd(paste(parentwd, "/scripts", sep=""))
        source("cuckoo_generate_random_buffered_polygons.R")
        #     randomtry[a] <- countouter
        
        #######################################################
        
        #--- CORINE ---#
        
        #       setwd(GISwd)
        #       writeOGR(randomstopSPDFs, dsn=".", layer = "cuckoo Skinner random stopovers 200km buffer 3", driver = "ESRI Shapefile")
        #       
        #       corine.values[[a]] <- extract(r, randomstopSPDFs)
        
        corine.values[[a]] <- newextract(r, randomstopSPDFs)
        
        check <- logical()
        
        for (i in 1:length(corine.values[[a]])) {
          if ( (length(which(corine.values[[a]][[i]]==44 | corine.values[[a]][[i]]==39 | corine.values[[a]][[i]]==40 | corine.values[[a]][[i]]==41 | corine.values[[a]][[i]]==42 | corine.values[[a]][[i]]==43 | is.na(corine.values[[a]][[i]])))/length(corine.values[[a]][[i]]) < 0.5 ) ) {
            check[i] <- TRUE
          } else {
            check[i] <- FALSE
          }
        }
        
        # checks whether all random stopovers have < 50% water
        if (all(check==TRUE)) break # breaks OUTER repeat loop, for corine extractions
        
      } # END outer repeat loop, for corine extractions
      
      corine.MCP <- list()
      
      for (i in 1:length(corine.values[[a]])) {
        corine.MCP[[i]] <- data.frame(
          name=levels(newdat@data$name), 
          mgroup=levels(newdat@data$mgroup)[i], 
          year=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "year"][1], 
          laststop=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "laststop"][1],
          strategy=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "strategy"][1],
          Sahara.success=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "Sahara.success"][1],
          corine.values=corine.values[[a]][[i]])
      }
      
    } # END generaterandom = TRUE loop
    
    #######################################################
    #######################################################
    
    #--- IF GENERATERANDOM = FALSE ---#
    if (!generaterandom) { # START generaterandom = FALSE loop
      
      #--- CLIMATE ---#
      
      allstop4326 <- spTransform(allstopSPDFs, CRS=epsg4326)
      
      buffer.metres <- c(25000,50000,75000,100000,125000,150000) # different than the random stopovers code, which only buffers up to 50km.  Less leeway with the actual stopovers, so buffer up to 150km in the hopes of finding a non-NAN raster cell, if necessary
      
      countinner <- 1 # counts the buffer levels, starts at smallest buffer (25 km)
      
      extracted <- lapply(climate.data, extract.climate, allstop4326)
      extract.clim.values <- convertclimate(extracted)
      
      # if all mean extracted climate values for all stopovers an individual makes have values (are not NA), then the climate data extraction process is complete
      # if some stopovers are NA, then enter repeat loop, first identifying the problem stopovers (the ones with NAs)
      
      if (all(!is.na(extract.clim.values$precip.mm))) {
        climate.complete <- extract.clim.values
      } else {
        
        # writes non-NA mgroups to compiled dataset
        climate.complete <- subset(extract.clim.values, !is.na(extract.clim.values$precip))
        
        # identify mgroups with NAs
        problemmgroups <- subset(extract.clim.values, is.na(extract.clim.values$precip.mm), mgroup, drop=TRUE)
        problemmgroups <- droplevels(problemmgroups)
        probmgroupids <- levels(problemmgroups)
        
        repeat { # REPEAT - cycles through buffer levels for problem mgroups
          
          # for mgroups with NAs, subset the allstopover SPDF to the problem mgroups, increase the buffer level by 1, and re-extract the climate raster
          probstops <- subset(allstop4326, mgroup %in% probmgroupids)
          
          countinner <- countinner + 1
          extracted <- lapply(climate.data, extract.climate, probstops)
          extract.clim.values <- convertclimate(extracted)
          
          # if all these problem mgroups now have values (ie. all are non-NA), then add the climate data for these stopovers to climate.complete and break out of the repeat loop
          if (all(!is.na(extract.clim.values$precip.mm))) {
            climate.complete <- rbind(climate.complete, extract.clim.values)
            break # breaks repeat loop when all mgroups have data
            
          } else {
            
            # if there is still > 1 mgroup with a problem, then need to continue loop; first, check if any mgroups that were a problem are now ok (ie. are non-NA) at the new buffer level. For ones that are ok, add its data to climate.complete, and repeat the loop with the mgroups that are still a problem
            okmgroups <- subset(extract.clim.values, !is.na(extract.clim.values$precip.mm))
            climate.complete <- rbind(climate.complete, okmgroups)
            
            # identify mgroups still with NAs
            problemmgroups <- subset(extract.clim.values, is.na(extract.clim.values$precip.mm), mgroup, drop=TRUE)
            problemmgroups <- droplevels(problemmgroups)
            probmgroupids <- levels(problemmgroups)
            
          } # end inner IF statement (if problem mgroups are all still NA)
        } # end REPEAT loop - cycles through buffer levels
      } # end outer IF statement (if some mgroups are NA)
      
      ### CLIMATE SUMMARY, ACTUAL STOPOVERS ###
      # 18 birds captured non-NAN cells at 25km, 5 at 50km, 4 at 75km, 1 at 100km, 1 at 150km
      # See per bird information summary file at C:\Users\samf\Documents\Git\cuckoos\output\all cuckoo stopovers and non-NAN climate raster buffer values
      
      #######################################################
      
      #--- CORINE ---#
      
      corine.values[[a]] <- extract(r, allstopSPDFs)
      
      corine.MCP <- list()
      
      for (i in 1:length(corine.values[[a]])) {
        corine.MCP[[i]] <- data.frame(
          name=levels(newdat@data$name), 
          mgroup=levels(newdat@data$mgroup)[i], 
          year=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "year"][1], 
          laststop=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "laststop"][1],
          strategy=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "strategy"][1],
          Sahara.success=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "Sahara.success"][1],
          corine.values=corine.values[[a]][[i]])
      }
      
      
    } # end IF !generaterandom
    
    #######################################################
    #######################################################
    
    ###--- COMPILE FULL CLIMATE AND CORINE DATASETS FOR INDIVIDUAL ---###
    
    ### CLIMATE: Create full dataset to output with extra variables, only output climate data for relevant years (years with migration info)  
    climate.values <- merge(data.frame(name=levels(newdat@data$name), climate.complete), mgroupyear.combn[[a]])
    
    fullclimate[[a]] <- climate.values[order(rev(climate.values$season), as.numeric(as.character(climate.values$mgroup))),]
    
    ### CORINE
    fullcorine[[a]] <- do.call(rbind,corine.MCP)
    
    #######################################################
    #######################################################
    #######################################################
    
    ###-------- REAL OR RANDOM STOPS --------
    
    if (generaterandom) {
      mcp <- randomstopSPDFs
    }
    
    if (!generaterandom) {
      mcp <- allstopSPDFs
    }
    
    ###-------------------------------###
    
    #######################################################
    #######################################################
    #######################################################
    
    ###---------------------------###
    #----   ELEVATION  ----
    ###---------------------------###
    
    mcp.elevtransform <- spTransform(mcp, CRS=elevation@crs)
    
    elevation.values[[a]] <- extract(elevation, mcp.elevtransform)
    
    elevation.dataset <- list()
    
    for (i in 1:length(elevation.values[[a]])) {
      elevation.dataset[[i]] <- data.frame(
        name=levels(newdat@data$name), 
        mgroup=levels(newdat@data$mgroup)[i], 
        year=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "year"][1], 
        laststop=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "laststop"][1],
        strategy=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "strategy"][1],
        Sahara.success=newdat@data[which(newdat@data$mgroup==levels(newdat@data$mgroup)[i]), "Sahara.success"][1],
        elevation.m=elevation.values[[a]][[i]])
    }
    
    fullelevation[[a]] <- do.call(rbind, elevation.dataset)
    
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
      MCPoverlap[[a]] <- data.frame(
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
      
      
      MCPoverlap[[a]] <- data.frame(
        name=levels(newdat@data$name), 
        MCPid=mcp.id, 
        MCP.area.km2, 
        PA.overlap.area.km2=PA.overlap.area2, 
        prop.MCP.overlap
      )
      
      # attribute info about overlapping PAs
      indPA.overlap <- gArea(PAintersection2, byid=TRUE)/1000000 
      #   names(indPA.overlap) <- names(PAintersection)
      
      
      indPA.overlap2 <- indPA.overlap
      
      PAoverlap.info[[a]] <- data.frame(birdname=levels(newdat@data$name),mgroup=parentMCP,PA.subsetintersect@data[,c("wdpaid","country","name","desig_eng","desig_type","iucn_cat","rep_area")], MCP.overlap.area.km2=indPA.overlap2)
      
    }
    
  }  # end FOR loop counting through individuals 1:31
  
  #######################################################
  #######################################################
  #######################################################
  
  #################################
  #
  #### WRITE DATASETS TO CSV  ####
  #
  #################################
  
  ###--- CLIMATE ---###
  
  setwd(paste(datawd,"/resampled data corine PA elevation climate extracted values/climate",sep=""))
  all.fullclimate <- do.call(rbind, fullclimate)
  
  if (!generaterandom) {
    write.csv(all.fullclimate, "climate all birds.csv", row.names=FALSE)
  }
  
  if (generaterandom) {
    write.csv(all.fullclimate, paste("climate random ", randomradius[z]/1000, " km radius all birds.csv", sep=""), row.names=FALSE)
  }
  
  ###--- CORINE ---###
  
  # for corine, each bird needs its own folder (can't merge all birds into one dataset, too many observations)
  
  if (!generaterandom) {
    for (i in 1:31){
      if (length(fullcorine[[i]]) == 0) {
        next
      } else {
        setwd(paste(datawd,"/resampled data corine PA elevation climate extracted values/corine",sep=""))
        write.csv(fullcorine[[i]], paste(fullcorine[[i]]$name[1],".csv", sep=""), row.names=FALSE)
      }
    }
  }
  
  if (generaterandom) {
    for (i in 1:31){
      if (length(fullcorine[[i]]) == 0) {
        next
      } else {
        
        setwd(paste(datawd,"/resampled data corine PA elevation climate extracted values/corine random", sep=""))
        
        if(!(paste(randomradius[z]/1000, " km radius", sep="") %in% list.files())){
          dir.create(paste(randomradius[z]/1000, " km radius", sep=""))
          setwd(paste(randomradius[z]/1000, " km radius", sep=""))
          
        } else {
          setwd(paste(randomradius[z]/1000, " km radius", sep=""))
        }
        
        write.csv(fullcorine[[i]], paste(fullcorine[[i]]$name[1],".csv", sep=""), row.names=FALSE)
      }
    }
  }
  
  ###--- ELEVATION ---###
  
  setwd(paste(datawd,"/resampled data corine PA elevation climate extracted values/elevation",sep=""))
  all.fullelevation <- do.call(rbind, fullelevation)
  
  if (!generaterandom) {
    write.csv(all.fullelevation, "elevation all birds.csv", row.names=FALSE)
  }
  
  if (generaterandom) {
    write.csv(all.fullelevation, paste("elevation random ", randomradius[z]/1000, " km radius all birds.csv", sep=""), row.names=FALSE)
  }
  
  
  ###--- PROTECTED AREAS ---###
  
  #---- Overlapping protected area site details ----#
  
  ### for all birds, combine details about PAs which overlap with stopover MCPs
  
  setwd(paste(datawd,"/resampled data corine PA elevation climate extracted values/protected areas",sep=""))
  allbirds.PAoverlap.info <- do.call(rbind, PAoverlap.info)
  
  if (!generaterandom){
    write.csv(allbirds.PAoverlap.info, "PA site details all birds.csv", row.names=FALSE)
  }
  
  if (generaterandom){
    write.csv(allbirds.PAoverlap.info, paste("PA site details random ", randomradius[z]/1000, " km radius all birds.csv", sep=""), row.names=FALSE)
  }
  
  #---- Stopover overlap details ----#
  
  ### combine all MCPoverlap datasets into one, and merge with other useful stopover information (year, laststop, Sahara success, etc)
  
  allbirds.MCPoverlap <- do.call(rbind, MCPoverlap)
  
  stopoverinfo <- do.call(rbind, mgroupyear.combn)
  
  allstopoverinfo <- merge(x=stopoverinfo, y=allbirds.MCPoverlap, by.x=c("name","mgroup"), by.y=c("name","MCPid"))
  
  allstopoverinfo$mgroup <- as.numeric(as.character(allstopoverinfo$mgroup))
  
  allstopoverinfo.ordered <- allstopoverinfo[order(allstopoverinfo$name,allstopoverinfo$mgroup),]
  
  if (!generaterandom){
    write.csv(allstopoverinfo.ordered, "stopover overlap area all birds.csv", row.names=FALSE)
  }
  
  if (generaterandom){  
    write.csv(allstopoverinfo.ordered, paste("stopover overlap area random ", randomradius[z]/1000, " km radius all birds.csv", sep=""), row.names=FALSE)
  }
  
}