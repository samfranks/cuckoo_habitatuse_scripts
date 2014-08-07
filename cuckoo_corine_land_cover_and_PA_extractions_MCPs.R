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

cluster <- FALSE

#library(adehabitat)
library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
#library(shapefiles) # allows you to read and write ESRI shapefiles
library(rgeos)
library(raster)
library(geosphere)
library(adehabitatHR)

# allPAobjects <- which(ls() == "PA.subset" | ls() == "r" | ls() == "Europeraster" | ls() == "corine.crs")
# 
# allPAobjects <- which(ls() == "PA.subset" | ls() == "r" | ls() == "Europeraster" | ls() == "corine.crs")
# 
# toremove <- ls()[-allPAobjects]
# 
# rm(list=toremove)

###-----------------------------------------------------------###
#         SET WORKING DIRECTORIES
###-----------------------------------------------------------###

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

corinewd <- paste(parentwd, "/data/corine raster 100x100", sep="")

if (!cluster) {
  GISwd <- c("C:/GIS/cuckoos")
}

if (cluster) {
  GISwd <- c("/users1/samf/cuckoos")
}

corine.crs <- CRS("+init=epsg:3035")

setwd(parentwd)

###--------------------------------------------------------###
#        LOAD BASE LAYERS (if needed)
###--------------------------------------------------------###

if (!("r" %in% ls())) { ### run prep code to load base layers if not in workspace already
  
  ### CORINE layer ###
  
  # projection and datum information for all maps
  # epsg: 3395 - +proj=merc +ellps=wGS84 +datum=WGS84
  # epsg: 4326 - +proj=latlong +ellps=wGS84 +datum=WGS84
  # epsg: 3035 - LAEA CRS, used for the Corine Land Cover 2006 data
  
  ### ------- CORINE land cover raster file ------- ###
  
  setwd(corinewd)
  r <- raster("g100_06.tif")
  
  Europebox <- extent(2483500,5890375,1276625,4286750) # bounding box extent for relevant parts of Europe (no Scandanavia or Turkey)
  #Europeraster <- crop(r,Europebox)
  
}

### ------- Europe/Africa national boundaries layer ------- ###

# if (!("Eur.Afr" %in% ls())) { ### run prep code to load base layers if not in workspace already
#   
#   dat.world <- getMap(resolution="high")
#   
#   Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
#   Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
#   Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]
#   
#   Eur.Afr <- spTransform(Eur.Afr,CRS("+proj=longlat +datum=WGS84")) # transform map to exactly the same projection as cuckoo newlongs/newlats data
#   
#   Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
#   Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]
#   
# }
  
### ------- Protected area GIS layer (raster or shapefile) ------- ###

if (!("PA.subset" %in% ls())) { ### run prep code to load base layers if not in workspace already
  
  #   GISwd <- c("C:/GIS/cuckoos")
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

##############################################################
##############################################################
#
# EXTRACT CORINE & PA VALUES FOR MCPs (resampled data only)
#
##############################################################
##############################################################

### TO RUN ORIGINAL or RESAMPLED data, change below line

original <- FALSE # or FALSE if resampled

### TO RUN RANDOM MCPs, CHANGE BELOW LINE

generaterandom <- FALSE
randomradius <- 50000 # change to desired search radius for randomized pseudoabsence polygon, in metres

### TO EXTRACT DATA FROM 500m BUFFERED POINT SPATIAL POLYGONS
buffered <- TRUE

# create list to hold datasets with extracted corine and protected area values for each cuckoo

corine.values <- list()
fullcorine.MCP <- list()
MCPoverlap <- list()
PAoverlap.info <- list()
fullcorine.randomMCP <- list()

# for loop for each cuckoo file starts here
#for (a in 1:31) { # START LOOP

a<-2
  
  
  ###--- LOAD DATA, ADD PROJECTION INFO, TRANSFORM ---###
  
  if (original) setwd(paste(datawd,origwd,sep=""))
  if (!original) setwd(paste(datawd,resampwd,sep=""))
  
  dataset <- read.csv(list.files()[a], header=T)
  
  # long/lat (original) are in CRS("+init=epsg:4326")
  # long/lat (resampled) are in CRS("+init=epsg:3395")  
  # newlongs/newlats (resampled) are in CRS("+init=epsg:4326")
  # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
  
  ### check that dataset is not an individual to exclude (Idemili & Karma); break out of current loop run and continue to next loop level
  if (dataset$name[1] == "Karma" | dataset$name[1] =="Idemili") {
    next
  }
  
  # ---------------------- RESAMPLED DATA -----------------------#
  if (!original) {
    coordinates(dataset) <- c("newlongs","newlats")
    proj4string(dataset) <- CRS("+init=epsg:4326") 
    newdataset <- spTransform(dataset, CRS=corine.crs) 
  }
  
  # ---------------------- ORIGINAL DATA -----------------------#
  if (original) {
    coordinates(dataset) <- c("long","lat")
    proj4string(dataset) <- CRS("+init=epsg:4326")
    newdataset <- spTransform(dataset, CRS=corine.crs)
    
  }
  
  # subset dataset by stopover sites (no breeding sites)
  newdat <- subset(newdataset, stopoversite=="Y")
  newdat@data <- droplevels(newdat@data)

###-----------------------------------------------###
#   CALCULATE BUFFERED POINTS USED FOR RESAMPLED STOPOVER LOCATIONS
###-----------------------------------------------###


if (nrow(newdat) == 0) {
  next
}

######## NOT USING THIS KERNEL METHOD ############
# # for each movement group of each cuckoo, calculate the UD using function kernelUD() in library(adehabitatHR)
# newdat@data$mgroup <- as.factor(newdat@data$mgroup)
# stoppoints <- newdat[,c("mgroup")]
# kud <- kernelUD(stoppoints, h="href")
# # kud1 <- kernelUD(stoppoints, h="LSCV") # this doesn't work - doesn't minimize
# homerange <- getverticeshr(kud)
##################################################

mgroupSPDFs <- list()

for (i in 1:length(levels(newdat$mgroup))) {
  newdatsub <- subset(newdat, mgroup==levels(newdat$mgroup)[i])
  newdatsub@data <- droplevels(newdatsub@data)
  mgroupSPDFs[[i]] <- gBuffer(newdatsub, width=500, id=levels(newdat$mgroup)[i])
}

allstopSPDFs <- do.call(spRbind, mgroupSPDFs)

corine.values <- extract(r, allstopSPDFs)

#sapply(slot(allstopSPDFs, 'polygons'), function(i) slot(i, 'ID'))

# plot(croppedraster)
# plot(allstopSPDFs, add=T, lwd=1.5)
# plot(newdat, add=T, pch=16, cex=0.6)

  
  ###-----------------------------------------------###
  #   CALCULATE MCP FOR RESAMPLED STOPOVER LOCATIONS
  ###-----------------------------------------------###
  
  ### skip to next bird if no stopover sites ###
  if (nrow(newdat) == 0) {
    next
  }
  
  # for each movement group of each cuckoo, draw a MCP around the mgroup's points using mcp() to estimate home range
  # choose points for MCP such that 5% of outliers are excluded
  newdat@data$mgroup <- as.factor(newdat@data$mgroup)
  stoppoints <- newdat[,c("mgroup")]
  MCP <- mcp(stoppoints, percent=100, unin=c("m"), unout=c("km2"))
  
  ###--- LAND-COVER (CORINE) EXTRACTIONS ---###
  
  #  INSERT SOURCE FILE HERE TO GENERATE RANDOM MCPs AND EXTRACT CORINE DATA
  
  if (generaterandom) {
    setwd(paste(parentwd, "/scripts", sep=""))
    source("cuckoo_generate_random_MCPs.R")
  }
  
  if (!generaterandom) {  
    corine.values[[a]] <- extract(r, MCP)
  }
  
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
  
  if (generaterandom) {
    fullcorine.randomMCP[[a]] <- do.call(rbind,corine.MCP)
  }
  
  if (!generaterandom) {
    fullcorine.MCP[[a]] <- do.call(rbind,corine.MCP)
  }
  
  
  ###--- PROTECTED AREA MCP EXTRACTIONS ---###
  
  if (generaterandom) {
    mcp <- randomMCP
  }
  
  if (!generaterandom) {
    mcp <- MCP
  }
    
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
      MCPid=mcp$id, 
      MCP.area.km2=mcp$area, 
      PA.overlap.area.km2=rep(0, nrow(mcp@data)), 
      prop.MCP.overlap=rep(0, nrow(mcp@data))
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
    
    PA.overlap.area2[which(mcp$id %in% names(PA.overlap.area))] <- PA.overlap.area
    PA.overlap.area2[which(!(mcp$id %in% names(PA.overlap.area)))] <- 0
    
    # max(gArea(PAintersection, byid=TRUE))/1000000
    
    # proportion area of MCP overlapping with PA (for MCPs which intersect with a PA)
    #MCP.area.km2 <- mcp@data[which(names(PA.overlap.area) %in% mcp$id),c("area")]
    MCP.area.km2 <- mcp$area
    prop.MCP.overlap <- PA.overlap.area2/MCP.area.km2
    
    
    MCPoverlap[[a]] <- data.frame(
      name=levels(newdat@data$name), 
      MCPid=mcp$id, 
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
    
    ### --- CROP Europeraster to box around stopover area --- ###
    
    # clip the raster to a bounding box with an extent that buffers each MCP by ~ 5km and plot the clipped rasters with the stopover points on it, and export to a TIFF
    
    # for each polygon in MCP, create a vector of the min(newlongs),max(newlongs),min(newlats),max(newlats), and use this to crop the raster Europeraster
    # clip the raster to a bounding box with an extent that buffers the MCP by ~ 10km
    
    MCPbox <- list()
    
    # ??? can convert below loop to a function if desired
    
    bufferradius <- 10000
    
    for (i in 1:length(mcp)){
      polycoords <- mcp@polygons[[i]]@Polygons[[1]]@coords
      MCPbox[[i]] <- c(min(polycoords[,1])-bufferradius, max(polycoords[,1])+bufferradius, min(polycoords[,2])-bufferradius, max(polycoords[,2])+bufferradius)
      MCPbox
    }
    
    clipMCP <- function(MCPminmax){
      MCPbounds <- extent(MCPminmax)
      return(MCPbounds)
    }
    
    mcp.bbox <- lapply(MCPbox,clipMCP)
    
      
      # test plot (change mcp$id to appropriate mgroup number)
    
    croppedraster <- crop(r,mcp.bbox[[2]])
    
    xlims <- c(mcp.bbox[[2]]@xmin, mcp.bbox[[2]]@xmax)
    ylims <- c(mcp.bbox[[2]]@ymin, mcp.bbox[[2]]@ymax)
  
      pa <- PA.subset[which(subsetintersect[,c("26")]==TRUE),]
    # plot(croppedraster)
      plot(mcp[mcp$id==26,], border="red", lwd=2, xlim=xlims, ylim=ylims)
      plot(pa, add=T, border="darkgreen", lwd=2, xlim=xlims, ylim=ylims)
      plot(newdat[newdat@data$mgroup==26,], add=T, pch=16, cex=0.4, xlim=xlims, ylim=ylims)
    plot(MCP[MCP$id==26,], border="black", add=T, lwd=1.5)
      plot(PAintersection, add=T, border="blue", col= "blue", xlim=xlims, ylim=ylims)
    
    setwd(mapoutputdir)
    write.csv(PAoverlap.info[[a]], "Chance mgroup 29 example.csv", row.names=FALSE)
      
    


} # END LOOP through each cuckoo file



###----------------------------------------------------###
#       WRITE DATASETS TO CSV
###----------------------------------------------------###

###--- corine values ---###

if (!generaterandom) {
  for (i in 1:31){
    if (length(fullcorine.MCP[[i]]) == 0) {
      next
    } else {
      setwd(paste(datawd,"/resampled data corine and PA extracted values/corine",sep=""))
      write.csv(fullcorine.MCP[[i]], paste(fullcorine.MCP[[i]]$name[1],".csv", sep=""), row.names=FALSE)
    }
  }
}

if (generaterandom) {
  for (i in 1:31){
    if (length(fullcorine.randomMCP[[i]]) == 0) {
      next
    } else {
      
      setwd(paste(datawd,"/resampled data corine and PA extracted values/corine random", sep=""))
      
      if(!(paste(randomradius/1000, " km radius", sep="") %in% list.files())){
        dir.create(paste(randomradius/1000, " km radius", sep=""))
        setwd(paste(randomradius/1000, " km radius", sep=""))
        
      } else {
        setwd(paste(randomradius/1000, " km radius", sep=""))
      }
      
      write.csv(fullcorine.randomMCP[[i]], paste(fullcorine.randomMCP[[i]]$name[1],".csv", sep=""), row.names=FALSE)
    }
  }
}


###--- overlapping protected areas ---###

### for all birds, combine details about PAs which overlap with stopover MCPs

setwd(paste(datawd,"/resampled data corine and PA extracted values",sep=""))

if (!generaterandom){
  allbirds.PAoverlap.info <- do.call(rbind, PAoverlap.info)
  write.csv(allbirds.PAoverlap.info, "allbirds PA overlap details.csv", row.names=FALSE)
}

if (generaterandom){
  allbirds.PAoverlap.info <- do.call(rbind, PAoverlap.info)
  write.csv(allbirds.PAoverlap.info, paste("random_", randomradius/1000, " km_", "allbirds PA overlap details.csv", sep=""), row.names=FALSE)
}

### combine all MCPoverlap datasets into one, and merge with other useful stopover information (year, laststop, Sahara success, etc)

allbirds.MCPoverlap <- do.call(rbind, MCPoverlap)
stopoverinfo <- read.csv(paste(datawd,"/cuckoo stopover data - mgroup_year_laststop_etc.csv", sep=""), header=T)

stopoverinfo$mgroup <- as.factor(stopoverinfo$mgroup)

allstopoverinfo <- merge(x=stopoverinfo, y=allbirds.MCPoverlap, by.x=c("name","mgroup"), by.y=c("name","MCPid"))

allstopoverinfo$mgroup <- as.numeric(as.character(allstopoverinfo$mgroup))

allstopoverinfo.ordered <- allstopoverinfo[order(allstopoverinfo$name,allstopoverinfo$mgroup),]

setwd(paste(datawd,"/resampled data corine and PA extracted values",sep=""))

if (!generaterandom){
  write.csv(allstopoverinfo.ordered, "allbirds MCP overlap data.csv", row.names=FALSE)
}

if (generaterandom){  
  write.csv(allstopoverinfo.ordered, paste("random_", randomradius/1000, " km_", "allbirds MCP overlap data.csv", sep=""), row.names=FALSE)
}

