##########################################################
#
#  Analysis examining:
#       *1) proportional habitat use of original vs resampled data
#       2) proportional habitat use vs habitat availability in actual stopover polygon
#       3a) proportional habitat availability in actual stopover polygon vs random stopover polygons (actual vs null) for birds using different migration strategies
#       3b) proportional habitat availability in actual stopover polygon vs random stopover polygons (actual vs null) for birds with failed vs successful Sahara crossing
#       4) proportional habitat use of birds using different migration strategies: a) all stopovers; b) last stopover
#       5) proportional habitat use of birds with failed vs successful Sahara crossing: a) all stopovers; b) last stopover
#       6) proportional habitat availability 
#
# 7) compare habitat availability at 1) last stopovers for a) birds with different migration strategies; b) successful/failed birds and 2) all other stopovers for a) different mig strategies; b) successful/failed birds
#
#  Samantha Franks
# 12 Dec 2013
#
##########################################################

######### NOTES ABOUT THE DATASETS

# README
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
# long/lat (resampled) = epsg: 3395 +proj=merc +ellps=wGS84 +datum=WGS84
# long/lat (original) = epsg: 3035 LAEA CRS, used for the Corine Land Cover 2006 data
# continent, country = continent and country an observation falls in, extracted previously
# newlongs/newlats = resampled location points, in projection epsg: 3035 LAEA CRS, used for the Corine Land Cover 2006 data
# corine.values = corine grid code cell values extracted from the raster underlying point locations, lookup meaning of value in clclegend.csv
# age = bird's age at capture
# breedsite = Y for mgroups falling at UK breeding site, otherwise N
# strategy = migration strategy used (SE or SW), Balkans is lumped in with SE
# Sahara.success = Y is a successful crossing in a given year, N is an unsuccessful crossing in a given year
# laststop = whether an mgroup was the last stopover in Europe in a given year (Y or N)
# corine.name = the clclegend description of the grid code value, at its most detailed levels (LABEL1 and LABEL2 or less detailed descriptions)

library(adehabitat) # 
library(sp) # converts R dataframes into spatially explicit dta
library(rgdal)
#library(shapefiles) # allows you to read and write ESRI shapefiles
library(rgeos)
library(raster)
library(adehabitatHR)
library(NCStats)
library(shape)
library(splancs)
library(maptools)

### --- Set-up data: bring in resampled and original data with corine data + other variables --- ###


### --- Set working directories --- ###

setwd("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos")
parentwd <- getwd()

datawd <- ("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/data")
outputwd <- ("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/output")

resampwd <- ("/resampled data corine + extra grouping variables")
origwd <- ("/original data corine + extra grouping variables")

corinewd <- c("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/data/corine raster 100x100")


### ------- CORINE land cover raster file ------- ###
setwd(corinewd)
r <- raster("g100_06.tif")

Europebox <- extent(2483500,5890375,1276625,4286750) # bounding box extent for relevant parts of Europe (no Scandanavia or Turkey)
Europeraster <- crop(r,Europebox)

### import corine legend
clcwd <- c("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/data/")

setwd(clcwd)
clclegend <- read.csv("clc_legend.csv", header=T)


##################################################################
#
# Analysis 1: Compare proportional habitat use of original vs resampled data
#
# Chi-sq test of original habitat class frequencies vs resampled habitat class frequencies at European stopovers (breedsite=N), for each individual
#
##################################################################

# could not directly compare observed vs expected distribution of proportional habitat use
# zeros in the expected frequency distribution made test-statistics and p-values meaningless

# instead, I tested whether the absolute difference in the proportion of a habitat used in the original vs bootstrapped data was the same across all habitats, ie. the goodness of fit test tests whether the absolute magnitude of the difference between the proportion of original vs resampled habitat used is the same across all habitat types, or if use of one (or more) habitats are disproportionately over- or under-represented by the bootstrapping procedure (if the difference departs from an expected equal difference between observed/expected values across all habitat types)

# goodness of fit test produces a warning (chi-squared approximation may be incorrect) because of the small observed and expected frequencies associated with the magnitude of the difference in proportions

for (a in 1:29) {
  
  ### select cuckoo dataset (leave out Idemili and Karma)
  
  setwd(paste(datawd,resampwd,sep=""))
  resampledcorine <- list.files()[-8] # remove Idemili
  resampledcorine <- resampledcorine[-11] # remove Karma
  rdat <- read.csv(resampledcorine[a], header=T) # resampled data
  
  setwd(paste(datawd,origwd,sep=""))
  originalcorine <- list.files()[-8]
  originalcorine <- originalcorine[-11]
  odat <- read.csv(originalcorine[a], header=T) # original data
  
  ### --- Goodness of fit test original habitat class frequences vs resampled --- ###
  
  ### all European stopovers pooled
  
  rstopovers <- Subset(rdat, breedsite=="N")
  ostopovers <- Subset(odat, breedsite=="N")
  
  rprop <- as.data.frame(prop.table(table(rstopovers$corine.name)))
  colnames(rprop) <- c("corine.name","freq")
  oprop <- as.data.frame(prop.table(table(ostopovers$corine.name)))
  colnames(oprop) <- c("corine.name","freq")
  
  omissing <- rprop[which(!rprop$corine.name %in% oprop$corine.name),"corine.name"]
  rmissing <- oprop[which(!oprop$corine.name %in% rprop$corine.name),"corine.name"]
  
  oprop2 <- rbind(oprop,data.frame(corine.name=omissing,freq=rep(0,length(omissing))))
  rprop2 <- rbind(rprop, data.frame(corine.name=rmissing, freq=rep(0,length(rmissing))))
  
  # reorder both so that corine.names are in same alphabetical order
  orderedcorinenames <- sort(as.character(oprop2$corine.name))
  
  oprop2$corine.name <- factor(oprop2$corine.name, orderedcorinenames, ordered=TRUE)
  oprop3 <- oprop2[order(oprop2$corine.name),]
  
  rprop2$corine.name <- factor(rprop2$corine.name, orderedcorinenames, ordered=TRUE)
  rprop3 <- rprop2[order(rprop2$corine.name),]
  
  diff <- data.frame(corine.name=rprop3$corine.name, difference=abs(rprop3$freq - oprop3$freq))
  
  chisqout <- chisq.test(diff$difference)
  
  setwd(paste(outputwd,"/Europe habitat use/",sep=""))
  
  cat(paste(levels(rdat$name), " : habitat use of resampled vs original location points", "\n", "Difference in use of land cover classes", "\n", sep=""), file="chisqtest.txt", sep="\n\n", append=TRUE)
  out<-capture.output(chisqout)
  cat(out, file="chisqtest.txt",sep="\n",append=TRUE)
  
  out<-capture.output(data.frame(corine.name=rprop3$corine.name, resampled=rprop3$freq, original=oprop3$freq, difference=(rprop3$freq-oprop3$freq)))
  cat(out, file="chisqtest.txt",sep="\n", append=TRUE)
  
  cat(c("########################################\n########################################\n########################################"), file="chisqtest.txt", sep="\n", append=TRUE)
  
}

############------------- SUMMARY OF ANALYSIS 1 -------------############

# no significant difference found in the proportional use of different habitats between the original location points and the bootstrapped points
# zeros were substituted into the original habitat use table for habitats picked up in the bootstrapping procedure that weren't sampled by the original location points

##################################################################
#
# For subsequent analyses using habitat availability: extract corine values for original stopover polygons and create random stopover polygons and extract corine values
#
#
##################################################################

for (a in 1:29) {
  
  setwd(paste(datawd,resampwd,sep=""))
  resampledcorine <- list.files()[-8] # remove Idemili
  resampledcorine <- resampledcorine[-11] # remove Karma
  rdat <- read.csv(resampledcorine[a], header=T) # resampled data
  corine.crs <- CRS("+init=epsg:3035")
  
  ###----------------------------------------------------------------###
  # STEP 1: subset to include only stopovers (breedsite=N), convert dataframe to spdf
  ###----------------------------------------------------------------###
  
  stopovers <- Subset(rdat, breedsite=="N")
  coordinates(stopovers) <- c("newlongs","newlats")
  proj4string(stopovers) <- corine.crs
  
  #   stopovers$mgroup <- as.factor(stopovers$mgroup)
#   bymgroup <- stopovers$mgroup
#   splitstop <- split(stopovers, list(bymgroup))
 
  
  ###----------------------------------------------------------------###
  # STEP 2: calculate MCP around each stopover, and extract centroid long/lat values to use in randomization procedure
  ###----------------------------------------------------------------###
    
  # for each movement group of each cuckoo, draw a MCP around the mgroup's points using mcp() to estimate home range
  # choose points for MCP such that 5% of outliers are excluded
  stopid <- stopovers[,c("mgroup")]
  MCP <- mcp(stopid, percent=100)
  
  ###----------------------------------------------------------------###
  # STEP 3: draw a circle with centroid coordinates that are the centroid of the MCP using library(shape), then use csr() from library(splancs) to generate a random point somewhere within that circle polygon that will represent the new centroid coordinates of the randomized polygon
  # calculate the difference (in metres) between the new randomized centroid coordinates and the original stopover polygon's centroid coordinates
  ###----------------------------------------------------------------###
  
  # for each Polygons object in MCP@polygons, create function that will randomly shift the original stopover polygon somewhere within a 50km radius of the original MCP midpoint
  
  getnewMCP <- function(spdf) {
  circleMCP <- getellipse(rx=50000, ry=50000, mid=spdf@labpt)
  
  randpts <- csr(circleMCP,1)
  centdiff <- spdf@labpt-randpts # difference in metres between original MCP centroid and new MCP centroid
  
  origMCP <- SpatialPolygons(list(spdf))
  proj4string(origMCP) <- corine.crs
  
  # use elide() and argument shift=c() in package(maptools) to shift centroid coordinates of polygon by c(x,y) amount (centdiff) - c(x,y) is the difference between the centroid of original MCP and random pt in 50km radius circle around centroid of original MCP
  newMCP <- elide(origMCP, shift=c(centdiff)) # this outputs a SpatialPolygons object, I need a list of Polygons objecs so I can stick all newMCP polygons together into a SpatialPolygons object  
  
  return(newMCP@polygons) # returns a list of Polygons objects
  }
  
  
  ###----------------------------------------------------------------###
  # STEP 4: extract corine values for original stopover polygon and new randomly placed stopover polygon
  ###----------------------------------------------------------------###
  
  # NOTE: for MCPs close to the coast or with points over water, inevitably I will get a randomized polygon (or several) with no corine values (NAs), or with the polygon mainly over sea grid cells
  # if any of the randomMCPs are comprised of > 50% NAs or values 39-44 (intertidal, coastal, sea and ocean, etc), then repeat polygon randomization function and raster value extraction until all randomMCPs sample predominantly from land

  
  repeat {
    # turn original MCPs into newly shifted ones
    newMCPs <- sapply(MCP@polygons,getnewMCP) # don't use lapply, because that effectively makes it a nested list, use sapply so that the returned object is a single list of Polygons objects from function
    newMCPs2 <- SpatialPolygons(newMCPs)
    randomMCP <- SpatialPolygonsDataFrame(newMCPs2, MCP@data)
    proj4string(randomMCP) <- corine.crs
    
    # plot randomized stopover polygons with original stopover polygons
    #   windows(12,12)
    #   plot(MCP, col='black',xlim=c(MCPbox[[i]][1]-50000, MCPbox[[i]][2]+50000), ylim=c(MCPbox[[i]][3]-50000, MCPbox[[i]][4]+50000))
    #   plot(randomMCP, col='blue',add=TRUE,xlim=c(MCPbox[[i]][1]-50000, MCPbox[[i]][2]+50000), ylim=c(MCPbox[[i]][3]-50000, MCPbox[[i]][4]+50000))
    
    # NOT ADDING A BUFFER ANYMORE bufferMCP <- gBuffer(MCP, width=50, byid=TRUE) # buffer each MCP by 50m and create new SPDF object
    
    corine.origMCP <- extract(Europeraster, MCP)
    
    corine.randomMCP <- extract(Europeraster, randomMCP)
    
    # check that > 50% of corine values in each random stopover polygon are not from sea
    check <- logical()
    
    for (i in 1:length(corine.randomMCP)) {
      if ( (length(which(corine.randomMCP[[i]]==44 | corine.randomMCP[[i]]==39 | corine.randomMCP[[i]]==40 | corine.randomMCP[[i]]==41 | corine.randomMCP[[i]]==42 | corine.randomMCP[[i]]==43 | is.na(corine.randomMCP[[i]])))/length(corine.randomMCP[[i]]) < 0.5 ) ) {
        check[i] <- TRUE
      } else {
        check[i] <- FALSE
      }
    }
    
    if (all(check==TRUE)) break
    
  }

  ###----------------------------------------------------------------###
  # STEP 5: write habitat availability within the original stopover polygon and the random stopover polygon to a csv file
  ###----------------------------------------------------------------###
  
  # variables: mgroup number, corine.origMCP, corine.randomMCP
  
  mgroup <- levels(as.factor(stopovers$mgroup))
  
  for (i in 1:length(corine.origMCP)) {
    corine.origMCP[[i]] <- data.frame(name=rep(levels(stopovers$name),length(corine.origMCP[[i]])), mgroup=rep(mgroup[i],length(corine.origMCP[[i]])), corine.values=corine.origMCP[[i]])
  }

  for (i in 1:length(corine.randomMCP)) {
    corine.randomMCP[[i]] <- data.frame(name=rep(levels(stopovers$name),length(corine.randomMCP[[i]])), mgroup=rep(mgroup[i],length(corine.randomMCP[[i]])), corine.values=corine.randomMCP[[i]])
  }
  

  fullcorine.origMCP <- do.call(rbind,corine.origMCP)
  fullcorine.randomMCP <- do.call(rbind,corine.randomMCP)
  
  setwd(paste(datawd, "/Europe habitat availability/original stopover", sep=""))
  write.csv(fullcorine.origMCP, paste(levels(stopovers$name), " - original stopover.csv", sep=""), row.names=FALSE)
  
  setwd(paste(datawd, "/Europe habitat availability/random stopover", sep=""))
  write.csv(fullcorine.randomMCP, paste(levels(stopovers$name), " - random stopover.csv", sep=""), row.names=FALSE)
  
  habsummary.origMCP <- with(fullcorine.origMCP, table(corine.values,mgroup))
  habsummary.randomMCP <- with(fullcorine.randomMCP, table(corine.values,mgroup))
  
  # output a summary text file of habitat availability tables for original vs random stopover polygons
  setwd(paste(outputwd,"/Europe habitat use/",sep=""))
  
  cat(paste(levels(rdat$name), " : habitat availability", "\n\n", "Original stopover polygon", "\n\n", sep=""), file="habitat availability summary2.txt", sep="\n\n", append=TRUE)
  out<-capture.output(habsummary.origMCP)
  cat(out, file="habitat availability summary2.txt",sep="\n\n",append=TRUE)
  
  cat(paste("Random stopover polygon", "\n\n", sep="\n"), file="habitat availability summary2.txt", sep="\n\n", append=TRUE)
  
  out<-capture.output(habsummary.randomMCP)
  cat(out, file="habitat availability summary2.txt",sep="\n\n", append=TRUE)
  
  cat(c("########################################\n########################################\n########################################"), file="habitat availability summary2.txt", sep="\n\n", append=TRUE)
  
}

##################################################################
#
# Analysis 2: Proportional habitat use vs habitat availability in actual stopover polygon of
#
#
##################################################################

for (a in 1:29) {
  
  ### import resampled data + corine habitat use variables
  setwd(paste(datawd,resampwd,sep=""))
  resampledcorine <- list.files()[-8] # remove Idemili
  resampledcorine <- resampledcorine[-11] # remove Karma
  rdat <- read.csv(resampledcorine[a], header=T) # resampled data
  rdat$mgroup <- as.factor(rdat$mgroup)
  
  ### import habitat availability data, add corine use names
  setwd(paste(datawd, "/Europe habitat availability/original stopover", sep=""))
  availfiles <- list.files()
  avail <- read.csv(availfiles[a], header=T)
  avail$mgroup <- as.factor(avail$mgroup)
  
  # add names of corine land types to cuckoo dataset
  corine.name <- rep(NA,nrow(dat))
  corine.nameLABEL1 <- rep(NA,nrow(dat))
  corine.nameLABEL2 <- rep(NA,nrow(dat))
  
  for (i in 1:nrow(avail)){
    corine.name[i] <- as.character(clclegend[which(clclegend$GRID_CODE %in% avail[i,"corine.values"]), "LABEL4"])
    corine.nameLABEL1[i] <- as.character(clclegend[which(clclegend$GRID_CODE %in% avail[i,"corine.values"]), "LABEL1"])
    corine.nameLABEL2[i] <- as.character(clclegend[which(clclegend$GRID_CODE %in% avail[i,"corine.values"]), "LABEL2"])
  }
  
  newavail <- data.frame(avail, corine.name, corine.nameLABEL1, corine.nameLABEL2)
  
  # order land cover levels in both datasets so they appear in order of grid codes in barplot
  orderedcorine <- clclegend$LABEL4[which(clclegend$LABEL4 %in% levels(rdat$corine.name))] 
  rdat$corine.name <- factor(rdat$corine.name, orderedcorine, ordered=TRUE)
  
  orderedcorine <- clclegend$LABEL4[which(clclegend$LABEL4 %in% levels(newavail$corine.name))] 
  newavail$corine.name <- factor(newavail$corine.name, orderedcorine, ordered=TRUE)
  

  ### STEP 1: compare proportion of different habitats used with proportion of different habitats available, for each stopover
  # library(adehabitatHS), compana() to look at compositional analysis of habitat use, and widesIII() to look at selection ratios
  
  # subset out only stopovers
  rstopovers <- Subset(rdat, breedsite=="N")
  
  rprop <- as.data.frame.matrix(prop.table(with(rstopovers, table(mgroup, corine.name))))
  availprop <- as.data.frame.matrix(prop.table(with(newavail, table(mgroup, corine.name))))
  
  availmissing <- names(rprop)[which(!names(rprop) %in% names(availprop))]
  rmissing <- names(availprop)[which(!names(availprop) %in% names(rprop))]
  
  radd <- as.data.frame(matrix(0, nrow=nrow(rprop), ncol=length(rmissing)))
  colnames(radd) <- rmissing
  
  availadd <- as.data.frame(matrix(0, nrow=nrow(availprop), ncol=length(availmissing)))
  colnames(availadd) <- availmissing
  
  rprop2 <- data.frame(rprop, radd)
  availprop2 <- data.frame(availprop, availadd)
  
  rprop3 <- rprop2[,order(names(rprop2))]
  availprop3 <- availprop2[,order(names(availprop2))]
  
  colnames(rprop3) <- c(1:length(names(rprop3)))
  colnames(availprop3) <- c(1:length(names(availprop3)))
  
  rownames(rprop3) <- c(1:4)
  rownames(availprop3) <- c(1:4)
  
  
  
  ######################
  
  
  
  bymgroup <- rstopovers$mgroup
  splitstopovers <- split(rstopovers, list(bymgroup))
  
  bymgroup <- newavail$mgroup
  splitavail <- split(newavail, list(bymgroup))
  
  # create function to apply over each element of splitstopovers list, takes as arguments the habitat use dataset, and the habitat availabilitiy dataset
  function {
  use <- splitstopovers[[i]]
  use <- droplevels(use)
  
  habavail <- splitavail[[i]]
  habavail <- droplevels(habavail)
  
  rprop <- as.data.frame(prop.table(table(use$corine.name)))
  colnames(rprop) <- c("corine.name","freq")
  availprop <- as.data.frame(prop.table(table(habavail$corine.name)))
  colnames(availprop) <- c("corine.name","freq")
  
  omissing <- rprop[which(!rprop$corine.name %in% availprop$corine.name),"corine.name"]
  rmissing <- availprop[which(!availprop$corine.name %in% rprop$corine.name),"corine.name"]
  
  availprop2 <- rbind(availprop,data.frame(corine.name=omissing,freq=rep(0,length(omissing))))
  rprop2 <- rbind(rprop, data.frame(corine.name=rmissing, freq=rep(0,length(rmissing))))
  
  # reorder both so that corine.names are in same alphabetical order
  orderedcorinenames <- sort(as.character(availprop2$corine.name))
  
  availprop2$corine.name <- factor(availprop2$corine.name, orderedcorinenames, ordered=TRUE)
  availprop3 <- availprop2[order(availprop2$corine.name),]
  
  rprop2$corine.name <- factor(rprop2$corine.name, orderedcorinenames, ordered=TRUE)
  rprop3 <- rprop2[order(rprop2$corine.name),]
  
  diff <- data.frame(corine.name=rprop3$corine.name, difference=abs(rprop3$freq - availprop3$freq))
  
  summary <- data.frame(corine.name=rprop3$corine.name, used=rprop3$freq, available=availprop3$freq, difference=(rprop3$freq-availprop3$freq))
  
  rprop3$freq*300
  
  chisqout <- chisq.test(rprop3$freq*300, p=availprop3$freq)
  chisqout <- chisq.test(diff$difference)
  
}