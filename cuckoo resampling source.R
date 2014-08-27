##########################################################
#
#  Cuckoo data - tidying source file, STOPOVER point resampling according to error distribution, and adding new variables (distance moved between points, bearing, movement groups, and movement types)
#
#	Samantha Franks
#	15 Nov 2013
# 28 Nov 2013
# 17 Dec 2013 - reran resampling and new code for resampled mgroups, mtypes based on original data trimmed down to just C and S mtypes (leave out migration movements to avoid resampling migration movements that occur within a transmission cycle)
#
##########################################################

rm(list=ls())

library(NCStats)
library(rgdal)
library(maptools)
library(shape)
library(splancs)
library(geosphere)
library(adehabitat) # 
library(sp) # converts R dataframes into spatially explicit dta
library(rgeos)
library(raster)
library(adehabitatHR)


###---------------------------------------------------------------###
#         Load cuckoo data from original source, trim out "M" type movements
###---------------------------------------------------------------###

setwd("C:/Users/samf/Documents/Git/cuckoos/scripts")
source("cuckoo original data tidying source.R")

# creates object called fulloriginal, which is the dataset with distance, bearing, mgroups, and mtypes added

setwd("C:/Users/samf/Documents/Git/cuckoos/data")
write.csv(fulloriginal, "cleaned original data.csv", row.names=F)

### --- Subset out migration movements --- ###


dat <- Subset(fulloriginal, mtype=="C" | mtype=="S" | is.na(mtype))

###-----------------------------------------------###
#         Point location resampling
###-----------------------------------------------###

# resample from UTM points so units are all in metres

# STEP 1: Within each transmission cycle, there will be a number of points with associated elliptical uncertainty.  Resample 1 point from the available points of the transmission cycle such that the probability of a point being sampled is weighted by 1/(area of the uncertainty ellipse).

# area of an ellipse = pi * r1 * r2, where r1 = length of semi-major axis (ie 1/2 the major axis) and r2 = length of the semi-minor axis

# STEP 2: Repeat Step 1 for each transmission cycle that is part of the stopover.

# STEP 3: For each point selected from each transmission cycle, randomly sample one point from the error ellipse of each observation selected.

# STEP 4: Draw an MCP or take the 95% kernel around that sample of points

# STEP 5: Repeat Step 1 + 2 some number of times (e.g. 100, 1000, etc), and then overlay the MCPs or kernels and the habitat that is considered "used" by a cuckoo will be everything defined by the overlapping MCPs or the 95% kernel.

# create new list of dataframes, with n number of elements (n is the number of resampling occasions, e.g. 100)
# each element of list has cuckoo names, tcycle numbers, and one set of resampled lats/longs

byname <- dat$name

newdat <- split(dat, list(byname))

nbootstraps <- 100 # SET NUMBER OF BOOTSTRAP RESAMPLINGS FOR EACH POINT FROM THE ERROR ELLIPSES ON THE POINT

allbird <- list(ls) # list in which to store resampled, tcycle-selected data for each bird, where a list element is a bird

totaltime <- system.time({
  
  for (a in 1:length(newdat)){ # BEGIN LOOP 1, where "a" counts through the number of list elements (ie. the number of individual birds)
    
    # for each newdat[[i]], count the number of levels of tcycle
    
    ind <- newdat[[a]] # counter a - number of elements in list (ie. number of cuckoo names)
    
    ### Convert long/lat coordinates into wGS 84/World Mercator projection coordinates (units are metres), in order to be able to sample from error ellipses, which are in metres
    # epsg: 3395 - +proj=merc +ellps=wGS84 +datum=WGS84
    
    coordinates(ind) <- c("long", "lat")
    
    # set projection and datum information for cuckoo data
    proj4string(ind) <- CRS("+proj=longlat +datum=WGS84")
    # summary(dat)
    
    # Transform cuckoo coordinates to UTM using spTransform() in library(rgdal)
    
    temp <- spTransform(ind, CRS("+init=epsg:3395"))
    # summary(temp)
    
    indbird.resampledlocations <- list()
    
    bootstraptime1 <- system.time({
      
      for (b in 1:nbootstraps) { # BEGIN LOOP 2, where "b" counts the bootstrap resamplings per individual bird. System time for 100 bootstraps = 
        
        ############ Can't get this code using the runifdisc/affine functions in library(spatstat) to produce the right size and shape of polygon
        # resample random point from ellipse for each row of temp - each point gets written to [[1]] of indbird.resampledlocations
        # random points in the ellipse with centre (x,y - lat,long) and major and minor axis
        
        
        #   for (i in 1:nrow(temp)){
        #   x <- runifdisc(42, centre=c(coordinates(temp)[i,2], coordinates(temp)[i,1]))
        #   y <- affine.ppp(x, mat=diag(c(temp$error.major[i]*2, temp$error.minor[i]*2)), unitname=metre))
        #   }
        #   
        
        ########### NOT USING THIS METHOD - create a spatial polygon ellipse shape, then generate coordinates with random long/lats bounded by max and min of the polygon ellipse
        # test whether point falls within ellipse
        # if yes, then select point
        # if no, then choose a new random point, and test again
        # test until point falls inside ellipse
        
        #   errorellipse <- getellipse(rx=temp$error.major[1], ry=temp$error.minor[1], mid=c(coordinates(temp)[1,1], coordinates(temp)[1,2]))
        #   
        #   ellpoly <- SpatialPolygons(list(Polygons(list(Polygon(errorellipse)),"errorellipse")))
        #   plot(ellpoly)
        
        ########### NEW METHOD - create ellipse polygon shape using getellipse() in library(shape), then use csr() in library(splancs) to generate a completely spatially random point within the polygon ellipse
        
        # for every observation in temp (ie. the dataset pertaining to the cuckoo of interest), create a new observation with long/lat position based on a resampled point within the error ellipse of the original observation
        
        # create two new variables called newlong and newlat in a dataframe
        # write each randomly generated point within the error ellipse of each observation to this new dataframe
        
        newlongs <- rep(NA,nrow(temp))
        newlats <- rep(NA,nrow(temp))
        
        for (i in 1:nrow(temp)){ # loop counts through each row in temp
          
          ###??? add in the orientation of the ellipse
        errorellipse <- getellipse(rx=temp$error.major[i], ry=temp$error.minor[i], mid=c(coordinates(temp)[i,1], coordinates(temp)[i,2]), angle=temp$ell.orientation[i])
          
          randpts <- csr(errorellipse,1)
          newlongs[i] <- randpts[1]
          newlats[i] <- randpts[2]
          
        }
  
    
          indbird.resampledlocations[[b]] <- data.frame(temp, newlongs, newlats) # list element for this variable will count up to 100 or however many bootstrap resamplings there are per individual
        
        
      } # END LOOP 2 - bootstrap resamplings "b"
      
    })
    
    bootstraptime1 <- c(bootstraptime1, nbootstraps)
    
    indbird.resampled.tcycleselect <- list()
    
    for (c in 1:length(indbird.resampledlocations)) { # BEGIN LOOP 3 - FOR EACH RESAMPLED LIST FROM LOOP 2, PICK ONE POINT FROM EACH TRANSMISSION CYCLE
      # feed results of LOOP 2 into LOOP 3 which picks one point for each transmission cycle
      
      temp3 <- indbird.resampledlocations[[c]]
      
      # below code samples one point for each transmission cycle (l = number of tcycles in list element) in the list element
      # sample probability is proportional to 1/area of the ellipse
      # builds a new list element (where an element = an individual cuckoo) with those sampled points for every transmission cycle
      
      l <- levels(droplevels(temp3$tcycle)) # number of tcycles in selected list element
      
      selectedeventrowid <- rep(NA,length(l)) # create blank vector in which to store the row id positions of selected events from each transmission cycle
      
      for (i in 1:length(l)){ # counts through number of t cycles dataset
        
        tcyclesample <- temp3[temp3$tcycle==l[i],] # identifies transmission cycle on which to apply resampling
        
        selectedeventrowid[i] <- ifelse(nrow(tcyclesample)==1, tcyclesample$id, with(tcyclesample,sample(tcyclesample$id, 1, replace=T, prob=1/(pi*tcyclesample$error.major*tcyclesample$error.minor))))
        
        # for all data in temp corresponding to the row of selectedeventrowid, write to a new dataframe which will compile the points used for every individual
        
      }
      
      indbird.resampled.tcycleselect[[c]] <- temp3[which(temp3$id %in% selectedeventrowid),]
      
    } # END LOOP 3 - indbird.resampled.tcycleselect is a list c elements long, where c = # of bootstraps per individual.  Each element of the list is comprised of one selected point per transmission cycle, where the point was randomly generated from within the error ellipse of the original point, and where the point was selected from all the points in a transmission cycle based on a probability of 1/(area of the error ellipse)
    
    # once last element (i = 100) is written to list, then unsplit the list (if possible) to reconcatenate rows into a single dataframe for that individual
    
    allbird[[a]] <- do.call(rbind, indbird.resampled.tcycleselect)
    
    #???????? write dataframe to a .csv to be used later
    
  } # END LOOP 1 - loop that counts through each list-subsetted individual bird
  
})

summary(allbird)


# save output in reusable format
# 100 boostraps for all raw data and all birds took 6675.55 seconds

allbirdfinal <- do.call(rbind, allbird) # rbind all list elements

setwd("C:/Users/samf/Documents/Git/cuckoos/data")

write.csv(allbirdfinal, "allbirdfinal - resampling.csv", row.names=FALSE) # write entire file to csv

# write each individual's bootstrapped data to an individual csv to be called again later

setwd("C:/Users/samf/Documents/Git/cuckoos/data/resampled data with uncertainty 100 bootstraps")

for (i in 1:length(allbird)){
  write.csv(allbird[[i]], paste(allbird[[i]]$name[1],".csv", sep=""), row.names=FALSE)
  }

# write data summary file to .csv
allsummary <- data.frame(name=rep(NA,length(allbird)), no.samples=rep(NA, length(allbird)))

for (i in 1:length(allbird)){
  allsummary[i,] <- c(as.character(allbird[[i]]$name[1]), nrow(allbird[[i]]))
}

setwd("C:/Users/samf/Documents/Git/cuckoos/data")
write.csv(allsummary, "summary of bootstrap resamplings.csv", row.names=FALSE)


####################################################################
####################################################################
####################################################################
#############     DO NOT NEED SECTION BELOW IF MGROUP & MTYPE ALREADY INCLUDED         ###################################################
####################################################################
####################################################################
####################################################################


###############################################################
#
# ADD DISTANCES, BEARINGS, AND MOVEMENT TYPES TO RESAMPLED DATA
#
###############################################################

###----------------------------------------------------------
#
#   NOTE
#
#   distances and bearings in following code use the original location points for a tcycle for calculations
#
###----------------------------------------------------------

newresample <- list()

for (a in 1:31) {
  
  setwd("C:/Users/samf/Documents/Git/cuckoos/data/resampled data with uncertainty 100 bootstraps")
  resamplefiles <- list.files()
  
  resample <- read.csv(resamplefiles[a], header=T)
  
  # removes observations with really large errors
  
  resample <- Subset(resample, error.major <= 5000)
  
  setwd("C:/Users/samf/Documents/Git/cuckoos/data/original data + distance, movement groups, etc")
  originalfiles <- list.files()
  
  original <- read.csv(originalfiles[a], header=T)
  
  # for each level of tcycle in original data, look up the corresponding mgroup and mtype values
  # code extracts the position number of each new tcycle, then looks up the corresponding mgroup and mtype value in that position in original, then merging with resample according to tcycle to add mgroup and mtype 
  firsts <- rep(NA,length(levels(as.factor(original$tcycle)))) # blank variable to hold row numbers of first row of data for each individual
  
  for (j in 1:length(levels(as.factor(original$tcycle)))){ # loop counts through tcycle levels
    firsts[j] <- min(which(original$tcycle==levels(as.factor(original$tcycle))[j])) # row number of first observation of each tcycle
  }
  
  orig.mvmts <- data.frame(tcycle=original$tcycle[firsts], mgroup=original$mgroup[firsts], mtype=original$mtype[firsts])
  
  newresample[[a]] <- merge(resample,orig.mvmts,by="tcycle",all=TRUE)
   
}


# write each individual's bootstrapped data with new variables to an individual csv to be called again later

setwd("C:/Users/samf/Documents/Git/cuckoos/data/resampled data with uncertainty 100 bootstraps + movement groups and types")

for (i in 1:length(newresample)){
  write.csv(newresample[[i]], paste(newresample[[i]]$name[1],".csv", sep=""), row.names=FALSE)
}




###################################################
###################################################
###################################################
###################################################
###################################################


###----------------------------------------------------------
#
#   NOTE
#
#   distances and bearings in following code use the centroid of the points for a tcycle for calculations
#
###----------------------------------------------------------


# add distance/bearing variables for distance moved between previous and current location, and bearing taken to get there
# first location will be NA (ie. no previous location, because the first location for an individual is the capture event)

# for all rows of each list element of allbird, create new variables: distance, bearing - calculate these based on the centroid of the resampled points for each tcycle

# remove reference points

# remove points where semi-major error ellipse is > 5km (5000m)

newallbird <- list()

nbirds <- 31 #number of individuals

totaltime <- system.time({
  
  
  # could probably write this for loop as a function, and lapply it over the elements of a list containing the resampled datafiles for each bird (to make things faster if it needs tweaking)
  
for (a in 1:nbirds){
  
  ### --- STEP 1: read in data file for bird, convert to appropriate coordinate system,  --- ###
  
  setwd("C:/Users/samf/Documents/Git/cuckoos/data/resampled data with uncertainty 100 bootstraps")
  cuckoofiles <- list.files()
  dat <- read.csv(cuckoofiles[a], header=T)
  
  # order dat according to transmission cycle
  datorder <- dat[order(dat$tcycle),]
  
  # convert long/lat of original locations to appropriate coordinate system
  coordinates(datorder) <- c("long","lat")
  
  proj4string(datorder) <- CRS("+init=epsg:3395") # set projection and datum information for cuckoo data, which are in wGS 84/World Mercator projection coordinates (set earlier in script)
  neworder <- spTransform(datorder, CRS("+proj=longlat +datum=WGS84")) # transform to long/lats which are required for the distance calculation
  
  neworder <- as.data.frame(neworder)

  neworder <- Subset(neworder, long > -19) # remove reference points in North America
  
  smallerror <- Subset(neworder, error.major <= 5000)
  
  # apply the distance/bearing etc calculation between tcycles: calculate the mean (centroid) point of the polygon described by the resampled points in tcycle "n", and calculate the distance between that centroid "n" and the centroid of the tcycle "n-1" to calculate average distance moved between tcycles
  
  # transform newlongs/newlats to coordinates
  coordinates(smallerror) <- c("newlongs","newlats")
  
  proj4string(smallerror) <- CRS("+init=epsg:3395") # set as World Mercator/WGS84
  
  
  
  
  ### --- STEP 2: CALCULATE CENTROIDS OF MCPs OF TCYCLES --- ###
  
  # use mcp() to calculate MCP around resampled points of the tcycle
  
  smallerror$tcycle <- as.factor(smallerror$tcycle)
  tcycles <- smallerror[,c("tcycle")]
  MCP <- mcp(tcycles, percent=100)
  
  extractcentroid <- function(spdf){ #function to extract centroid values from the polygons (list components) of a SpatialPolygonsDataFrame
    centroid <- spdf@labpt
    return(centroid)
  }
  
  extractid <- function(spdf){ #function to extract id name (tcycle no.)
    id <- spdf@ID
    return(id)
  }
  
  centroids <- t(sapply(MCP@polygons,extractcentroid)) # represents centroids of MCPs describing the resampled points of each tcycle
  ids <- t(sapply(MCP@polygons,extractid)) # represents centroids of MCPs describing the resampled points of each tcycle
  
  centroids <- as.data.frame(centroids)
  ids <- as.factor(ids)
  
  centroid.id <- data.frame(ids,centroids)
  colnames(centroid.id) <- c("tcycleid","centroidlong","centroidlat")
  
  
  
  ### --- STEP 3: DISTANCE & BEARING --- ###
  
  # calculate distances and bearings based on centroid longs/lats for each tcycle
  
  # distbearmvmt <- function(input){ # requires 2 variables called newlongs & newlats or long & lat or centroidlong & centroidlat
  # use newlongs & newlats to calculate distances based on resampled points
  # use long & lat to calculate distances based on original measured location point
  

  #dist.cuckoo <- input[,c("newlongs","newlats")]
  #coordinates(dist.cuckoo) <- c("newlongs", "newlats")
  #dist.cuckoo <- input[,c("long","lat")]
  #coordinates(dist.cuckoo) <- c("long", "lat")
  dist.cuckoo <- centroid.id[,c("centroidlong","centroidlat")] # change data.frame name to "input" if using as a function
  coordinates(dist.cuckoo) <- c("centroidlong", "centroidlat")
  proj4string(dist.cuckoo) <- CRS("+init=epsg:3395") # set coordinates as long/lats which are required for the distance calculation
  
  dist.cuckoo <- spTransform(dist.cuckoo, CRS("+proj=longlat +datum=WGS84"))
  
  d1 <- dist.cuckoo@coords
  d2 <- dist.cuckoo@coords[-1,] # create new matrix minus the first row so that d2 starts at the second observation
  last.d2 <- matrix(c(0,0), nrow=1, dimnames=list(1, c("centroidlong","centroidlat"))) # create a placeholder last row in the second distance matrix
  d2 <- rbind(d2, c(0,0))
  
  dist <- distCosine(d1,d2)/1000 # distance between points, in km
  bear <- bearing(d1,d2) # bearing between points
  dist <- c(NA,dist[-length(dist)])
  bear <- c(NA,bear[-length(bear)])
  
  distbear <- data.frame(centroid.id, distance=dist, bearing=bear)
  
  
  
  ### --- MOVEMENT GROUPS --- ###
  
  mgroup <- rep(NA,nrow(distbear))
  
  for (n in 1:nrow(distbear)){
    if (is.na(distbear$distance[n])) {
      mgroup[n] <- 1
    } else if (distbear$distance[n] <= 30) {
      mgroup[n] <- mgroup[n-1]
    } else {mgroup[n] <- mgroup[n-1] + 1}
  }
  
  ### --- MOVEMENT TYPES --- ###
  
  mtype <- c("C", rep(NA, nrow(distbear)-1))
  
  for (n in 2:(nrow(distbear)-1)){
    if (mgroup[n] != mgroup[n-1] & mgroup[n] != mgroup[n+1]) {
      mtype[n] <- "M"
    } else {
      mtype[n] <- "S"
    }
  }
  
  completedata <- data.frame(distbear,mgroup,mtype)
  # return(completedata)
  
  # } remove # if using this section of code as a function
  
  migvar <- completedata
  
  # add centroid dataframe to rest of bird info such that centroids[1,] is cbinded to all rows of smallerror where tcycle=1, etc (use merge to merge by tcycle)
  
  mergedat <- merge(smallerror, migvar, by.x="tcycle", by.y="tcycleid") # newlong/newlat and centroidlong/centroidlat are in EPSG:3395, long/lat are in +longlat WGS84
  
  newallbird[[a]] <- mergedat
  
}

})



# write data files with new variables (distance, movement groups, etc) to .csv - one large with all birds, and separate files for each bird

new.allbirdfinal <- do.call(rbind, newallbird) # rbind all list elements

setwd("C:/Users/samf/Documents/Git/cuckoos/data/")
write.csv(new.allbirdfinal, "allbirdfinal - resampling with distances and movement groups.csv", row.names=FALSE) # write entire file to csv

# write each individual's bootstrapped data to an individual csv to be called again later

setwd("C:/Users/samf/Documents/Git/cuckoos/data/resampled data with uncertainty 100 bootstraps + distance, movement groups, etc")

for (i in 1:length(newallbird)){
  write.csv(newallbird[[i]], paste(newallbird[[i]]$name[1],".csv", sep=""), row.names=FALSE)
}

# long/lat are in CRS("+proj=longlat +datum=WGS84")
# newlongs/newlats are in CRS("+init=epsg:3395")

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

# ###-----------------------------------------------###
# #         Load resampled datasets to add distances, bearings, and movement types
# ###-----------------------------------------------###

# ###########################################################
# #
# #   NOTE
# #
# #   distances and bearings in following code are calculated based on resampled values for each resampling
# #
# ###########################################################

# ### --- Load data --- ###

# setwd("C:/Users/samf/Documents/Git/cuckoos/data/")
# dd <- read.csv("allbirdfinal - resampling.csv", header=T)

# byname <- dd$name

# allbird <- split(dd, list(byname))


# ### ----------- Calculate distance moved, bearing, movement types, movement groups ---------- ###

# ### --- Distance + Bearing --- ###

# # add distance/bearing variables for distance moved between previous and current location, and bearing taken to get there
# # first location will be NA (ie. no previous location)

# # for all rows of each list element of allbird, create new variables: distance, bearing - calculate these based on the original lat/long variables

# # remove reference points

# # remove points where semi-major error ellipse is > 5km (5000m)

# newallbird <- list()

# for (i in 1:length(allbird)){
  
  # ### --- STEP 1: split list according to resample no. --- ###
  # dat <- allbird[[i]]
  
  # # remove reference points and any other weird points outside of the range we'd expect
  
  # # use newlongs & newlats if we want to calculate new variables based on resampled points
  # # use long & lat to calculate new variables based on original location points
  
  # #coordinates(dat) <- c("newlongs","newlats")
  # coordinates(dat) <- c("long","lat")
  
  # proj4string(dat) <- CRS("+init=epsg:3395") # set projection and datum information for cuckoo data, which are in wGS 84/World Mercator projection coordinates (set earlier in script)
  # newdat <- spTransform(dat, CRS("+proj=longlat +datum=WGS84")) # transform to long/lats which are required for the distance calculation
  
  # newdat <- as.data.frame(newdat)
  # newdat <- Subset(newdat, long > -19) # remove reference points in North America
  
  # # add resample no. variable to be able to split dataset according to sample number
  # newdat$tcycle <- as.factor(newdat$tcycle)
  # l <- droplevels(newdat$tcycle)
  # resample.no <- c()
  
  # for (j in 1:100){
    # x <- rep(j,length(levels(l)))
    # resample.no <- c(resample.no, x)
  # }
  
  # newdat <- data.frame(newdat, resample.no)
  # newdat$resample.no <- as.factor(resample.no)
  
  # # remove observations where semi-major axis error is > 5000km
  # newdat.smallerror <- Subset(newdat, error.major <= 5000)
  
  # byresampleno <- newdat.smallerror$resample.no
  # datbysample <- split(newdat.smallerror, byresampleno)

  # ### --- STEP 2: lapply function distbearmvmt across all 100 resamples in list datbysample, to calculate distance moved, bearing, and movement groups and types between tcycles --- ###
  
  # ### --- FUNCTION to measure distance moved between each transmission cycle, bearing direction between each tcycle, and whether movements between tcycles were stopovers or migratory --- ###
  
  # distbearmvmt <- function(input){ # requires 2 variables called newlongs & newlats or long & lat
    # # use newlongs & newlats to calculate distances based on resampled points
    # # use long & lat to calculate distances based on original measured location point
    
    # ### DISTANCE & BEARING
    # #dist.cuckoo <- input[,c("newlongs","newlats")]
    # #coordinates(dist.cuckoo) <- c("newlongs", "newlats")
    # dist.cuckoo <- input[,c("long","lat")]
    # coordinates(dist.cuckoo) <- c("long", "lat")
    # proj4string(dist.cuckoo) <- CRS("+proj=longlat +datum=WGS84") # set coordinates as long/lats which are required for the distance calculation
    
    # d1 <- dist.cuckoo@coords
    # d2 <- dist.cuckoo@coords[-1,] # create new matrix minus the first row so that d2 starts at the second observation
    # last.d2 <- matrix(c(0,0), nrow=1, dimnames=list(1, c("newlongs","newlats"))) # create a placeholder last row in the second distance matrix
    # d2 <- rbind(d2, c(0,0))
    
    # dist <- distCosine(d1,d2)/1000 # distance between points, in km
    # bear <- bearing(d1,d2) # bearing between points
    # dist <- c(NA,dist[-length(dist)])
    # bear <- c(NA,bear[-length(bear)])
    
    # distbear <- data.frame(input, distance=dist, bearing=bear)
    
    # ### MOVEMENT GROUPS
    # mgroup <- rep(NA,nrow(distbear))
    
    # for (n in 1:nrow(distbear)){
      # if (is.na(distbear$distance[n])) {
        # mgroup[n] <- 1
      # } else if (distbear$distance[n] <= 100) {
        # mgroup[n] <- mgroup[n-1]
      # } else {mgroup[n] <- mgroup[n-1] + 1}
    # }
    
    # ### MOVEMENT TYPES
    # mtype <- c("C", rep(NA, nrow(distbear)-1))
    
    # for (n in 2:(nrow(distbear)-1)){
      # if (mgroup[n] != mgroup[n-1] & mgroup[n] != mgroup[n+1]) {
        # mtype[n] <- "M"
      # } else {
        # mtype[n] <- "S"
      # }
    # }
    
    # completedata <- data.frame(distbear,mgroup,mtype)
    # return(completedata)
    
  # }
  
  # dat.allvar <- lapply(datbysample, distbearmvmt)
  # dat.allvar.onebird <- do.call(rbind, dat.allvar)
  # newallbird[[i]] <- dat.allvar.onebird
  
# }


# # write data files with new variables (distance, movement groups, etc) to .csv - one large with all birds, and separate files for each bird

# new.allbirdfinal <- do.call(rbind, newallbird) # rbind all list elements

# write.csv(new.allbirdfinal, "allbirdfinal - resampling with distances and movement groups.csv", row.names=FALSE) # write entire file to csv

# # write each individual's bootstrapped data to an individual csv to be called again later

# setwd("C:/Users/samf/Documents/Git/cuckoos/data/resampled data with uncertainty 100 bootstraps + distance, movement groups, etc")

# for (i in 1:length(newallbird)){
  # write.csv(newallbird[[i]], paste(newallbird[[i]]$name[1],".csv", sep=""), row.names=FALSE)
# }

# # long/lat are in CRS("+proj=longlat +datum=WGS84")
# # newlongs/newlats are in CRS("+init=epsg:3395")
    
# ######################################################
# ######################################################
# ######################################################
# ######################################################
# ######################################################

  
# # ### --- Movement group + Movement type --- ###
# # 
# # # for each movement, classify movement type as migratory or stopover
# # # if distance < 25km, movement is a stopover (S), if >25km, movement is migration (M)
# # # ?????????? work out what the best distance is to use for categorising movements - somewhere between 20-30km?
# # # for every group of movements (ie. SSSSMMMSSS are 3 groups), give those observations the same movement group number
# # 
# # ### Histogram of cuckoo distance movements
# # 
# # hist(cuckoo.use$distance)
# # 
# # par(mfrow=c(2,2))
# # hist(subset(cuckoo.use, distance < 100)$distance, main=c("d < 100km"), xlab=c("Distance (km)"))
# # hist(subset(cuckoo.use, distance > 100 & distance < 500)$distance, main=c("100km < d < 500km"), xlab=c("Distance (km)"))
# # hist(subset(cuckoo.use, distance > 500)$distance, main=c("d > 500km"), xlab=c("Distance (km)"))
# # hist(subset(cuckoo.use, distance > 20 & distance < 200)$distance, main=c("20km < d < 200km"), breaks=20, xlab=c("Distance (km)"))
# # 
# # 
# # # use <= 30km as the threshold for distinguishing between a stopover (S) and a migration (M), and between movement groups (e.g. 2 stopovers that might not be separated by a migration)
# # 
# # ### Movement numbers
# # 
# # # make a blank variable in which to put movement groups
# # mgroup <- rep(NA,nrow(cuckoo.use))
# # 
# # for (n in 1:nrow(cuckoo.use)){
# #   if (is.na(cuckoo.use$distance[n])) {
# #     # or alternatively, for loop counts from 2:nrow... and if(is.na(cuckoo.use$distance[n]))
# #     mgroup[n] <- 1
# #   } else if (cuckoo.use$distance[n] <= 30) {
# #     mgroup[n] <- mgroup[n-1]
# #   } else {mgroup[n] <- mgroup[n-1] + 1}
# # }
# # 
# # ###################################################
# # 
# # ##################### OLD CODE FOR MOVEMENT TYPE ##########################
# # # create new variable for movement type
# # # if distance is NA, then this is the capture occasion ("C")
# # # if bearing is NA, then could be capture occasion or bird didn't move between subsequent observations
# # 
# # # mtype <- ifelse(cuckoo.use$distance <= 30, c("S"), c("M"))
# # # mtype[is.na(mtype)] <- c("C")
# # # mtype <- as.factor(mtype)
# # 
# # ####################################################
# # 
# # ######### OLD CODE FOR MOVEMENT GROUP USING "M" and "S" QUALIFIERS, WHICH WILL GIVE INCORRECT RESULTS IF TWO STOPOVERS IN SEQUENCE ARE IN DIFFERENT LOCATIONS #########
# # # for each individual cuckoo (call using $name), get row number [x] of first observation for that bird (the capture observation) - that is starting point, and mgroup[x] <- 1, then proceed into for loop, where i counts (x+1):nrow(cuckoo$name[individual name])
# # # 
# # # for (j in 1:31){ # loop counts through cuckoo names
# # #   a <- min(which(cuckoo.use$name==levels(cuckoo.use$name)[j])) # row number of first observation of cuckoo j
# # #   b <- max(which(cuckoo.use$name==levels(cuckoo.use$name)[j])) # row number of last observation of cuckoo j
# # # 
# # #   mgroup[a] <- 1  # set capture occations as movement group 1
# # #   
# # #   for(i in (a+1):b){  # loops through rows of each cuckoo
# # #     if(mtype[i] == c("S") & mtype[i-1] == c("S")){          # stopover
# # #       mgroup[i] <- mgroup[i-1]
# # #     } else if(mtype[i] == c("S") & mtype[i-1] == c("M")){   # stopover
# # #       mgroup[i] <- mgroup[i-1]
# # #     } else if(mtype[i] == c("M") & mtype[i-1] == c("S")){   # migration
# # #       mgroup[i] <- mgroup[i-1] + 1
# # #     } else if(mtype[i] == c("M") & mtype[i-1] == c("M")){   # migration
# # #       mgroup[i] <- mgroup[i-1] + 1
# # #     } else if(mtype[i] == c("M") & mtype[i-1] == c("C")){   # if migration directly follows the capture occasion
# # #       mgroup[i] <- mgroup[i-1] + 1
# # #     } else {
# # #       mgroup[i] <- mgroup[i-1]                          # stopover, if capture is followed by "S"
# # #     }
# # #   }
# # #   
# # # }
# # 
# # ####################################################
# # 
# # # add a variable, mtype2, to be able to subset data to include only movement groups where the # of rows for that group number are > 2
# # 
# # # for each cuckoo, if mgroup value is the same as the one before it and after it, then mtype2 = S, otherwise mtype2 = M, EXCEPT for the first observation for every cuckoo, which was the capture occasion
# # # aka, if a bird has > 2 observations in a single movement group, then it constitutes as a stopover, even if the first observation was categorised as a migration (which brought it to that stopover)
# # # in order to change the number of consecutive locations that is required to equal a stopover from 2 to something else, then will need to modify if/else statement
# # 
# # firsts <- rep(NA,31) # blank variable to hold row numbers of first row of data for each individual
# # lasts <- rep(NA,31) # blank variable to hold row numbers of last row of data for each individual
# # 
# # for (j in 1:31){ # loop counts through cuckoo names
# #   firsts[j] <- min(which(cuckoo.use$name==levels(cuckoo.use$name)[j])) # row number of first observation of cuckoo j
# #   lasts[j] <- max(which(cuckoo.use$name==levels(cuckoo.use$name)[j])) # row number of last observation of cuckoo j
# # }
# # 
# # # loop to create mtype2 based on mgroup variable
# # mtype <- rep(NA, nrow(cuckoo.use))
# # 
# # for (n in 1:nrow(cuckoo.use)){
# #   
# #   if (n %in% firsts == TRUE) {
# #     mtype[n] <- "C"
# #   } else {
# #     if (mgroup[n] != mgroup[n-1] & mgroup[n] != mgroup[n+1]) {
# #       mtype[n] <- "M"
# #     } else {
# #       mtype[n] <- "S"
# #     }
# #   }
# # }
# # 
# # newcuckoo <- data.frame(cuckoo.use, mgroup, mtype)
# # newcuckoo[1:999,c("name","mgroup","mtype")]     # check output
# # 
# # write.csv(newcuckoo, "cuckoo movements for analysis.csv", row.names=FALSE)



