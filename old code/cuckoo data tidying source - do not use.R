##########################################################
#
#  Cuckoo data - tidying source file and point resampling according to error distribution
#
#	Samantha Franks
#	15 Nov 2013
#
#
##########################################################

library(NCStats)
library(sp)
library(rgdal)
library(maptools)
#library(spatstat)
library(shape)
library(splancs)
library(geosphere)

###-----------------------------------------------###
#         Load cuckoo data & tidying
###-----------------------------------------------###

# cuckoo coordinate locations are in Lat-Long and UTM coordinates and WGS84 datum

setwd("C:/Users/samf/Documents/R/projects/cuckoos/data")
#d <- read.csv("BTO Cuckoo migration study - raw data + new variables 15112013.csv", header=T) # cuckoo movement

d <- read.csv("cuckoos raw data 20131118_reformatted.csv", header=T)

# cuckoo <- read.csv("cuckoo movements.csv", header=T) # cuckoo movement
# 
# reqcuckoo <- with(cuckoo, data.frame(timestamp, julian=julian.date, year, month, day, id, lat, long, duration, days.since.tagging, name, country, continent, movement.number))

# $filter refers to the values identified as best in the transmission cycle by the Douglas filter applied in Movebank: TRUE values are the best locations
# $loc.class refers to the Argos location class: from best to worst, 3,2,1,0,A,B,Z. Alpha location classes have no associated errors, so remove. Location class=0 have gigantic major axis errors, so also remove these. Use only location classes 1:3
# error.major = error radius along the major axis of the ellipse
# error.minor = error radius along the minor axis of the ellipse

# d2 <- with(d, data.frame(name=individual.local.identifier, id=event.id, filter=visible, tcycle=Transmission.cycle, timestamp, julian, year, month, day, long=location.long, lat=location.lat, loc.class=argos.lc, error.major=argos.semi.major, error.minor=argos.semi.minor))

d2 <- with(d, data.frame(name, id, tcycle=transmission.cycle, timestamp, julian, year, month, day, hour, long=longitude, lat=latitude, loc.class=location_class, error.major=semi_major_axis, error.minor=semi_minor_axis, ell.orientation=ellipse_orientation, sensor1=X_1, sensor2=X_2, sensor3=X_3, sensor4=X_4))

# d2$filter <- as.factor(d2$filter) # define the filter as a factor rather than logical

d2$tcycle <- as.factor(d2$tcycle)

d3 <- Subset(d2, loc.class=="1" | loc.class=="2" | loc.class=="3") # subset by best location classes

d4 <- na.omit(d3) # complete cases only

d5 <- Subset(d4, name!="" & name!="115592" & name!="115601" & name!=121792)

d6 <- Subset(d5, error.major>0 & error.minor>0)

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

dat <- d6

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
          # convert semi-major and semi-minor axes of error ellipse measures into radians
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

write.csv(allbirdfinal, "allbirdfinal - resampling.csv", row.names=FALSE) # write entire file to csv

# write each individual's bootstrapped data to an individual csv to be called again later

setwd("C:/Users/samf/Documents/R/projects/cuckoos/data/resampled data with uncertainty 100 bootstraps")

for (i in 1:length(allbird)){
  write.csv(allbird[[i]], paste(allbird[[i]]$name[1],".csv", sep=""), row.names=FALSE)
  }

# write data summary file to .csv
allsummary <- data.frame(name=rep(NA,length(allbird)), no.samples=rep(NA, length(allbird)))

for (i in 1:length(allbird)){
  allsummary[i,] <- c(as.character(allbird[[i]]$name[1]), nrow(allbird[[i]]))
}

write.csv(allsummary, "summary of bootstrap resamplings.csv", row.names=FALSE)

###-----------------------------------------------###
#         Load resampled datasets to add distances, bearings, and movement types
###-----------------------------------------------###

### --- Load data --- ###



### ----------- Calculate distance moved, bearing, movement types, movement groups ---------- ###

### --- Distance + Bearing --- ###

# add distance/bearing variables for distance moved between previous and current location, and bearing taken to get there
# first location will be NA (ie. no previous location)

# for all rows of each list element of allbird, create new variables: distance, bearing

# reference points 

newallbird <- list()

for (i in 1:length(allbird)){
  
  ### --- STEP 1: split list according to resample no. --- ###
  dat <- allbird[[i]]
  
  # remove reference points and any other weird points outside of the range we'd expect
  
  coordinates(dat) <- c("newlongs","newlats")
  proj4string(dat) <- CRS("+init=epsg:3395") # set projection and datum information for cuckoo data, which are in wGS 84/World Mercator projection coordinates (set earlier in script)
  newdat <- spTransform(dat, CRS("+proj=longlat +datum=WGS84")) # transform to long/lats which are required for the distance calculation
  
  newdat <- as.data.frame(newdat)
  newdat <- Subset(newdat, newlongs > -19) # remove reference points in North America
  
  l <- droplevels(newdat$tcycle)
  resample.no <- c()
  
  for (j in 1:100){
    x <- rep(j,length(levels(l)))
    resample.no <- c(resample.no, x)
  }
  
  newdat <- data.frame(newdat, resample.no)
  datbysample <- split(newdat, resample.no)
  
  ### --- STEP 2: lapply function distbearmvmt across all 100 resamples in list datbysample, to calculate distance moved, bearing, and movement groups and types between tcycles --- ###
  
  ### --- FUNCTION to measure distance moved between each transmission cycle, bearing direction between each tcycle, and whether movements between tcycles were stopovers or migratory --- ###
  
  distbearmvmt <- function(input){ # requires 2 variables called newlongs & newlats
    
    ### DISTANCE & BEARING
    dist.cuckoo <- input[,c("newlongs","newlats")]
    coordinates(dist.cuckoo) <- c("newlongs", "newlats")
    proj4string(dist.cuckoo) <- CRS("+proj=longlat +datum=WGS84") # set coordinates as long/lats which are required for the distance calculation
    
    d1 <- dist.cuckoo@coords
    d2 <- dist.cuckoo@coords[-1,] # create new matrix minus the first row so that d2 starts at the second observation
    last.d2 <- matrix(c(0,0), nrow=1, dimnames=list(1, c("newlongs","newlats"))) # create a placeholder last row in the second distance matrix
    d2 <- rbind(d2, c(0,0))
    
    dist <- distCosine(d1,d2)/1000 # distance between points, in km
    bear <- bearing(d1,d2) # bearing between points
    dist <- c(NA,dist[-length(dist)])
    bear <- c(NA,bear[-length(bear)])
    
    distbear <- data.frame(input, distance=dist, bearing=bear)
    
    ### MOVEMENT GROUPS
    mgroup <- rep(NA,nrow(distbear))
    
    for (n in 1:nrow(distbear)){
      if (is.na(distbear$distance[n])) {
        mgroup[n] <- 1
      } else if (distbear$distance[n] <= 30) {
        mgroup[n] <- mgroup[n-1]
      } else {mgroup[n] <- mgroup[n-1] + 1}
    }
    
    ### MOVEMENT TYPES
    mtype <- c("C", rep(NA, nrow(distbear)-1))
    
    for (n in 2:(nrow(distbear)-1)){
      if (mgroup[n] != mgroup[n-1] & mgroup[n] != mgroup[n+1]) {
        mtype[n] <- "M"
      } else {
        mtype[n] <- "S"
      }
    }
    
    completedata <- data.frame(distbear,mgroup,mtype)
    return(completedata)
    
  }
  
  dat.allvar <- lapply(datbysample, distbearmvmt)
  dat.allvar.onebird <- do.call(rbind, dat.allvar)
  newallbird[[i]] <- dat.allvar.onebird
  
}


# write data files with new variables (distance, movement groups, etc) to .csv - one large with all birds, and separate files for each bird

new.allbirdfinal <- do.call(rbind, newallbird) # rbind all list elements

write.csv(new.allbirdfinal, "allbirdfinal - resampling with distances and movement groups.csv", row.names=FALSE) # write entire file to csv

# write each individual's bootstrapped data to an individual csv to be called again later

setwd("C:/Users/samf/Documents/R/projects/cuckoos/data/resampled data with uncertainty 100 bootstraps + distance, movement groups, etc")

for (i in 1:length(newallbird)){
  write.csv(newallbird[[i]], paste(newallbird[[i]]$name[1],".csv", sep=""), row.names=FALSE)
}
    
######################################################
######################################################
######################################################
######################################################
######################################################

  
# ### --- Movement group + Movement type --- ###
# 
# # for each movement, classify movement type as migratory or stopover
# # if distance < 25km, movement is a stopover (S), if >25km, movement is migration (M)
# # ?????????? work out what the best distance is to use for categorising movements - somewhere between 20-30km?
# # for every group of movements (ie. SSSSMMMSSS are 3 groups), give those observations the same movement group number
# 
# ### Histogram of cuckoo distance movements
# 
# hist(cuckoo.use$distance)
# 
# par(mfrow=c(2,2))
# hist(subset(cuckoo.use, distance < 100)$distance, main=c("d < 100km"), xlab=c("Distance (km)"))
# hist(subset(cuckoo.use, distance > 100 & distance < 500)$distance, main=c("100km < d < 500km"), xlab=c("Distance (km)"))
# hist(subset(cuckoo.use, distance > 500)$distance, main=c("d > 500km"), xlab=c("Distance (km)"))
# hist(subset(cuckoo.use, distance > 20 & distance < 200)$distance, main=c("20km < d < 200km"), breaks=20, xlab=c("Distance (km)"))
# 
# 
# # use <= 30km as the threshold for distinguishing between a stopover (S) and a migration (M), and between movement groups (e.g. 2 stopovers that might not be separated by a migration)
# 
# ### Movement numbers
# 
# # make a blank variable in which to put movement groups
# mgroup <- rep(NA,nrow(cuckoo.use))
# 
# for (n in 1:nrow(cuckoo.use)){
#   if (is.na(cuckoo.use$distance[n])) {
#     # or alternatively, for loop counts from 2:nrow... and if(is.na(cuckoo.use$distance[n]))
#     mgroup[n] <- 1
#   } else if (cuckoo.use$distance[n] <= 30) {
#     mgroup[n] <- mgroup[n-1]
#   } else {mgroup[n] <- mgroup[n-1] + 1}
# }
# 
# ###################################################
# 
# ##################### OLD CODE FOR MOVEMENT TYPE ##########################
# # create new variable for movement type
# # if distance is NA, then this is the capture occasion ("C")
# # if bearing is NA, then could be capture occasion or bird didn't move between subsequent observations
# 
# # mtype <- ifelse(cuckoo.use$distance <= 30, c("S"), c("M"))
# # mtype[is.na(mtype)] <- c("C")
# # mtype <- as.factor(mtype)
# 
# ####################################################
# 
# ######### OLD CODE FOR MOVEMENT GROUP USING "M" and "S" QUALIFIERS, WHICH WILL GIVE INCORRECT RESULTS IF TWO STOPOVERS IN SEQUENCE ARE IN DIFFERENT LOCATIONS #########
# # for each individual cuckoo (call using $name), get row number [x] of first observation for that bird (the capture observation) - that is starting point, and mgroup[x] <- 1, then proceed into for loop, where i counts (x+1):nrow(cuckoo$name[individual name])
# # 
# # for (j in 1:31){ # loop counts through cuckoo names
# #   a <- min(which(cuckoo.use$name==levels(cuckoo.use$name)[j])) # row number of first observation of cuckoo j
# #   b <- max(which(cuckoo.use$name==levels(cuckoo.use$name)[j])) # row number of last observation of cuckoo j
# # 
# #   mgroup[a] <- 1  # set capture occations as movement group 1
# #   
# #   for(i in (a+1):b){  # loops through rows of each cuckoo
# #     if(mtype[i] == c("S") & mtype[i-1] == c("S")){          # stopover
# #       mgroup[i] <- mgroup[i-1]
# #     } else if(mtype[i] == c("S") & mtype[i-1] == c("M")){   # stopover
# #       mgroup[i] <- mgroup[i-1]
# #     } else if(mtype[i] == c("M") & mtype[i-1] == c("S")){   # migration
# #       mgroup[i] <- mgroup[i-1] + 1
# #     } else if(mtype[i] == c("M") & mtype[i-1] == c("M")){   # migration
# #       mgroup[i] <- mgroup[i-1] + 1
# #     } else if(mtype[i] == c("M") & mtype[i-1] == c("C")){   # if migration directly follows the capture occasion
# #       mgroup[i] <- mgroup[i-1] + 1
# #     } else {
# #       mgroup[i] <- mgroup[i-1]                          # stopover, if capture is followed by "S"
# #     }
# #   }
# #   
# # }
# 
# ####################################################
# 
# # add a variable, mtype2, to be able to subset data to include only movement groups where the # of rows for that group number are > 2
# 
# # for each cuckoo, if mgroup value is the same as the one before it and after it, then mtype2 = S, otherwise mtype2 = M, EXCEPT for the first observation for every cuckoo, which was the capture occasion
# # aka, if a bird has > 2 observations in a single movement group, then it constitutes as a stopover, even if the first observation was categorised as a migration (which brought it to that stopover)
# # in order to change the number of consecutive locations that is required to equal a stopover from 2 to something else, then will need to modify if/else statement
# 
# firsts <- rep(NA,31) # blank variable to hold row numbers of first row of data for each individual
# lasts <- rep(NA,31) # blank variable to hold row numbers of last row of data for each individual
# 
# for (j in 1:31){ # loop counts through cuckoo names
#   firsts[j] <- min(which(cuckoo.use$name==levels(cuckoo.use$name)[j])) # row number of first observation of cuckoo j
#   lasts[j] <- max(which(cuckoo.use$name==levels(cuckoo.use$name)[j])) # row number of last observation of cuckoo j
# }
# 
# # loop to create mtype2 based on mgroup variable
# mtype <- rep(NA, nrow(cuckoo.use))
# 
# for (n in 1:nrow(cuckoo.use)){
#   
#   if (n %in% firsts == TRUE) {
#     mtype[n] <- "C"
#   } else {
#     if (mgroup[n] != mgroup[n-1] & mgroup[n] != mgroup[n+1]) {
#       mtype[n] <- "M"
#     } else {
#       mtype[n] <- "S"
#     }
#   }
# }
# 
# newcuckoo <- data.frame(cuckoo.use, mgroup, mtype)
# newcuckoo[1:999,c("name","mgroup","mtype")]     # check output
# 
# write.csv(newcuckoo, "cuckoo movements for analysis.csv", row.names=FALSE)



