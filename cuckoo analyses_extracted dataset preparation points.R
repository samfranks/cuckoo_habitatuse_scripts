##########################################################
#
#   CUCKOO ANALYSES: DATASET PREPARATION
#
#   Code to get extracted environmental data into format suitable for analysis: observation is an mgroup, associated variables are name, laststop, sahara.success, year, PAoverlap, mean proportion of each habitat within stopover, mean elevation, mean climate
#   Data are written to a .csv to be used for analyses
#
#  Samantha Franks
#  24 May 2014
#  11 July 2014: prepares new extracted data, including SPEI and multiple random absences
#  Aug 2014: using points instead of polygon extracted data
#
##########################################################

library(plyr)
library(reshape)
library(dummies)

#analysisdate <- 20140516

###----------------------------------------------###
####         SET WORKING DIRECTORIES		####
###----------------------------------------------###

Mac <- FALSE

if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

extractionwd <- c("/corine PA elevation spei extracted values points/")

#### SET OTHER VARIABLES ####

#randomradius <- c(50, 100, 200, 500)
randomradius <- c(500)

present <- FALSE

#### LOOP THROUGH RANDOM RADIUS OPTIONS ####

if (!present) {
  numberrepeats <- length(randomradius)
}

if (present) {
  numberrepeats <- 1
}

for (z in 1:numberrepeats) {
  
  ###------------------------------------------------------###
  ####         SET-UP CORINE DATA & CONVERT ITS FORMAT		####
  ###------------------------------------------------------###
  
  
  ############## FOR POINTS ###############
  # for loop for each cuckoo file starts here
  corine.data <- list()
  climate.data <- list()
  elevation.data <- list()
  PAoverlap.data <- list()
  presentdata <- list()
  absentdata.byscale <- list()
  
  for (a in 1:28) { # START LOOP
    
    ###--- LOAD CORINE DATASETS ---###
    if (present) {
      setwd(paste(datawd, extractionwd, "/corine/presences", sep=""))
    }
    
    if (!present) {
      setwd(paste(datawd, extractionwd, "/corine/", randomradius[z], " km radius", sep=""))
    }
    
    corine <- read.csv(list.files()[a], header=TRUE)
    
    ### import corine legend
    setwd(datawd)
    clclegend <- read.csv("clc_legend.csv", header=T)
    
    ###--- Add names of corine land types to dataset---###
    
    # function to merge info in clclegend dataframe with dataset
    
    extractcorinedata <- function(datasetname){
      landclass <- merge(x=datasetname, y=clclegend[,c("GRID_CODE", "LAND.CLASS", "LABEL4")], by.x=c("corine.values"), by.y=c("GRID_CODE"))
      #       landclass.order <- landclass[order(landclass$mgroup, landclass$corine.values),]
      #       if (present) {
      #         landclass.order2 <- with(landclass.order, data.frame(name, mgroup, year, laststop, strategy, Sahara.success, LAND.CLASS, LABEL4))
      #       }
      #       if (!present) {
      #         landclass.order2 <- with(landclass.order, data.frame(name, mgroup, year, laststop, strategy, Sahara.success, nabsence, LAND.CLASS, LABEL4))
      #       }
      return(landclass)    
    }
    
    newcorine <- extractcorinedata(corine)
    x <- dummy(newcorine$LAND.CLASS, sep="_", drop=FALSE)
    corine.data[[a]] <- data.frame(newcorine, x)
    
    ###----------------------------------###
    ####         ADD CLIMATE DATA		####
    ###----------------------------------###
    
    
    ###--- LOAD DATA ACTUAL & RANDOM MCP CORINE DATASETS ---###
    
    if (present) {
      setwd(paste(datawd, extractionwd, "/climate/presences", sep=""))
    }
    
    if (!present) {
      setwd(paste(datawd, extractionwd, "/climate/", randomradius[z], " km radius", sep=""))
    }
    
    climate.data[[a]] <- read.csv(list.files()[a], header=TRUE)
    
    
    ###----------------------------------------###
    ####         ADD ELEVATION DATA		####
    ###----------------------------------------###
    
    if (present) {
      setwd(paste(datawd, extractionwd, "/elevation/presences", sep=""))
    }
    
    if (!present) {
      setwd(paste(datawd, extractionwd, "/elevation/", randomradius[z], " km radius", sep=""))
    }
    
    elevation <- read.csv(list.files()[a], header=TRUE)
    newelevation <- na.omit(elevation)
    newelevation[newelevation$elevation.m <0, "elevation.m"] <- 0
    
    elevation.data[[a]] <- newelevation
    
    ###----------------------------------------------###
    ####         ADD PROTECTED AREA DATA		####
    ###----------------------------------------------###
    
    if (present) {
      setwd(paste(datawd, extractionwd, "/protected areas overlap/presences", sep=""))
    }
    
    if (!present) {
      setwd(paste(datawd, extractionwd, "/protected areas overlap/", randomradius[z], " km radius", sep=""))
    }
    
    PAoverlap.data[[a]] <- read.csv(list.files()[a], header=TRUE)
    
    
    ###--------------------------------###
    ####        MERGE ALL DATA  	    ####
    ###--------------------------------###
    
    
    
    # merge corinecombine, climate, newelevation, PAdesigcombine
    
    if (present) {
      
      mergeby <- names(climate.data[[a]])[-which(names(climate.data[[a]]) %in% c("spei.Mar","spei.Aug"))]
      mergeby.noresampcoor <- names(elevation.data[[a]])[-which(names(elevation.data[[a]]) %in% c("elevation.m","newlongs","newlats"))]
      
      alldata1 <- merge(corine.data[[a]],climate.data[[a]], by.x=mergeby, by.y=mergeby)
      
      alldata2 <- merge(alldata1, PAoverlap.data[[a]], by.x=mergeby, by.y=mergeby)
      
      alldata3 <- merge(alldata2, elevation.data[[a]], by.x=mergeby.noresampcoor, by.y=mergeby.noresampcoor)
      
      alldata4 <- data.frame(alldata3, presence=1)
      
      presentdata[[a]] <- alldata4[order(alldata4$mgroup, alldata4$tcycle),]
      
      presentdata[[a]] <- rename(presentdata[[a]], c("newlongs.x"="newlongs.epsg3035", "newlats.x"="newlats.epsg3035", "newlongs.y"="newlongs.epsg4326", "newlats.y"="newlats.epsg4326"))
      
    }
    
    if (!present) {
      
      mergeby <- names(climate.data[[a]])[-which(names(climate.data[[a]]) %in% c("spei.Mar","spei.Aug"))]
      mergeby.noresampcoor <- names(elevation.data[[a]])[-which(names(elevation.data[[a]]) %in% c("elevation.m","randomlong","randomlat"))]
      
      alldata1 <- merge(corine.data[[a]],climate.data[[a]], by.x=mergeby, by.y=mergeby)
      
      alldata2 <- merge(alldata1, PAoverlap.data[[a]], by.x=mergeby, by.y=mergeby)
      
      alldata3 <- merge(alldata2, elevation.data[[a]], by.x=mergeby.noresampcoor, by.y=mergeby.noresampcoor)
      
      alldata4 <- data.frame(alldata3, presence=0, random.scale=randomradius[z])
      
      absentdata.byscale[[a]] <- alldata4[order(alldata4$mgroup, alldata4$tcycle, alldata4$nabsence),]
      
      absentdata.byscale[[a]] <- rename(absentdata.byscale[[a]], c("randomlong.x"="randomlong.epsg3035", "randomlat.x"="randomlat.epsg3035", "randomlong.y"="randomlong.epsg4326", "randomlat.y"="randomlat.epsg4326"))
      
    }
    
    
  } # END loop through individuals
  
}
  
  
  ###--------------------------------###
  ####        WRITE DATA       ####
  ###--------------------------------###
  
  
  if (present) {
    all.presentdata <- do.call(rbind, presentdata)
    
    setwd(paste(datawd, "/data for analysis/", sep=""))
    write.csv(all.presentdata, "presence data all variables points.csv", row.names=FALSE)
  }
  
  if (!present) {
    all.absentdata <- do.call(rbind, absentdata.byscale)
    
    setwd(paste(datawd, "/data for analysis/", sep=""))
    write.csv(all.absentdata, paste("absence data all variables points ", randomradius[z], " km.csv", sep=""), row.names=FALSE)
  }