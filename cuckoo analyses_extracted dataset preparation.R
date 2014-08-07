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
#
##########################################################

library(plyr)
library(reshape)

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

randomradius <- c(50, 100, 200, 500)

absentdata <- list()

present <- TRUE

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
  
  ############## FOR POINTS ################
  
  ###--- LOAD CORINE DATASETS ---###
  if (present) {
    setwd(paste(datawd, extractionwd, "/corine", sep=""))
    corine <- read.table("corine all birds.csv", sep=",")
    corine <- read.csv("corine all birds.csv", header=TRUE, quote = "")
    corine2 <- read.csv("corine all birds.csv", header=TRUE)
    
      
  }
  
  if (!present) {
    #setwd(paste(datawd, extractionwd, "/corine random/", randomradius[z], " km radius", sep=""))
  }
  
  
  ### import corine legend
  setwd(datawd)
  clclegend <- read.csv("clc_legend.csv", header=T)
  
  
  
  ############## FOR POLYGONS ###############
  # for loop for each cuckoo file starts here
  corinepropdata <- list()
  
  for (a in 1:28) { # START LOOP
    
    ###--- LOAD CORINE DATASETS ---###
    if (present) {
      setwd(paste(datawd, extractionwd, "/corine", sep=""))
    }
    
    if (!present) {
      setwd(paste(datawd, extractionwd, "/corine random/", randomradius[z], " km radius", sep=""))
    }
    
    corine <- read.csv(list.files()[a], header=TRUE)
    
    ### import corine legend
    setwd(datawd)
    clclegend <- read.csv("clc_legend.csv", header=T)
    
    ###--- Add names of corine land types to dataset---###
    
    # function to merge info in clclegend dataframe with dataset
    
    extractcorinedata <- function(datasetname){
      landclass <- merge(x=datasetname, y=clclegend[,c("GRID_CODE", "LAND.CLASS", "LABEL4")], by.x=c("corine.values"), by.y=c("GRID_CODE"))
      landclass.order <- landclass[order(landclass$mgroup, landclass$corine.values),]
      if (present) {
        landclass.order2 <- with(landclass.order, data.frame(name, mgroup, year, laststop, strategy, Sahara.success, LAND.CLASS, LABEL4))
      }
      if (!present) {
        landclass.order2 <- with(landclass.order, data.frame(name, mgroup, year, laststop, strategy, Sahara.success, nabsence, LAND.CLASS, LABEL4))
      }
      return(landclass.order2)    
    }
    
    newcorine <- extractcorinedata(corine)
    
    if (present) {
      x <- as.data.frame(prop.table(table(newcorine$LAND.CLASS, newcorine$mgroup),2))
      corine.prop <- reshape(x, timevar="Var1", idvar="Var2", direction="wide")
      colnames(corine.prop) <- c("mgroup", levels(newcorine$LAND.CLASS))
      
      corine.prop2 <- unique(merge(newcorine[,c("name","mgroup","year","laststop","strategy","Sahara.success")], corine.prop))
      corinepropdata[[a]] <- data.frame(name=corine.prop2$name, mgroup=corine.prop2$mgroup, corine.prop2[,c("year","laststop","strategy","Sahara.success","agriculture","forest","scrub.grassland","unsuitable","wetland.water")])
    }
    
    if (!present) {
      split.nabsence <- as.factor(newcorine$nabsence)
      y <- split(newcorine, split.nabsence)
      
      calc.corineprop <- function(corine.nabsence) {
        x <- as.data.frame(prop.table(table(corine.nabsence$LAND.CLASS, corine.nabsence$mgroup), 2))
        corine.prop <- reshape(x, timevar="Var1", idvar="Var2", direction="wide")
        colnames(corine.prop) <- c("mgroup", levels(corine.nabsence$LAND.CLASS))
        corine.prop2 <- unique(merge(corine.nabsence[,c("name","mgroup","year","laststop","strategy","Sahara.success","nabsence")], corine.prop))
        corineoutput <- data.frame(name=corine.prop2$name, mgroup=corine.prop2$mgroup, corine.prop2[,c("year","laststop","strategy","Sahara.success","nabsence", "agriculture","forest","scrub.grassland","unsuitable","wetland.water")])
        return(corineoutput)
      }
      
      corinepropdata.byabsence <- lapply(y, calc.corineprop)
      corinepropdata[[a]] <- do.call(rbind, corinepropdata.byabsence)
      
    }
    
  } # END loop through individuals
  
  corinecombine <- do.call(rbind, corinepropdata)
  
  ###----------------------------------###
  ####         ADD CLIMATE DATA		####
  ###----------------------------------###
  
  
  ###--- LOAD DATA ACTUAL & RANDOM MCP CORINE DATASETS ---###
  setwd(paste(datawd, extractionwd, "/climate", sep=""))
  
  if (present) {
    climate <- read.csv("climate all birds.csv", header=TRUE)
  }
  
  if (!present) {
    climate <- read.csv(paste("climate random ", randomradius[z], " km radius all birds.csv", sep=""), header=TRUE)
  }
  
  ### old code when extracting mean winter and spring precip & temp
  #climatesub <- climate[,c(1:2,4:6)]
  #climatecombine <- reshape(climatesub, timevar="season", idvar=c("name","mgroup"), direction="wide")
  
  ###----------------------------------------###
  ####         ADD ELEVATION DATA		####
  ###----------------------------------------###
  
  setwd(paste(datawd, extractionwd, "/elevation", sep=""))
  
  if (present) {
    elevation <- read.csv("elevation all birds.csv", header=TRUE)
    elevationmeans <- aggregate(elevation$elevation.m, by=list(name=elevation$name, mgroup=elevation$mgroup), mean, na.rm=TRUE)
    elevationcombine <- elevationmeans[order(elevationmeans$name,elevationmeans$mgroup),]
    colnames(elevationcombine) <- c("name","mgroup","elevation.m")
    
    # for all mean elevations below sea level (only 3 of them - Chris6, Chris52, Iolo4), make values a 0
    elevationcombine[elevationcombine$elevation.m < 0, "elevation.m"] <- 0
  }
  
  if (!present) {
    elevation <- read.csv(paste("elevation random ", randomradius[z], " km radius all birds.csv", sep=""), header=TRUE)
    elevationmeans <- aggregate(elevation$elevation.m, by=list(nabsence=elevation$nabsence, name=elevation$name, mgroup=elevation$mgroup), mean, na.rm=TRUE)
    x <- elevationmeans[order(elevationmeans$name,elevationmeans$mgroup, elevationmeans$nabsence),]
    elevationcombine <- data.frame(name=x$name, mgroup=x$mgroup, nabsence=x$nabsence, elevation.m=x$x)
    
    # some polygons have failed to extract any elevation (e.g. Wallace6_9 at 500km) - this removes those polygons
    NArows <- which(is.na(elevationcombine$elevation.m))
    NAids <- data.frame(elevationcombine[NArows,c("name","mgroup","nabsence")])
    NAids
    newelevation <- na.omit(elevationcombine)
    #elevationcombine[elevationcombine$elevation.m < 0, "elevation.m"] <- 0
    newelevation[newelevation$elevation.m < 0, "elevation.m"] <- 0
  }
  
  
  
  #   ### RUN 1 from 10/07/2014: check observations which output NAs for elevation; Iolo points are in IJsselmeer in NLD, will have to rerun absence generation and reduce amount of allowed corine category 41 (water bodies)
    # setwd(paste(outputwd,"/bootstrapped absence polygons/", sep=""))
    # prob50 <- readOGR(dsn=".", "Wallace 500 km random bootstrapped polygons")
    # #Iolo4_4, Iolo4_5, Iolo4_10
  # #Wallace6_9
    
    # #problemstops <- subset(prob50, ID=="Iolo4_4" | ID=="Iolo4_5" | ID=="Iolo4_10")
  # problemstops <- subset(prob50, ID=="Wallace6_9")
  
    
    # sapply(slot(problemstops, "polygons"), function(i) (slot(i,"labpt")))
    
    # cropextent <- extent(3387) # bounding box extent for relevant parts of Europe (no Scandanavia or Turkey)
    # Ioloraster <- crop(r,cropextent)
    
    # prob50.4326 <- spTransform(prob50, epsg4326)
    # plotKML(prob50.4326, paste(levels(prob50.4326$name), "polygons",sep=""), alpha=0.75)
  
  
  
  ###----------------------------------------------###
  ####         ADD PROTECTED AREA DATA		####
  ###----------------------------------------------###
  
  setwd(paste(datawd, extractionwd, "/protected areas", sep=""))
  
  if (present) {
    PAoverlap <- read.csv("PA overlap all birds.csv", header=TRUE)
    PAdetails <- read.csv("PA site details all birds.csv", header=TRUE)
  }
  
  if (!present) {
    PAoverlap <- read.csv(paste("stopover overlap area random ", randomradius[z], " km radius all birds.csv", sep=""), header=TRUE)
    PAdetails <- read.csv(paste("PA site details random ", randomradius[z], " km radius all birds.csv", sep=""), header=TRUE)
  }
  
  
  PAdetails <- rename(PAdetails, c("name"="sitename", "birdname"="name"))
  PAoverlap <- rename(PAoverlap, c("MCPid"="mgroup"))
  
  overlap <- factor(levels=c("Y","N"))
  for (i in 1:nrow(PAoverlap)) {
    if (PAoverlap$prop.MCP.overlap[i] > 0) {
      overlap[i] <- "Y"
    } else {overlap[i] <- "N"}
  }
  PAoverlap2 <- cbind(PAoverlap, overlap)
  
  
  if (present) {
    
    largeststop <- ddply(PAdetails, c("name","mgroup"), summarize, MCP.overlap.area.km2=max(MCP.overlap.area.km2), largestMCP="Y")
    
    PAdetails2 <- merge(largeststop, PAdetails, all.y=TRUE)
    
    PAdetails2$largestMCP[which(is.na(PAdetails2$largestMCP))] <- "N"
    
    islargeststop <- subset(PAdetails2, largestMCP=="Y")
    islargeststop <- droplevels(islargeststop)
    
    #   # Iolo, Skinner, and Whortle each have 1 stopover with the largest overlap being with 2 PAs (an SPA and an SAC) - given they have the same designation, will just use one of them
    #   islargeststop[islargeststop$name=="Iolo",]
    #   islargeststop[islargeststop$name=="Skinner",]
    #   islargeststop[islargeststop$name=="Whortle",]
    
    dupmgroups <- which(duplicated(islargeststop[,c("name", "mgroup")]))
    
    islargeststop2 <- islargeststop[-dupmgroups,]
    
    mergePAdata <- merge(PAoverlap2, islargeststop2, by.x=c("name","mgroup"), by.y=c("name","mgroup"), all.x=TRUE)
    
    PAdatacombine <- with(mergePAdata, data.frame(name, mgroup, stoparea.km2=MCP.area.km2, overlap, overlaparea.km2=PA.overlap.area.km2, prop.overlap=prop.MCP.overlap, desig_type))
  }
  
  if (!present) {
    
    split.nabsence <- as.factor(PAdetails$nabsence)
    PAdetailssplit <- split(PAdetails, split.nabsence)
    split.nabsence <- as.factor(PAoverlap2$nabsence)
    PAoverlapsplit <- split(PAoverlap2, split.nabsence)
    
    PAdesigcombine <- list()
    
    for (i in 1:10) {
      
      dd.details <- PAdetailssplit[[i]]
      dd.overlap <- PAoverlapsplit[[i]]
      
      largeststop <- ddply(dd.details, c("name","mgroup"), summarize, MCP.overlap.area.km2=max(MCP.overlap.area.km2), largestMCP="Y")
      PAdetails2 <- merge(largeststop, dd.details, all.y=TRUE)
      PAdetails2$largestMCP[which(is.na(PAdetails2$largestMCP))] <- "N"
      islargeststop <- subset(PAdetails2, largestMCP=="Y")
      islargeststop <- droplevels(islargeststop)
      
      dupmgroups <- which(duplicated(islargeststop[,c("name", "mgroup")]))
      
      islargeststop2 <- islargeststop[-dupmgroups,]
      
      mergePAdata <- merge(dd.overlap, islargeststop2, by.x=c("name","mgroup","nabsence"), by.y=c("name","mgroup","nabsence"), all.x=TRUE)
      
      PAdesigcombine[[i]] <- with(mergePAdata, data.frame(name, mgroup, nabsence, stoparea.km2=MCP.area.km2, overlap, overlaparea.km2=PA.overlap.area.km2, prop.overlap=prop.MCP.overlap, desig_type)) # variables omitted: largestMCP, wdpaid, country, sitename, desig_eng, iucn_cat
      
    }
    
    PAdatacombine <- do.call(rbind, PAdesigcombine)
    
  }
  
  ###--------------------------------###
  ####        MERGE ALL DATA  	    ####
  ###--------------------------------###
  
  # merge corinecombine, climate, newelevation, PAdesigcombine
  
  if (present) {
  	alldata1 <- merge(corinecombine,climate, by.x=c("name","mgroup","year","laststop","strategy","Sahara.success"), by.y=c("name","mgroup","year","laststop","strategy","Sahara.success"))
  	
  	alldata2 <- merge(alldata1, elevationcombine, by.x=c("name","mgroup"), by.y=c("name","mgroup"))
  	
  	alldata3 <- merge(alldata2, PAdatacombine, by.x=c("name","mgroup"), by.y=c("name","mgroup"))
  	
  	alldata4 <- data.frame(alldata3, presence=1)

    presentdata <- alldata4[order(alldata4$name,alldata4$mgroup),]
    presentdata <- rename(presentdata, c("year.x"="year", "laststop.x"="laststop", "strategy.x"="strategy","Sahara.success.x"="Sahara.success"))
	}
	
if (!present) {
	alldata1 <- merge(corinecombine,climate, by.x=c("name","mgroup","year","laststop","strategy","Sahara.success","nabsence"), by.y=c("name","mgroup","year","laststop","strategy","Sahara.success","nabsence"))
	
	alldata2 <- merge(alldata1, newelevation, by.x=c("name","mgroup","nabsence"), by.y=c("name","mgroup","nabsence"))
	
	alldata3 <- merge(alldata2, PAdatacombine, by.x=c("name","mgroup","nabsence"), by.y=c("name","mgroup","nabsence"))
  
    alldata4 <- data.frame(alldata3, presence=0)
      
    absentdatatemp <- alldata4[order(alldata4$name,alldata4$mgroup,alldata4$nabsence),]
    absentdata[[z]] <- data.frame(absentdatatemp, random.scale=randomradius[z])

	}
  
}

###--------------------------------###
####        WRITE DATA       ####
###--------------------------------###

if (present) {
  setwd(paste(datawd, "/data for analysis/", sep=""))
  write.csv(presentdata, "presence data all variables.csv", row.names=FALSE)
}

if (!present) {
  allabsentdata <- do.call(rbind, absentdata)
  setwd(paste(datawd, "/data for analysis/", sep=""))
  write.csv(allabsentdata, "absence data all variables.csv", row.names=FALSE)
}