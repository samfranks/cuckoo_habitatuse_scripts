##########################################################
#
#	Cuckoo African stopovers
#
# Code extracts country & continent geographic info for cuckoo resampled location points, clips dataset to only points occurring in Africa, creates a MCP around each stopover
#
# Combines all stopover polygons for cuckoos with Africa locations and exports to a shapefile (also can export separate shapefiles for each individual cuckoo)
#
#	Samantha Franks
# 10 Dec 2013
# 23 Dec 2013 - rerun with re-edited resampled dat
# 2 Sep 2014 - for Phil for Leventis visit
#
#
##########################################################

rm(list=ls())

###--- LOAD PACKAGES ---###

library(adehabitat) # 
library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
#library(shapefiles) # allows you to read and write ESRI shapefiles
library(rgeos)
library(raster)
library(geosphere)
library(adehabitatHR)
library(rworldmap) # getMap function
library(rworldxtra)
library(chron)
library(plyr)

##################################################
##################################################

###-----------------------------------------------------------###
#     PARENT INFORMATION ON VARIABLES (passed from PARENT SCRIPT)
###-----------------------------------------------------------###

### NOTE: if using as source code, then can leave out LOAD DATASET section and use the data passed from parent script.  Loop should be set up in the parent script to run through all datasets in a folder

original <- TRUE # whether original or resampled data
# dataset <- x, where x is the name of the dataset you want to extract info from

# Clip dataset to only southward migration breeding May 10 to Dec 21 (julian >=131 and <=365)
jdate1 <- 131 #minimum julian date bracketing period of interest
jdate2 <- 365 #maximum julian date bracketing period of interest

continentname <- "Africa" #name of the continent you want to extract info from

# countryname <- y, where y is the name of the country/countries you want to extract info from

select.mtype <- "S"

##################################################
##################################################

###-----------------------------------------------------------###
#         LOAD MAP BASE LAYERS
###-----------------------------------------------------------###

### ------- Europe/Africa national boundaries layer ------- ###

dat.world <- getMap(resolution="high")

Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]

Eur.Afr <- spTransform(Eur.Afr,CRS("+init=epsg:4326")) # transform map to exactly the same projection as cuckoo newlongs/newlats data

Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]


###-----------------------------------------------------------###
#         SET WORKING DIRECTORIES
###-----------------------------------------------------------###

parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")

datawd <- ("C:/Users/samf/Documents/Git/cuckoos/data")
outputwd <- ("C:/Users/samf/Documents/Git/cuckoos/output")

resampwd <- ("/resampled data with uncertainty 100 bootstraps")
origwd <- ("/Movebank original data new birds + distance, movement groups, etc")

orig.data.outputwd <- ("/original data + extra grouping variables Africa")
resamp.data.outputwd <- ("/resampled data + extra grouping variables")

##################################################
##################################################

allbirddata <- list()

for (a in 1:45) {
  
  ###-----------------------------------------------------------###
  #         LOAD DATASET - left out if used as source code
  ###-----------------------------------------------------------###
  
  # change to resampwd or origwd depending on whether corine values are being extracted for resampled or original data
  
  if (original) setwd(paste(datawd,origwd,sep=""))
  if (!original) setwd(paste(datawd,resampwd,sep=""))
  
  cuckoofiles <- list.files()
  
  # a <- 3 #test
  
  # if this is being used as source code, then "dataset" should be the name of the dataset passed to this code
  dataset <- read.csv(cuckoofiles[a], header=T)
  dataset <- dataset[-which(names(dataset) == c("datetime"))]
  dataset$mgroup <- as.factor(dataset$mgroup)
  
  #if (levels(dataset$name) == "Karma" | levels(dataset$name) == "Wistman") next
  
  ##################################################
  ##################################################
  
  
  ###-----------------------------------------------------------###
  #         EXTRACT GEOINFO AND ADD TO DATASET
  ###-----------------------------------------------------------###
  
  # long/lat (original) are in CRS("+proj=longlat +datum=WGS84"), or CRS("+init=epsg:4326")
  # long/lat (resampledl) are in CRS("+init=epsg:3395")
  # newlongs/newlats are in CRS("+init=epsg:3395")
  # centroidlong/lat are in CRS("+init=epsg:3395")
  
  if (original) {
    coordinates(dataset) <- c("long","lat")
    proj4string(dataset) <- CRS("+init=epsg:4326")
  }
  
  if (!original) {
    coordinates(dataset) <- c("newlongs","newlats")
    proj4string(dataset) <- CRS("+init=epsg:3395")
    dataset <- spTransform(dataset,CRS("+init=epsg:4326")) # transform newlongs/newlats from World Mercator/WGS84 to longlat WGS84
  }
  
  ### over() retrieves values of Eur.Afr at points described by the dataset
  # extracts continent and administrative country names from Eur.Afr
  geoinfo <- over(dataset, Eur.Afr)[,c("REGION","ADMIN")]
  colnames(geoinfo) <- c("continent", "country")
  geoinfo <- droplevels(geoinfo)
  
  ### add a new factor for continent/country NAs, which should be points at sea
  atsea.pts <- which(is.na(geoinfo$continent))
  
  continent.factors <- c(levels(geoinfo$continent), "at sea")
  levels(geoinfo$continent) <- continent.factors
  geoinfo$continent[atsea.pts] <- "at sea"
  
  country.factors <- c(levels(geoinfo$country),"at sea")
  levels(geoinfo$country) <- country.factors
  geoinfo$country[atsea.pts] <- "at sea"
  
  ### add geographic info to dataset and redefine projection
  
  fulldataset <- data.frame(dataset, geoinfo)
  coordinates(fulldataset) <- c("long","lat")
  proj4string(fulldataset) <- CRS("+init=epsg:4326")
  
  
  ### create continent- or country-specific geoinfo datasets
  
  # if Europe, then also include at sea points
  
  if (continentname == "Europe") {
    geodataset <- subset(fulldataset, continent==continentname | continent=="at sea")
    geodataset@data <- droplevels(geodataset@data)
  }
  
  if (continentname != "Europe") {
    geodataset <- subset(fulldataset, continent==continentname)
    geodataset@data <- droplevels(geodataset@data)
  }
  
  if ("countryname" %in% ls()) {
    geodataset <- subset(fulldataset, country==countryname)
  }
  
  ### break out of current bird to next if no points in dataset of interest
  if (nrow(geodataset)==0) next
  
  ###-----------------------------------------------------------###
  #         SUBSET DATA OF INTEREST (DATES, MTYPES, STOPOVERS WITH > 2 TRANSMISSION CYCLES)
  ###-----------------------------------------------------------###
  
  #   # SUBSET 1 - SEASON (DATES)
  #   subset1 <- subset(geodataset, julian >= jdate1 & julian <= jdate2)
  #   subset1@data <- droplevels(subset1@data)
  
  subset1 <- geodataset
  
  # SUBSET 2 - MOVEMENT TYPE (STOPOVERS ONLY)
  subset2 <- subset(subset1, mtype==select.mtype)
  subset2@data <- droplevels(subset2@data)
  
  #   # SUBSET 3 - MGROUPS WITH # of TCYCLES < 2    ### no tcycles for new 2014 data
  #   subset3 <- subset2
  #   subset3$mgroup <- as.factor(subset3$mgroup)
  #   
  #   for (i in 1:length(levels(as.factor(subset2$mgroup)))) {
  #     
  #     mgroup.number <- as.integer(levels(as.factor(subset2$mgroup))[i])
  #     
  #     temp <- subset(subset2, mgroup==mgroup.number)
  #     temp@data <- droplevels(temp@data)
  #     
  #     if (length(levels(as.factor(temp$tcycle))) < 2) {
  #       subset3 <- subset(subset3, mgroup!=mgroup.number)
  #       subset3@data <- droplevels(subset3@data)
  #     } else {
  #       subset3 <- subset3
  #     }
  #   }
  
  # SUBSET 3 - MGROUPS WITH # of LOCATIONS < 2  ### no tcycles for new 2014 Movebank data, this is similar
  
  # identify which mgroups have only 1 observation
  
  mgroups.keep <- names(which(summary(subset2$mgroup) != 1))
  subset3 <- subset2[which(subset2$mgroup %in% mgroups.keep),]
  
  ### ADD MORE SUBSETTING VARIABLES IF DESIRED
  
  nextdat <- subset3
  
  ###-----------------------------------------------------------###
  #         ADD MIGRATION STRATEGY
  ###-----------------------------------------------------------###
  
  setwd(datawd)
  othervar <- read.csv("cuckoo migratory strategy and Sahara crossing success 2014.csv", header=TRUE)
  
  strategy <- rep(othervar[othervar$name==levels(nextdat$name),"migratory.strategy"], nrow(nextdat))
  
  newdataset <- nextdat
  
  newdataset@data <- data.frame(newdataset@data,strategy)
  
  
  ###-----------------------------------------------------------###
  #         ADD STOPOVER DURATION
  ###-----------------------------------------------------------###
  
  newpresent <- newdataset
  x <- do.call(rbind,(strsplit(as.character(newpresent$timestamp), " ")))
  
  dat.datetime <- chron(x[,1], x[,2], format=c(dates="d/m/y", times="h:m:s"))
  
  present.date <- newpresent
  
  present.date@data <- data.frame(present.date@data, datetime=dat.datetime)
  present.date <- present.date[order(present.date$datetime),]
  
  calculate.LOS <- function(fulldata) {
    x <- tapply(fulldata$datetime, fulldata$mgroup, diff)
    LOS <- do.call(rbind, lapply(x, sum))
    LOStable <- data.frame(mgroup=levels(as.factor(fulldata$mgroup)), LOS)
    newdat.LOS <- merge(fulldata, LOStable, by.x="mgroup", by.y="mgroup")
    newdat.LOS <- newdat.LOS[order(newdat.LOS$datetime),]
    return(newdat.LOS)
  }
  
  makeLOStable <- function(fulldata) {
    x <- tapply(fulldata$datetime, fulldata$mgroup, diff)
    LOS <- do.call(rbind, lapply(x, sum))
    LOStable <- data.frame(name=fulldata$name[1], mgroup=levels(as.factor(fulldata$mgroup)), LOS)
    return(LOStable)
  }
  
  present.LOS <- present.date
  
  present.LOS@data <- calculate.LOS(present.LOS@data)
  
  
  ###-----------------------------------------------------------###
  #         ADD SEASON - autumn/winter, spring - southernmost point
  ###-----------------------------------------------------------###
  
  dat.coor <- data.frame(present.LOS@data, coordinates(present.LOS))
  
  ### seasons in Africa span across the year: need to add a season variable to indicate a new winter time period ###
  
  # calculate time elapsed between each subsequent location and count a new period when elapsed time is > 2 months (60 days - breeding season time in Europe, approximately)
  y <- unclass(round(diff(dat.coor$datetime), 2))
  new_per <- y > 60
  
  # creates a dummy variable for new periods, cumsum increases by 1 when a time gap > 60 days is found (handy trick from Ben Bolker on a Stackoverflow post)
  Africaperiod <- 1+c(0,cumsum(new_per))
  
  # adds the period to the datarame
  datwithperiod <- data.frame(dat.coor, Africaperiod)
  
  # identifies the southernmost point in each period in Africa
  southernmost <- datwithperiod$lat %in% tapply(datwithperiod$lat, datwithperiod$Africaperiod, min)
  
  # creates a dummy variable separating periods separated by the southernmost point: periods before (by default, odd periods since period 1 will always be an autumn/winter period) are the autumn/winter period, periods after (by default, even periods) are the spring periods
  seasonperiod <- 1 + cumsum(southernmost)
  
  if (length(levels(as.factor(seasonperiod))) == 1) {
    datwithseason <- data.frame(datwithperiod, season="autumn.winter")
  } else {
    season <- (Africaperiod + seasonperiod) %% 2 == 0
    season <- factor(season, labels=c("spring", "autumn.winter"))
    
    datwithseason <- data.frame(datwithperiod, season)
  }
  

  ###-----------------------------------------------------------###
  #         OUTPUT NEW DATASETS
  ###-----------------------------------------------------------###
  
  ### outputs are without projection information
  # original data (longs, lats) are in epsg: 4326
  # resampled data (newlongs/newlats) are in epsg: 4326
  # resampled data (longs/lats) are in epsg: 3395
  
  
  
  if (original) setwd(paste(datawd,orig.data.outputwd,sep=""))
  if (!original) setwd(paste(datawd,resamp.data.outputwd,sep=""))
  
  outputfile <- file(paste(levels(datwithseason$name), ".csv", sep=""), "w")
  write.table(datwithseason, file=outputfile, row.names=FALSE, sep=",")
  
  close(outputfile)
  
  allbirddata[[a]] <- datwithseason
  
}

finalallbirddata <- do.call(rbind, allbirddata)


############################################################
############################################################

###------------------------------------------------###
####     CREATING SHAPEFILE OF AFRICA DATAPOINTS  ####
###------------------------------------------------###

allbirdsAfrica <- finalallbirddata

allbirdsAfrica$season <- revalue(finalallbirddata$season, c("autumn.winter"="autumnwinter"))


epsg4326 <- CRS("+init=epsg:4326")

coordinates(allbirdsAfrica) <- c("long","lat")
proj4string(allbirdsAfrica) <- epsg4326

allbirdsAfrica$datetime <- as.character(allbirdsAfrica$datetime)

setwd(paste(parentwd,"/output","/stopover points and polygons/", sep="")) 

writeOGR(allbirdsAfrica, dsn=".", layer="all_stopovers_Africa_2014birds_epsg4326", driver="ESRI Shapefile")


############################################################
############################################################
############################################################
############################################################

# ###-----------------------------------------------------------###
# #         Load base layers
# ###-----------------------------------------------------------###
# 
# ### ------- Europe/Africa national boundaries layer ------- ###
# 
# dat.world <- getMap(resolution="high")
# 
# Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
# Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
# Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]
# 
# Eur.Afr <- spTransform(Eur.Afr,CRS("+proj=longlat +datum=WGS84")) # transform map to exactly the same projection as cuckoo newlongs/newlats data
# 
# Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
# Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]
# 
# # clipAfr <- drawExtent()
# # Afrmap <- crop(Afr.only,clipAfr)
# 
# #writeOGR(Afrmap, dsn = ".", layer = "Africa map", driver = "ESRI Shapefile")
# 
# ###------------------------------------------------------------###
# #         Clip cuckoo stopovers to Europe map
# ###------------------------------------------------------------###
# 
# # create list to hold centroids of stopover points for all cuckoos
# 
# AfrMCPstopovers <- list()
# Afrdata <- list()
# 
# # for loop for each cuckoo file starts here
# for (a in 1:31) {
#   
#   ### --- RESAMPLED cuckoo data: Load and add geoinfo --- ###
#   
#   setwd("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/data/resampled data with uncertainty 100 bootstraps")
#   cuckoofiles <- list.files()
#   
#   cuckoo <- read.csv(cuckoofiles[a], header=T)
#   
#   # long/lat are in CRS("+proj=longlat +datum=WGS84")
#   # newlongs/newlats are in CRS("+init=epsg:3395")
#   # centroidlong/lat are in CRS("+init=epsg:3395")
#   
#   coordinates(cuckoo) <- c("newlongs","newlats")
#   proj4string(cuckoo) <- CRS("+init=epsg:3395")
#   cuckoo <- spTransform(cuckoo,CRS("+proj=longlat +datum=WGS84")) # transform newlongs/newlats from World Mercator/WGS84 to longlat WGS84
#   
#   ### over() retrieves values of Eur.Afr at points described by cuckoo
#   # extracts continent and administrative country names from Eur.Afr
#   geoinfo <- over(cuckoo, Eur.Afr)[,c("REGION","ADMIN")]
#   colnames(geoinfo) <- c("continent", "country")
#   geoinfo <- droplevels(geoinfo)
#   
#   ### add a new factor for continent/country NAs, which should be points at sea
#   atsea.pts <- which(is.na(geoinfo$continent))
#   
#   continent.factors <- c(levels(geoinfo$continent), "at sea")
#   levels(geoinfo$continent) <- continent.factors
#   geoinfo$continent[atsea.pts] <- "at sea"
#   
#   country.factors <- c(levels(geoinfo$country),"at sea")
#   levels(geoinfo$country) <- country.factors
#   geoinfo$country[atsea.pts] <- "at sea"
#   
#   ### add geographic info to cuckoo dataset and redefine projection of newlongs/newlats
#   fullcuckoo <- data.frame(cuckoo, geoinfo)
#   coordinates(fullcuckoo) <- c("newlongs","newlats")
#   proj4string(fullcuckoo) <- CRS("+proj=longlat +datum=WGS84")
#   
#   ### create Europe only (+ at sea points) and Africa only data sets
#   Eurcuckoo <- fullcuckoo[which(fullcuckoo$continent=="Europe" | fullcuckoo$continent=="at sea"),]
#   Afrcuckoo <- fullcuckoo[which(fullcuckoo$continent=="Africa"),]
#   
#   if (nrow(Afrcuckoo@data) != 0) {
#     
#     ###------------------------------------------------------------###
#     #         Trim Africa data - stopovers only
#     ###------------------------------------------------------------###
#     
#     # subset out only stopovers in Africa (includes winter "stopovers")
#     
#     ### RESAMPLED CUCKOO DATA
#     
#     Afrcuckoo$mgroup <- as.factor(Afrcuckoo$mgroup)
#     stopovers <- Afrcuckoo[which(Afrcuckoo$mtype=="S" | is.na(Afrcuckoo$mtype)),]
#     
#     stopovers@data <- droplevels(stopovers@data)
#     
#     # can add other subsetting variables here, and change whatever "dat" is called
#     
#     # remove stopovers (mgroups) with < 2 transmission cycles
#     
#     # for each level of mgroup, check the number of levels of tcycle - if < 2, then Subset data to remove that mgroup
#     
#     newstopovers <- stopovers
#     
#     for (i in 1:length(levels(stopovers$mgroup))) {
#       temp <- newstopovers[newstopovers$mgroup==levels(stopovers$mgroup)[i],]
#       temp@data <- droplevels(temp@data)
#       if (length(levels(as.factor(temp$tcycle))) < 2) {
#         newstopovers <- newstopovers[which(newstopovers$mgroup!=levels(stopovers$mgroup)[i]),]
#         newstopovers@data <- droplevels(newstopovers@data)
#       } else {
#         newstopovers <- newstopovers
#       }
#     }
#     
#     newdat <- spTransform(newstopovers, CRS("+init=epsg:3395"))
#     
#     # write African data to list, to write to csv later
#     
#     Afrdata[[a]] <- newdat
#     
#     ###------------------------------------------------------------###
#     #         Plot individual cuckoo stopovers on Africa map 
#     ###------------------------------------------------------------###
#     
#     
#     
#     
#     ###### set directory for export of all to one folder
#     
#     outputdir <- c("C:/Users/samf/Documents/Sam's Stuff/R/projects/cuckoos/output/Africa stopovers/")
#     setwd(outputdir)
#     
#     # export Africa stopovers to a tiff file, with different colours for different years
#     
#     if (length(levels(as.factor(newdat@data$year)))==1) {
#       colpoints <- factor(newdat@data$year, labels=c("black"))  
#     } else if (length(levels(as.factor(newdat@data$year)))==2) {
#       colpoints <- factor(newdat@data$year, labels=c("darkmagenta","black"))  
#     } else {
#       colpoints <- factor(newdat@data$year, labels=c("royalblue","darkmagenta","black"))  
#     }
#     
#     #   tiff(paste(levels(newdat$name)," - Africa.tiff"),res=300,height=1000,width=1200,units="px")
#     #   
#     #   par(mar=c(1,1,1,1))
#     #   plot(Afrmap)
#     #   plot(newdat, pch=16, cex=0.8, col=as.character(colpoints), add=T)
#     #   coords <- cbind(-10,0)
#     #   textpts <- SpatialPoints(coords)
#     #   proj4string(textpts) <- CRS("+proj=longlat +datum=WGS84")
#     #   text(textpts, levels(newdat$name))
#     #   
#     #   
#     #   dev.off()
#     
#     ###-----------------------------------------------###
#     #         Draw MCP around stopover points
#     ###-----------------------------------------------###
#     
#     # for each movement group of each cuckoo, draw a MCP around the mgroup's points using mcp() to estimate home range
#     # choose points for MCP such that all points are included (percent = 100)
#     newdat@data <- droplevels(newdat@data)
#     stoppoints <- newdat[,c("mgroup")]
#     MCP <- mcp(stoppoints, percent=100)
#     
#     MCP@data <- data.frame(name=rep(levels(newdat$name),nrow(MCP@data)),MCP@data)
#     colnames(MCP@data) <- c("name","mgroupid","MCParea")
#     
#     AfrMCPstopovers[[a]] <- MCP
#     
#     #   writeOGR(MCP, dsn = ".", layer = paste(levels(newdat$name)," - Africa stopover MCPs", sep=""), driver = "ESRI Shapefile")
#     #   list.files(pattern = "stopover MCPs")
#   } else {}
#   
# }
# 
# AfrMCPstopovers2 <- AfrMCPstopovers[!unlist(lapply(AfrMCPstopovers, is.null))]
# 
# Afrdata.complete <- Afrdata[!unlist(lapply(Afrdata, is.null))]
# 
# Afrdata.full <- do.call(rbind, Afrdata.complete)
# 
# write.csv(Afrdata.full, "Africa resampled location points.csv", row.names=F)
# 
# for (i in 1:length(AfrMCPstopovers2)){
#   AfrMCPstopovers2[[i]] <- spChFIDs(AfrMCPstopovers2[[i]], paste(AfrMCPstopovers2[[i]]$name, AfrMCPstopovers2[[i]]$mgroupid))
# }
# 
# allAfrstopovers <- do.call(rbind, AfrMCPstopovers2)
# 
# allAfrstopovers2 <- spTransform(allAfrstopovers, CRS("+init=epsg:4326"))
# 
# writeOGR(allAfrstopovers2, dsn = ".", layer = "all birds - Africa stopover MCPs - 20131223 - m2", driver = "ESRI Shapefile")
# 
# #############################
# #############################
# #############################
# 
# 
