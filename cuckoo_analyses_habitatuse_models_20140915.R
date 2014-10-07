##########################################################
#
#   CUCKOO ANALYSES: 1) OCCUPANCY MODEL
#
#  Samantha Franks
#  11 April 2014
#  24 May 2014
#  17 Jul 2014: new random polygons generating 10 pseudoabsences/presence
#  5 Aug 2014: analysis using randomly generated points instead of polygons for absences, still 10 absent stopovers per presence
#  10 Sep 2014: added country data to absence points and added country data to at sea points
#
##########################################################

######## NOTES ON ANALYSIS ########

###--- ANALYSIS 1 - OCCUPANCY MODEL ---###

### DO HABITAT, SIZE OF PROTECTED AREA, LEVEL OF PROTECTED AREA PROTECTION, CLIMATE, AND ELEVATION PREDICT CUCKOO PRESENCE/ABSENCE?

### MODEL STRUCTURE ###

# logistic regression, sampling unit = points, nested random effects of individual(stopover(pointid))

# presence/absence (sampling unit = the stopover MCP) ~
#     HABITAT (5 land cover categories)
#     PROTECTED AREA (Y/N)
#     LEVEL OF PROTECTION (national or international)
#     ELEVATION
#     WINTER CLIMATE
#     SPRING CLIMATE
#     MIGRATION STRATEGY (E/W)

#### REVISED ANALYSIS 5 AUG 2014 - sampling unit = individual simulated points ####
# 
# A. PSEUDOREPLICATION
#
# Using points instead of stopover polygons introduces pseudoreplication at various levels:
#   
# 1) number of points in a stopover: stopovers with long duration and/or lots of points will introduce spatial bias by oversampling certain areas
# - control by including mgroup as a nested random effect within individual
# - length of stopover may in fact be due to the presence or quality of particular habitat
# 
# 2) resampled points: points simulated by incorporating location error are spatially autocorrelated
# - control by including point id as a nested random effect within mgroup
# 
#
# **B. WEIGHT OBSERVATIONS IN MODEL ACCORDING TO LOCATION UNCERTAINTY - Rob Rob thinks this is unnecessary because we've already accounted for the error in the simulation bootstrap - also, you can't do this in R (weights argument is when you need to indicate the number of trials that leads to proportional response data) but apparently you can in SAS
# - weight = 1/area of the error ellipse of a point
#
#
# C. HABITAT VARIABLES ARE CATEGORICAL
# - dummy variables created for 5 different habitat categories
#
#
# D. SPEI CLIMATE VARIABLES
# - 6-month SPEI values are calculated based on the current month's plus 5 previous months data
# - this timescale should provide a good representation of vegetation growing conditions over time periods likely to be relevant for determining the quality of habitats that cuckoos use
# - SPEI-month values to use will be March and August
# - SPEI-Mar represents the conditions over the winter (Oct-Mar period), which is of greatest relevance in the Mediterranean where winter rainfall is most important for vegetation growth
# - SPEI-Aug represents the conditions over the spring and summer (Mar-Aug period), which is of greater relevance in the more northern parts of Europe where rainfall early in the growing season sets the conditions for vegetation growth
# - from Vicente-Serrano et al. (2012 PNAS) a 6 month SPEI may be the best trade-off in using a timescale that will be relatively well correlated with vegetation activity (as measured by NDVI and other indices)
#
#
# D. MODEL VARIABLES & CANDIDATE MODEL SET
# - Response = presence/absence
# - Protected area (categorical): does point fall in a protected area, Y/N
# - Winter climate (continuous): SPEI 6-month SPEI in March, which summarizes the relative precipitation/evapotranspiration that has occurred over the last 6 months (Oct-Mar) compared to the norm (standardized against the background 30 year dataset)
# - spring climate (continuous): SPEI 6-month SPEI in August, which summarizes the relative precipitation/evapotranspiration that has occurred over the last 6 months (Mar-Aug) compared to the norm (standardized against the background 30 year dataset)
# - elevation

#### COMMAND ARGUMENTS ####

# passed from PATTERN in Condor script
# route (2 characters), spatial scale (3 characters), analysis type (2 characters)

# Analysis type:
# FM = full model, all data
# DL = designation level model
# SS = short stopover
# LS = long stopover
# FS = final stopover
# DW = designation level model for SW lumping N's with National PAs
# FR = France only
# SP = Spain only
# IT = Italy only
# DE = Germany + Low countries only

Mac <- FALSE

if(.Platform$OS =='windows') cluster <- FALSE
if(.Platform$OS=='unix' & !Mac) cluster <- TRUE

if (cluster) {
  Args <- commandArgs()[3]
  
  route <- substr(Args,1,2)	# SW/SE
  randomradius <- as.numeric(substr(Args,3,5)) # 50, 100, 200, or 500km
  analysistype <- substr(Args,6,7)
}

if (!cluster) {
  route <- "DE"
  randomradius <- 200
  analysistype <- "DE"
}

if (cluster) {
  library(plyr)
  library(reshape)
  #library(glmmML)
  library(lme4)
  #library(AICcmodavg)
  library(arm)
}

if (!cluster) {
  library(plyr)
  library(reshape)
  #library(glmmML)
  library(lme4)
  #library(AICcmodavg)
  library(arm)
  library(ggplot2)
  library(scales)
  library(sp)
  library(rgeos)
  library(rgdal)
}



###----------------------------------------------###
####         SET WORKING DIRECTORIES		####
###----------------------------------------------###


if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  
}

  
# BTO HPC cluster
# if (cluster) parentwd <- c("/users1/samf/cuckoos")

# Wales HPC cluster
if (cluster) parentwd <- c("/home/samantha.franks/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")
workspacewd <- paste(parentwd, "/workspaces", sep="")


####==== IMPORT POINT DATA ====####

setwd(paste(datawd, "/data for analysis/", sep=""))
present <- read.csv("presence data all variables points all with country data.csv", header=T)
absent <- read.csv(paste("absence data all variables points all with country data ", randomradius, " km.csv", sep=""), header=T)


### without at sea country data and correct country points for absences

# setwd(paste(datawd, "/data for analysis/", sep=""))
# present2 <- read.csv("presence data all variables points.csv", header=T)
# absent2 <- read.csv(paste("absence data all variables points ", randomradius, " km.csv", sep=""), header=T)

# present <- rename(present, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))
# 
# absent <- rename(absent, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))



####################################################
#
####   SIMPLE ANALYSIS EXPLORATIONS - PLOTS		####
#
####################################################

# setwd(paste(parentwd, "/scripts/", sep=""))
# source("cuckoo_analyses_data exploration_20140805.R")


###########################################################################
#
####    LOGISTIC REGRESSION MIXED MODEL ANALYSIS & MODEL SELECTION     ####
#
###########################################################################

#### NOTES ####
# Check for co-linearity among continuous explanatory variables using pairplots pairs() (habitat features, climate, elevation), and cor(x, method="spearman") 
#   **** co-linearity is less of an issue now because habitat is coded as a categorical variable

### ----------------------------- ###
####        DATA PREPARATION     ####
### ------------------------------###

#### ---- CORRECT COLUMN NAMES OF DATASET SO CAN MERGE ---- ####

newpresent <- rename(present, c("newlongs.epsg3035"="long.epsg3035", "newlats.epsg3035"="lat.epsg3035", "newlongs.epsg4326"="long.epsg4326", "newlats.epsg4326"="lat.epsg4326"))
newabsent <- rename(absent, c("randomlong.epsg3035"="long.epsg3035", "randomlat.epsg3035"="lat.epsg3035", "randomlong.epsg4326"="long.epsg4326", "randomlat.epsg4326"="lat.epsg4326"))

#### ---- ADD STOPOVER DURATION TO DATA ---- ####

### creates 2 files: present.LOS and absent.LOS with length of stay data added onto present and absent datasets
setwd(paste(parentwd, "/scripts/", sep=""))
source("source code to add stopover duration data for analysis.R")

############# subset out only one round of absences
#newabsent <- subset(newabsent, nabsence==1)
newabsent <- newabsent[,-which(names(newabsent) %in% c("nabsence","random.scale"))]
absent.LOS <- absent.LOS[,-which(names(absent.LOS) %in% c("nabsence","random.scale"))]

present.LOS <- present.LOS[,-which(names(present.LOS) %in% c("datetime"))]
present.LOS <- present.LOS[c(2,1,3:length(present.LOS))]

alldata <- rbind(newpresent, newabsent)
alldata.LOS <- rbind(present.LOS, absent.LOS)

# plotLOSdat <- alldata.LOS[alldata.LOS$presence==1,c("name","mgroup","LOS")]
# 
# plotLOSdat2 <- unique(plotLOSdat)
# 
# hist(plotLOSdat2$LOS, breaks=35)
# 
# plotdistdat <- alldata.LOS[alldata.LOS$presence==1, c("name","mgroup","id","distance")]
# plotdistdat <- unique(plotdistdat)
# 
# hist(plotdistdat$distance[plotdistdat$distance<50], breaks=50)

#### ----add 3-level designation level variable ---- ####

alldata.LOS$desig <- factor(alldata.LOS$desiglevel, levels=c(levels(alldata.LOS$desiglevel), "N"))
alldata.LOS$desig[which(is.na(alldata.LOS$desiglevel))] <- "N"


#### ---- ADD NEW 'STRATEGY' VARIABLE, GROUPED BY COUNTRY ---- ####

strategy2 <- factor(levels=c("SW","SE"))
strategy2[which(alldata.LOS$usecountry=="France" | alldata.LOS$usecountry=="Spain")] <- "SW"
strategy2[-which(alldata.LOS$usecountry=="France" | alldata.LOS$usecountry=="Spain")] <- "SE"

alldata.withnewstrategy <- data.frame(alldata.LOS, strategy2)

#### ---- MODEL DATASET ---- ####

fulldat <- with(alldata.withnewstrategy, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, country=usecountry, desig, laststop, Sahara.success, LOS))

# if full dataset model
if (analysistype=="FM") {
  dat <- subset(fulldat, strat==route)
  dat <- droplevels(dat)
}

# if long stopover
if (analysistype=="LS"){
  longstops <- subset(fulldat, LOS >= 10)
  longstops <- droplevels(longstops)
  dat <- subset(longstops, strat==route)
  dat <- droplevels(dat)
}


# if short stopover
if (analysistype=="SS"){
  shortstops <- subset(fulldat, LOS < 10)
  shortstops <- droplevels(shortstops)
  dat <- subset(shortstops, strat==route)
  dat <- droplevels(dat) 
}

# if final stopover
if (analysistype=="FS"){
  ## add a variable for whether Sahara crossing was attempted or not (ie. remove birds that failed before attempting)
  setwd(datawd)
  othervar <- read.csv("cuckoo migratory strategy and Sahara crossing success.csv", header=TRUE)
  
  attemptSahara <- with(othervar, data.frame(name, attemptSahara=attempted.Sahara.crossing))
  attemptdat <- merge(fulldat, attemptSahara)
  
  finalstops <- subset(attemptdat, laststop=="Y")
  finalstops <- droplevels(finalstops)
  
  dat <- subset(finalstops, strat==route)
  dat <- droplevels(dat) 
}

# if testing designation level, call desig variable "PA" so same models can be used
if (analysistype=="DL") {
  dat <- subset(fulldat, strat==route)
  dat <- droplevels(dat)
  
  dat <- rename(fulldat, c("PA"="PAoverlap", "desig"="PA"))
}

# if testing designation level for SW, lump designation level = N with desig level = National and call it No Desigation
if (analysistype=="DW") {
  fulldat$newdesig <- fulldat$desig
  fulldat$newdesig[which(fulldat$desig=="National")] <- "N"
  fulldat <- droplevels(fulldat)
  
  SWdat <- subset(fulldat, strat==route)
  SWdat <- droplevels(SWdat)
  
  dat <- rename(SWdat, c("PA"="PAoverlap", "newdesig"="PA"))
}

# if testing France only
if (analysistype=="FR") {
  dat <- subset(fulldat, country=="France")
  dat <- droplevels(dat)
}

# if testing Spain only
if (analysistype=="SP") {
  dat <- subset(fulldat, country=="Spain")
  dat <- droplevels(dat)
}

# if testing Italy only
if (analysistype=="IT") {
  dat <- subset(fulldat, country=="Italy")
  dat <- droplevels(dat)
}

# if testing Germany/Low Countries only
if (analysistype=="DE") {
  dat <- subset(fulldat, country=="Germany" | country=="Belgium" | country=="Netherlands")
  dat <- droplevels(dat)
}

############################################################
############################################################
############################################################
############################################################

### ----------------------------------------------------------------- ###
####      MODELS - PROTECTED AREA AND HABITAT OVERALL MODELS         ####
### ----------------------------------------------------------------- ###

###--- NOTES ---###
# Reacher's single stopover is likely driving a strong interaction between PA and unsuitable habitat
# to change the habitat reference level from agriculture to unsuitable - a reference level of agriculture probably makes sense in Europe since that is the predominant land cover type
# dat$habitat <- factor(dat$habitat, levels=c("unsuitable","agriculture","scrub.grassland","wetland.water","forest"))

# add a description of which set of models is running to the output
cat("\n\n\n\n\n################################### Starting to run ", analysistype, " models for ", randomradius, " km ", route, " data ###################################\n\n\n\n\n")

m <- list()

  # NULL MODEL
  system.time({
    m[[1]] <- glmer(presence ~ 1 + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  })
  
  ####------------------- PROTECTED AREA MODELS -------------------####
  
  # PA
  # PA + elevation
  # PA + elevation + winter climate
  # PA + elevation + spring climate
  # PA + elevation + winter climate + spring climate
  # PA*elevation
  # PA*elevation + winter climate
  # PA*elevation + spring climate
  # PA*elevation + winter climate + spring climate
  
  ### main effects only
  system.time({
    m[[2]] <- glmer(presence ~ PA + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[3]] <- glmer(presence ~ PA + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[4]] <- glmer(presence ~ PA + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[5]] <- glmer(presence ~ PA + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[6]] <- glmer(presence ~ PA + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  })
  
  # ### interaction between PA x elevation
  # m[[7]] <- glmer(presence ~ PA*elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  # m[[8]] <- glmer(presence ~ PA*elevation + winter climate + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  # m[[9]] <- glmer(presence ~ PA*elevation + spring climate + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  # m[[10]] <- glmer(presence ~ PA*elevation + winter climate + spring climate + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  ####------------------- HABITAT MODELS -------------------####
  # habitat
  # habitat + elevation
  # habitat + elevation + winterclimate
  # habitat + elevation + springclimate
  # habitat + elevation + winterclimate + springclimate
  # habitat*elevation
  # habitat*elevation + winterclimate
  # habitat*elevation + springclimate
  # habitat*elevation + winterclimate + springclimate
  # habitat*elevation
  # habitat*elevation + winterclimate
  # habitat*elevation + springclimate
  # habitat*elevation + winterclimate + springclimate
  
  ### main effects only
  system.time({
    m[[7]] <- glmer(presence ~ habitat + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[8]] <- glmer(presence ~ habitat + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[9]] <- glmer(presence ~ habitat + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[10]] <- glmer(presence ~ habitat + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[11]] <- glmer(presence ~ habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  })
   
  # ### interaction between habitat x elevation
  # m[[16]] <- glmer(presence ~ habitat*elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  # m[[17]] <- glmer(presence ~ habitat*elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  # m[[18]] <- glmer(presence ~ habitat*elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  # m[[19]] <- glmer(presence ~ habitat*elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  
  ####------------------- PROTECTED AREA +/x HABITAT MODELS -------------------####
  # PA + habitat
  # PA + habitat + elevation
  # PA + habitat + elevation + winterclimate
  # PA + habitat + elevation + springclimate
  # PA + habitat + elevation + winterclimate + springclimate
  # PA*elevation
  # PA*elevation + winter climate
  # PA*elevation + spring climate
  # PA*elevation + winter climate + spring climate
  
  ### main effects only
  system.time({
    m[[12]] <- glmer(presence ~ PA + habitat + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[13]] <- glmer(presence ~ PA + habitat + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[14]] <- glmer(presence ~ PA + habitat + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[15]] <- glmer(presence ~ PA + habitat + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[16]] <- glmer(presence ~ PA + habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  })
  
  
  ### interaction between PA x habitat
  system.time({
    m[[17]] <- glmer(presence ~ PA*habitat + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[18]] <- glmer(presence ~ PA*habitat + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[19]] <- glmer(presence ~ PA*habitat + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[20]] <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[21]] <- glmer(presence ~ PA*habitat + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  })
  
  ####------------------- ELEVATION & CLIMATE MODELS ONLY -------------------####
  system.time({
    m[[22]] <- glmer(presence ~ elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[23]] <- glmer(presence ~ elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[24]] <- glmer(presence ~ spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[25]] <- glmer(presence ~ spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[26]] <- glmer(presence ~ elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[27]] <- glmer(presence ~ elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
    m[[28]] <- glmer(presence ~ spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  })


setwd(workspacewd)
rm(list=(ls()[-which(ls()=="m" | ls()=="dat" | ls()=="route" | ls()=="randomradius" | ls()=="analysistype")]))
save.image(paste(analysistype, " models ", randomradius, " km ", route, ".RData", sep=""))

cat("\n\n\n\n\n################################### Finished running ", analysistype, " models for ", randomradius, " km ", route, " data ###################################\n", "###################################\n###################################\n###################################\n\n\n\n\n")


#####################################
#####################################
#####################################
#####################################