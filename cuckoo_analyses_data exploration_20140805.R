##########################################################
#
#   CUCKOO ANALYSES: 1) OCCUPANCY MODEL, 2) SUCCESS/FAILURE DIFFERENCES
#
#  Samantha Franks
#  11 April 2014
#  24 May 2014
#  17 Jul 2014: new random polygons generating 10 pseudoabsences/presence
#  5 Aug 2014: analysis using randomly generated points instead of polygons for absences, still 10 absent stopovers per presence
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
# - spring climate (continuous): sPEI 6-month SPEI in August, which summarizes the relative precipitation/evapotranspiration that has occurred over the last 6 months (Mar-Aug) compared to the norm (standardized against the background 30 year dataset)
# - elevation

###--- ANALYSIS 2 - SUCCESS OF SAHARA CROSSING MODEL ---###

### DO HABITAT, SIZE OF PROTECTED AREA, LEVEL OF PROTECTED AREA PROTECTION, CLIMATE, AND ELEVATION OF THE *FINAL EUROPEAN STOPOVER* PREDICT SUCCESS/FAILURE OF SAHARA CROSSING?

# logistic regression, with individual as a random effect

# success/failure (sampling unit = the final stopover MCP) ~
#     HABITAT (proportion of 6 land cover categories)
#     AREA PROTECTED AREA (amount of stopover MCP that is PA, by proportion or absolute area in km2)
#     LEVEL OF PROTECTION (national or international)
#     ELEVATION
#     WINTER RAINFALL (average local rainfall for grid cell Nov-Feb)
#     WINTER TEMP (average local temp for grid cell Nov-Feb)
#     SPRING RAINFALL (average local rainfall for grid cell Mar-Jun)
#     SPRING TEMP (average local temp for grid cell Mar-Jun)
#     MIGRATION STRATEGY (E/W)
#     random effect of stopover? or individual? or stopover nested within individual?


library(plyr)
library(reshape)
#library(glmmML)
library(lme4)
#library(AICcmodavg)
library(arm)


###----------------------------------------------###
####         SET WORKING DIRECTORIES		####
###----------------------------------------------###


points <- TRUE

cluster <- FALSE

Mac <- TRUE

randomradius <- 200

if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  
}

  
if (cluster) parentwd <- c("/users1/samf/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

#extractionwd <- c("/corine PA elevation spei extracted values/")


####################################################
#
####   SIMPLE ANALYSIS EXPLORATIONS - PLOTS		####
#
####################################################

####==== IMPORT POINT DATA ====####

if (points) {
  
setwd(paste(datawd, "/data for analysis/", sep=""))
present <- read.csv("presence data all variables points.csv", header=T)
absent <- read.csv(paste("absence data all variables points ", randomradius, " km.csv", sep=""), header=T)

present <- rename(present, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))

absent <- rename(absent, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))

}


####==== IMPORT POLYGON DATA ====####

if (!points) {
  
  setwd(paste(datawd, "/data for analysis/", sep=""))
  present <- read.csv("presence data all variables.csv", header=T)
  absent <- read.csv("absence data all variables.csv", header=T)
  
  # # combining habitats artificial, bareland, water into "unsuitable" category
  # unsuitable <- rowSums(present[,c("artificial","bare_land","wetland_water")])
  # present.new <- data.frame(present[,1:12],unsuitable,present[,13:23])
  # present.new <- present.new[,-(c(8,9,12))]
  
  # unsuitable <- rowSums(absent[,c("artificial","bare_land","wetland_water")])
  # absent.new <- data.frame(absent[,1:12],unsuitable,absent[,13:24])
  # absent.new <- absent.new[,-(c(8,9,12))]
  
  # absentsplit <- split(absent, as.factor(absent$random.scale))
  # 
  
  habitatvarnames <- c("agriculture","forest","scrub.grassland","unsuitable","wetland.water")
  
}


###-----------------------------------------------------###
####   HABITAT USE - plot habitat types used vs random ####
###-----------------------------------------------------###

####=========== HABITAT USE: POINTS ============####

if (points) {
  
  corine.present <- prop.table(table(present$LAND.CLASS))
  corine.absent <- prop.table(table(absent$LAND.CLASS))
  corine.summary <- rbind(corine.present, corine.absent)
  
  ### barplot showing all random scales
  par(mar=c(5,5,2,1))
  
  barplot(as.matrix(corine.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present", paste("random ", randomradius, " km", sep="")), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
  
}


####=========== HABITAT USE: POLYGONS ============####

if (!points) {
  
  ################# MULTIPLE RANDOM SCALES ###################
  
  meancorine.present <- aggregate(present[,c(habitatvarnames)], list(presence=present$presence), mean)
  meancorine.present <- meancorine.present[,-1]
  rownames(meancorine.present) <- "present"
  
  corinesummary.absent <- list()
  for (i in 1:length(absentsplit)) {
    absentsubset <- absentsplit[[i]]
    corinesummary.absent[[i]] <- aggregate(absentsubset[,habitatvarnames], list(presence=absentsubset$presence), mean)
  }
  
  names(corinesummary.absent) <- paste("absent.",names(absentsplit), sep="")
  
  meancorine.absent <- do.call(rbind,corinesummary.absent)
  meancorine.absent <- meancorine.absent[,-1]
  
  meancorine <- rbind(meancorine.present, meancorine.absent)
  
  # meancorine <- reshape(corinesummary, times=c("agriculture","artificial","bare_land","forest","scrub_grassland","wetland_water"), timevar="landcover", varying=list(c("agriculture","artificial","bare_land","forest","scrub_grassland","wetland_water")), direction="long")
  
  ### barplot showing all random scales
  par(mar=c(5,5,2,1))
  
  barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present","random 50 km", "random 100km", "random 200km", "random 500km"), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
  
  ### pooling across random scales ###
  absentpooled <- apply(meancorine.absent, 2, mean)
  meancorine <- rbind(meancorine.present, absentpooled)
  rownames(meancorine) <- c("present", "absent")
  
  ### barplot showing pooled random scales
  par(mar=c(5,5,2,1))
  
  barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present","random"), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
  
  ################# SINGLE RANDOM SCALE ###################
  
  meancorine.absent <- aggregate(absent[,c(habitatvarnames)], list(presence=absent$presence), mean)
  meancorine.absent <- meancorine.absent[,-1]
  rownames(meancorine.absent) <- "absent"
  
  meancorine <- rbind(meancorine.present, meancorine.absent)
  
  ### barplot showing single random scale
  par(mar=c(5,5,2,1))
  
  barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present", paste("random ", randomradius, " km", sep="")), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
  
  #########################################################
    
  ###!!!!!!!!!!!!!!!!! CAUTIONARY NOTE !!!!!!!!!!!!!!!!!!!!!!###
  
  # mtext("random wetland-water not necessarily representative due to sampling methodology", cex=0.8, side=3, line=1)
  
  # more wetland/water than random, but random habitat extraction specifically tries to avoid stopover polygons with too much water, so possibly an artifact of random stopover sampling procedure
  
  # may need to re-run with a lower avoidance rate of wetland-water sites (sea only?), or increase the proportion to 0.6-0.7
  
  # latest revision of random stopover generation code to generate 10 pseudoabsences per absence allowed more wetland.water habitat (particularly inland marshes, peat bogs, and water courses) - 0.5 threshold for choosing absence stopovers without too much other water categories (marine and water bodies)
  
}


###-----------------------------------------###
####            PROTECTED AREAS            ####
###-----------------------------------------###

#### PROPORTION OF ALL POINTS THAT OVERLAP WITH PROTECTED AREA - POINTS ####

if (points) {
  
  overlap.present <- prop.table(table(present$PAoverlap))
  overlap.absent <- prop.table(table(absent$PAoverlap))
  overlap.summary <- rbind(overlap.present, overlap.absent)
  
  par(mar=c(5,5,2,1))
  barplot(as.matrix(overlap.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("N", "Y"), cex.names=0.8, xpd=TRUE, legend=c("present", paste("random ", randomradius, " km", sep="")), args.legend=list(bty="n", cex=0.8), xlab="In a protected area?", ylab="proportion of all stopovers \n overlapping protected areas")
  
}

#### PROPORTION OF ALL STOPOVER POLYGONS THAT OVERLAP WITH PROTECTED AREA ####

if (!points) {
  
  overlap.present <- prop.table(table(present$overlap))
  overlap.absent <- prop.table(table(absent$overlap))
  overlap.summary <- rbind(overlap.present, overlap.absent)
  
  par(mar=c(5,5,2,1))
  barplot(as.matrix(overlap.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("N", "Y"), cex.names=0.8, xpd=TRUE, legend=c("present", paste("random ", randomradius, " km", sep="")), args.legend=list(bty="n", cex=0.8), xlab="In a protected area?", ylab="proportion of all points \n in protected areas")
  
  ####==== Proportion by area of stopover that overlaps with PA overlap ====####
  
  meanpropPA.present <- mean(present$prop.overlap)
  
  summarypropPA.absent <- list()
  for (i in 1:length(absentsplit)) {
    absentsubset <- absentsplit[[i]]
    summarypropPA.absent[[i]] <- mean(absentsubset$prop.overlap)
  }
  
  names(summarypropPA.absent) <- paste("absent.",names(absentsplit), sep="")
  
  meanpropPA.absent <- do.call(rbind,summarypropPA.absent)
  
  meanpropPA <- rbind(meanpropPA.present, meanpropPA.absent)
  rownames(meanpropPA) <- c("present",rownames(meanpropPA.absent))
  
  par(mar=c(5,5,2,1))
  barplot(meanpropPA, beside=TRUE, space=0.1, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), names.arg=c("present","random 50 km", "random 100km", "random 200km", "random 500km"), cex.names=0.7, xpd=TRUE, ylab="mean proportion of stopover \n overlapping protected areas", ylim=c(0,0.4))
  
  ### pooled random stopovers ###
  meanpropPA <- rbind(meanpropPA.present, mean(meanpropPA.absent))
  
  ### pooled barplot ###
  par(mar=c(5,5,2,1))
  barplot(meanpropPA, beside=TRUE, space=0.1, col=c("darkblue","lightblue"), names.arg=c("present","random"), cex.names=0.7, xpd=TRUE, ylab="mean proportion of stopover \n overlapping protected areas", ylim=c(0,0.4))
  
}


####==== Proportion of international vs national PA use ====####

#### POINTS ####

if (points) {
  
  x <- subset(present, PAoverlap=="Y")
  y <- subset(absent, PAoverlap=="Y")
  desig.P <- prop.table(table(x$desiglevel))
  desig.A <- prop.table(table(y$desiglevel))
  desig.summary <- rbind(desig.P, desig.A)
  
  par(mar=c(5,5,2,1), mfrow=c(1,1))
  barplot(as.matrix(desig.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("International", "National"), cex.names=0.8, xpd=TRUE, legend=c("present", paste("random ", randomradius, " km", sep="")), args.legend=list(bty="n", cex=0.8), xlab="Protected Area designation level", ylab="proportion of all points \n in protected areas")
  
  
}

#### POLYGONS ####

if (!points) {
  
  #### MULTIPLE SCALES ###
  
  desigprop.P <- prop.table(table(present$overlap, present$desig_type))[2,]
  
  desigprop.A <- list()
  for (i in 1:4) {
    absentsubset <- absentsplit[[i]]
    overlapsubset <- subset(absentsubset, overlap=="Y")
    desigprop.A[[i]] <- prop.table(table(overlapsubset$desig_type))
  }
  
  designation <- rbind(desigprop.P, do.call(rbind, desigprop.A))
  rownames(designation) <- c("present","absent.50","absent.100","absent.200","absent.500")
  
  par(mar=c(5,5,2,8), xpd=TRUE)
  barplot(designation, beside=TRUE, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), xlab="Protected area designation level", ylab="Proportion of overlapping protected \n area at each designation level")
  legend("right", xpd=TRUE, c("present","random 50 km", "random 100km", "random 200km", "random 500km"), bty="n", cex=0.8, fill=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), inset=c(-0.3,-0.2))
  
  ### pooled random stopovers ###
  designation <- rbind(desigprop.P, apply(do.call(rbind, desigprop.A), 2, mean))
  
  ### pooled random stopovers barplot ###
  par(mar=c(5,5,2,5), xpd=TRUE)
  barplot(designation, beside=TRUE, col=c("darkblue","lightblue"), xlab="Protected area designation level", ylab="Proportion of overlapping protected \n area at each designation level", ylim=c(0,0.7))
  legend("right", xpd=TRUE, c("present","random"), bty="n", cex=0.8, fill=c("darkblue","lightblue"), inset=c(-0.15,-0.2))
  
}



###--------------------------------###
####   PROTECTED AREAS & HABITAT ####
###--------------------------------###

if (points) {
  
  PAhab.present <- prop.table(table(present$PAoverlap, present$LAND.CLASS))
  PAhab.absent <- prop.table(table(absent$PAoverlap, absent$LAND.CLASS))
  overlap.summary <- rbind(PAhab.present, PAhab.absent)
  
  par(mar=c(5,5,2,1), mfrow=c(1,2))
  barplot(as.matrix(PAhab.present), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("Present, Outside PA","Present, Inside PA"), args.legend=list(bty="n", cex=0.8), xlab="Land cover class", ylab="proportion of all points")
  
  barplot(as.matrix(PAhab.absent), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c(paste("Absent ", randomradius, " km, Outside PA", sep=""), paste("Absent ", randomradius , " km, Inside PA", sep="")), args.legend=list(bty="n", cex=0.8), xlab="Land cover class", ylab="proportion of all points")
  
  
}