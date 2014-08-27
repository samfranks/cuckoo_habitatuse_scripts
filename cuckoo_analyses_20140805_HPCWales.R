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

#### REVISED ANALYSIS 5 AUG 2014 - sampling unit = individual simulated points ####
# 
# A. PSEUDOREPLICATION
#
# Using points instead of polygons introduces pseudoreplication at various levels:
#   
# 1) number of points in a stopover: stopovers with long duration and/or lots of points will introduce spatial bias by oversampling certain areas
# - control by including mgroup as a nested random effect within individual
# - length of stopover may in fact be due to the presence or quality of particular habitat
# 
# 2) resampled points: points simulated by incorporating location error are spatially autocorrelated
# - control by including point id as a nested random effect within mgroup
# 
#
# B. WEIGHT OBSERVATIONS IN MODEL ACCORDING TO LOCATION UNCERTAINTY
# - weight = 1/area of the error ellipse of a point
#
#
# C. HABITAT VARIABLES ARE CATEGORICAL
# - dummy variables created for 5 different habitat categories

###--- ANALYSIS 1 - OCCUPANCY MODEL ---###

### DO HABITAT, SIZE OF PROTECTED AREA, LEVEL OF PROTECTED AREA PROTECTION, CLIMATE, AND ELEVATION PREDICT CUCKOO PRESENCE/ABSENCE?

### MODEL STRUCTURE ###

# logistic regression, with individual as a random effect

# presence/absence (sampling unit = the stopover MCP) ~
#     HABITAT (proportion of 6 land cover categories)
#     PROTECTED AREA (Y/N)
#     (AREA PROTECTED AREA (amount of stopover MCP that is PA, by proportion or absolute area in km2))
#     LEVEL OF PROTECTION (national or international)
#     ELEVATION
#     WINTER RAINFALL (average local rainfall for grid cell Nov-Feb)
#     WINTER TEMP (average local temp for grid cell Nov-Feb)
#     SPRING RAINFALL (average local rainfall for grid cell Mar-Jun)
#     SPRING TEMP (average local temp for grid cell Mar-Jun)
#     MIGRATION STRATEGY (E/W)
#     random effect of stopover? or individual? or stopover nested within individual?

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

HPCWales <- TRUE

cluster <- FALSE

if (cluster) parentwd <- c("/users1/samf/cuckoos")
if (HPCWales) parentwd <- c("/home/rob.robinson/sam")

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
absent <- read.csv("absence data all variables points 50 km.csv", header=T)

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

# if (points) {
#   
#   corine.present <- prop.table(table(present$LAND.CLASS))
#   corine.absent <- prop.table(table(absent$LAND.CLASS))
#   corine.summary <- rbind(corine.present, corine.absent)
#   
#   ### barplot showing all random scales
#   par(mar=c(5,5,2,1))
#   
#   barplot(as.matrix(corine.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present","random 50 km"), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
#   
# }
# 
# 
# ####=========== HABITAT USE: POLYGONS ============####
# 
# if (!points) {
#   
#   ################# MULTIPLE RANDOM SCALES ###################
#   
#   meancorine.present <- aggregate(present[,c(habitatvarnames)], list(presence=present$presence), mean)
#   meancorine.present <- meancorine.present[,-1]
#   rownames(meancorine.present) <- "present"
#   
#   corinesummary.absent <- list()
#   for (i in 1:length(absentsplit)) {
#     absentsubset <- absentsplit[[i]]
#     corinesummary.absent[[i]] <- aggregate(absentsubset[,habitatvarnames], list(presence=absentsubset$presence), mean)
#   }
#   
#   names(corinesummary.absent) <- paste("absent.",names(absentsplit), sep="")
#   
#   meancorine.absent <- do.call(rbind,corinesummary.absent)
#   meancorine.absent <- meancorine.absent[,-1]
#   
#   meancorine <- rbind(meancorine.present, meancorine.absent)
#   
#   # meancorine <- reshape(corinesummary, times=c("agriculture","artificial","bare_land","forest","scrub_grassland","wetland_water"), timevar="landcover", varying=list(c("agriculture","artificial","bare_land","forest","scrub_grassland","wetland_water")), direction="long")
#   
#   ### barplot showing all random scales
#   par(mar=c(5,5,2,1))
#   
#   barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present","random 50 km", "random 100km", "random 200km", "random 500km"), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
#   
#   ### pooling across random scales ###
#   absentpooled <- apply(meancorine.absent, 2, mean)
#   meancorine <- rbind(meancorine.present, absentpooled)
#   rownames(meancorine) <- c("present", "absent")
#   
#   ### barplot showing pooled random scales
#   par(mar=c(5,5,2,1))
#   
#   barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present","random"), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
#   
#   ################# SINGLE RANDOM SCALE ###################
#   
#   meancorine.absent <- aggregate(absent[,c(habitatvarnames)], list(presence=absent$presence), mean)
#   meancorine.absent <- meancorine.absent[,-1]
#   rownames(meancorine.absent) <- "absent"
#   
#   meancorine <- rbind(meancorine.present, meancorine.absent)
#   
#   ### barplot showing single random scale
#   par(mar=c(5,5,2,1))
#   
#   barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present","random 50 km"), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
#   
#   #########################################################
#     
#   ###!!!!!!!!!!!!!!!!! CAUTIONARY NOTE !!!!!!!!!!!!!!!!!!!!!!###
#   
#   # mtext("random wetland-water not necessarily representative due to sampling methodology", cex=0.8, side=3, line=1)
#   
#   # more wetland/water than random, but random habitat extraction specifically tries to avoid stopover polygons with too much water, so possibly an artifact of random stopover sampling procedure
#   
#   # may need to re-run with a lower avoidance rate of wetland-water sites (sea only?), or increase the proportion to 0.6-0.7
#   
#   # latest revision of random stopover generation code to generate 10 pseudoabsences per absence allowed more wetland.water habitat (particularly inland marshes, peat bogs, and water courses) - 0.5 threshold for choosing absence stopovers without too much other water categories (marine and water bodies)
#   
# }
# 
# 
# ###-----------------------------------------###
# ####            PROTECTED AREAS            ####
# ###-----------------------------------------###
# 
# #### PROPORTION OF ALL POINTS THAT OVERLAP WITH PROTECTED AREA - POINTS ####
# 
# if (points) {
#   
#   overlap.present <- prop.table(table(present$PAoverlap))
#   overlap.absent <- prop.table(table(absent$PAoverlap))
#   overlap.summary <- rbind(overlap.present, overlap.absent)
#   
#   par(mar=c(5,5,2,1))
#   barplot(as.matrix(overlap.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("N", "Y"), cex.names=0.8, xpd=TRUE, legend=c("present","random 50 km"), args.legend=list(bty="n", cex=0.8), xlab="In a protected area?", ylab="proportion of all stopovers \n overlapping protected areas")
#   
# }
# 
# #### PROPORTION OF ALL STOPOVER POLYGONS THAT OVERLAP WITH PROTECTED AREA ####
# 
# if (!points) {
#   
#   overlap.present <- prop.table(table(present$overlap))
#   overlap.absent <- prop.table(table(absent$overlap))
#   overlap.summary <- rbind(overlap.present, overlap.absent)
#   
#   par(mar=c(5,5,2,1))
#   barplot(as.matrix(overlap.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("N", "Y"), cex.names=0.8, xpd=TRUE, legend=c("present","random 50 km"), args.legend=list(bty="n", cex=0.8), xlab="In a protected area?", ylab="proportion of all points \n in protected areas")
#   
#   ####==== Proportion by area of stopover that overlaps with PA overlap ====####
#   
#   meanpropPA.present <- mean(present$prop.overlap)
#   
#   summarypropPA.absent <- list()
#   for (i in 1:length(absentsplit)) {
#     absentsubset <- absentsplit[[i]]
#     summarypropPA.absent[[i]] <- mean(absentsubset$prop.overlap)
#   }
#   
#   names(summarypropPA.absent) <- paste("absent.",names(absentsplit), sep="")
#   
#   meanpropPA.absent <- do.call(rbind,summarypropPA.absent)
#   
#   meanpropPA <- rbind(meanpropPA.present, meanpropPA.absent)
#   rownames(meanpropPA) <- c("present",rownames(meanpropPA.absent))
#   
#   par(mar=c(5,5,2,1))
#   barplot(meanpropPA, beside=TRUE, space=0.1, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), names.arg=c("present","random 50 km", "random 100km", "random 200km", "random 500km"), cex.names=0.7, xpd=TRUE, ylab="mean proportion of stopover \n overlapping protected areas", ylim=c(0,0.4))
#   
#   ### pooled random stopovers ###
#   meanpropPA <- rbind(meanpropPA.present, mean(meanpropPA.absent))
#   
#   ### pooled barplot ###
#   par(mar=c(5,5,2,1))
#   barplot(meanpropPA, beside=TRUE, space=0.1, col=c("darkblue","lightblue"), names.arg=c("present","random"), cex.names=0.7, xpd=TRUE, ylab="mean proportion of stopover \n overlapping protected areas", ylim=c(0,0.4))
#   
# }
# 
# 
# ####==== Proportion of international vs national PA use ====####
# 
# #### POINTS ####
# 
# if (points) {
#   
#   x <- subset(present, PAoverlap=="Y")
#   y <- subset(absent, PAoverlap=="Y")
#   desig.P <- prop.table(table(x$desiglevel))
#   desig.A <- prop.table(table(y$desiglevel))
#   desig.summary <- rbind(desig.P, desig.A)
#   
#   par(mar=c(5,5,2,1), mfrow=c(1,1))
#   barplot(as.matrix(desig.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("International", "National"), cex.names=0.8, xpd=TRUE, legend=c("present","random 50 km"), args.legend=list(bty="n", cex=0.8), xlab="Protected Area designation level", ylab="proportion of all points \n in protected areas")
#   
#   
# }
# 
# #### POLYGONS ####
# 
# if (!points) {
#   
#   #### MULTIPLE SCALES ###
#   
#   desigprop.P <- prop.table(table(present$overlap, present$desig_type))[2,]
#   
#   desigprop.A <- list()
#   for (i in 1:4) {
#     absentsubset <- absentsplit[[i]]
#     overlapsubset <- subset(absentsubset, overlap=="Y")
#     desigprop.A[[i]] <- prop.table(table(overlapsubset$desig_type))
#   }
#   
#   designation <- rbind(desigprop.P, do.call(rbind, desigprop.A))
#   rownames(designation) <- c("present","absent.50","absent.100","absent.200","absent.500")
#   
#   par(mar=c(5,5,2,8), xpd=TRUE)
#   barplot(designation, beside=TRUE, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), xlab="Protected area designation level", ylab="Proportion of overlapping protected \n area at each designation level")
#   legend("right", xpd=TRUE, c("present","random 50 km", "random 100km", "random 200km", "random 500km"), bty="n", cex=0.8, fill=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), inset=c(-0.3,-0.2))
#   
#   ### pooled random stopovers ###
#   designation <- rbind(desigprop.P, apply(do.call(rbind, desigprop.A), 2, mean))
#   
#   ### pooled random stopovers barplot ###
#   par(mar=c(5,5,2,5), xpd=TRUE)
#   barplot(designation, beside=TRUE, col=c("darkblue","lightblue"), xlab="Protected area designation level", ylab="Proportion of overlapping protected \n area at each designation level", ylim=c(0,0.7))
#   legend("right", xpd=TRUE, c("present","random"), bty="n", cex=0.8, fill=c("darkblue","lightblue"), inset=c(-0.15,-0.2))
#   
# }
# 
# 
# 
# ###--------------------------------###
# ####   PROTECTED AREAS & HABITAT ####
# ###--------------------------------###
# 
# if (points) {
#   
#   PAhab.present <- prop.table(table(present$PAoverlap, present$LAND.CLASS))
#   PAhab.absent <- prop.table(table(absent$PAoverlap, absent$LAND.CLASS))
#   overlap.summary <- rbind(PAhab.present, PAhab.absent)
#   
#   par(mar=c(5,5,2,1), mfrow=c(1,2))
#   barplot(as.matrix(PAhab.present), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("Outside PA","Inside PA"), args.legend=list(bty="n", cex=0.8), xlab="Land cover class", ylab="proportion of all points")
#   
#   barplot(as.matrix(PAhab.absent), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("Outside PA","Inside PA"), args.legend=list(bty="n", cex=0.8), xlab="Land cover class", ylab="proportion of all points")
#   
#   
# }


#################################################################
#################################################################
#################################################################

###############################################
#
####    LOGISTIC REGRESSION ANALYSIS       ####
#
###############################################

#### NOTES ####
# Check for co-linearity among continuous explanatory variables using pairplots pairs() (habitat features, climate, elevation), and cor(x, method="spearman") 
#   **** co-linearity is less of an issue now because habitat is coded as a categorical variable
# Have not yet included designation level as an explanatory variable since analysis would need to be conducted on a subset of data (only stopovers which overlap with PAs); otherwise NAs are introduced and models do not include the same number of observations)
# Initial analysis to include all sites, subsequent analysis may investigate the effect of designation level on use and interaction with habitat using only sites which intersect with Protected Areas

### ---------------------------------------- ###
####        DATA PREPARATION and PCA        ####
### ------------------------ ----------------###

#### ---- CHOOSE LANDSCALE SCALE FOR RANDOM DATASET ---- ####

### models to be run using the 50km random absent data
# with unsuitable category (combines artificial, bareland, water)

# absent.data <- subset(absent, random.scale==50, select=-c(nabsence,random.scale))
newpresent <- rename(present, c("newlongs.epsg3035"="long.epsg3035", "newlats.epsg3035"="lat.epsg3035", "newlongs.epsg4326"="long.epsg4326", "newlats.epsg4326"="lat.epsg4326"))
newabsent <- rename(absent, c("randomlong.epsg3035"="long.epsg3035", "randomlat.epsg3035"="lat.epsg3035", "randomlong.epsg4326"="long.epsg4326", "randomlat.epsg4326"="lat.epsg4326"))
newabsent <- newabsent[,-which(names(newabsent) %in% c("nabsence","random.scale"))]

alldata <- rbind(newpresent, newabsent)

### ALTERNATIVE RANDOM ABSENT DATASET ###
# instead of using absences from the 50km landscape scale, could pool across all landscape levels? if pooled, can't ask specific questions about habitat choice at the landscape level (it may defer at the different scales)

#### ---- COLLINEARITY ---- ####

# if (!points) {
#   
#   ### explore collinearity of continuous explanatory variables
#   exp.cont <- alldata[,c("elevation.m", "spei.Mar","spei.Aug")]
#   
#   #exp.cont <- alldata[,c("agriculture","forest","scrub_grassland","unsuitable","precip.mm.winter","temp.C.winter","precip.mm.spring","temp.C.spring","elevation.m")]
#   
#   source(paste(parentwd,"/scripts/Zuur_functions.r", sep=""))
#   # for 6 habitat categories
#   
#   pairs(exp.cont, lower.panel=panel.cor, panel=panel.smooth)
#   
#   # pairs(exp.cont[,7:11], lower.panel=panel.cor, panel=panel.smooth)
#   # pairs(exp.cont[,1:6], lower.panel=panel.cor, panel=panel.smooth)
#   # pairs(exp.cont[,1:11], lower.panel=panel.cor, panel=panel.smooth)
#   # 
#   # # for 4 habitat categories
#   # pairs(exp.cont[,5:9], lower.panel=panel.cor, panel=panel.smooth)
#   # pairs(exp.cont[,1:4], lower.panel=panel.cor, panel=panel.smooth)
#   # pairs(exp.cont[,1:9], lower.panel=panel.cor, panel=panel.smooth)
#   
#   # unsurprisingly, winter and spring temperatures are highly correlated (r = 0.75), so should drop one of these - will drop WINTER TEMP, but keep winter rainfall
#   # agriculture is unsurprisingly negative correlated with forest and scrub_grassland, but might leave this for the moment (correlations are < 0.7)
#   # see Zuur Mixed Effects Models Ch 21 GLMM and koala distribution for possible remedies for correlated habitat variables
#   
#   
# }

#### ---- MODEL DATASET ---- ####

#dat <- with(alldata, data.frame(presence, name, strat=strategy, hab.ag=agriculture, hab.for=forest, hab.scrub=scrub_grassland, hab.un=unsuitable, w.temp=temp.C.winter, w.precip=precip.mm.winter, s.precip=precip.mm.spring, s.temp=temp.C.spring, elevation=elevation.m, PA.cat=overlap, PA.cont=prop.overlap))
# dat <- with(alldata, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, error.major, error.minor, error.weight=1/(pi*error.major*error.minor), rescale.weights=rescale(1/(pi*error.major*error.minor))))

dat <- with(alldata, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS))

# 
# error.weight <- with(dat, pi*error.major*error.minor)
# error.weight.small <- error.weight/10000
# 
# inverse.sqrt.errors <- 1/sqrt(error.weight)
# inverse.log.errors <- 1/log(error.weight)
# 
# log.errors <- log(error.weight)
# sqrt.errors <- sqrt(error.weight)
# 
# summary(error.weight/10000)
# 
# hist(error.weight/10000)
# 
# m2 <- glmer(presence ~ PA + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"), weights=error.weight.small)
# m2.noweights <- glmer(presence ~ PA + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))


### ---------------------- ###
####        MODELS        ####
### -----------------------###

#### ---- CANDIDATE MODEL SET ---- ####

# presence/absence ~

# # protected area models
# strategy + propPA
# strategy + propPA + elevation
# strategy + propPA + elevation + winterclimate
# strategy + propPA + elevation + springclimate
# strategy + propPA + elevation + winterclimate + springclimate 

# # habitat models
# strategy + habitat
# strategy + habitat + elevation
# strategy + habitat + elevation + winterclimate
# strategy + habitat + elevation + springclimate
# strategy + habitat + elevation + winterclimate + springclimate

# # protected area + habitat models
# strategy + propPA + habitat
# strategy + propPA + habitat + elevation
# strategy + propPA + habitat + elevation + winterclimate
# strategy + propPA + habitat + elevation + springclimate
# strategy + propPA + habitat + elevation + winterclimate + springclimate

# RANDOM EFFECT: all models include nested random effects of individual
# ????? do I need some sort of spatial random effect e.g. geographic region needs to be nested within individual??????

####################################


# modelnames <- c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","m13","m14","m15","m16")

# NULoverlap MODEL
m1 <- glmer(presence ~ 1 + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# # protected area models
# strategy + propPA
# strategy + propPA + elevation
# strategy + propPA + elevation + spei.Mar
# strategy + propPA + elevation + spei.Aug
# strategy + propPA + elevation + spei.Mar + spei.Aug

m2 <- glmer(presence ~ PA + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"), method="ML")
m3 <- glmer(presence ~ PA + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m4 <- glmer(presence ~ PA + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m5 <- glmer(presence ~ PA + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m6 <- glmer(presence ~ PA + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# # habitat models
# strategy + habitat
# strategy + habitat + elevation
# strategy + habitat + elevation + winterclimate
# strategy + habitat + elevation + springclimate
# strategy + habitat + elevation + winterclimate + springclimate

m7 <- glmer(presence ~ habitat + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m8 <- glmer(presence ~ habitat + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m9 <- glmer(presence ~ habitat + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m10 <- glmer(presence ~ habitat + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m11 <- glmer(presence ~ habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))



# # protected area + habitat models
# strategy + propPA + habitat
# strategy + propPA + habitat + elevation
# strategy + propPA + habitat + elevation + winterclimate
# strategy + propPA + habitat + elevation + springclimate
# strategy + propPA + habitat + elevation + winterclimate + springclimate

m12 <- glmer(presence ~ PA + habitat + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m13 <- glmer(presence ~ PA + habitat + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m14 <- glmer(presence ~ PA + habitat + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m15 <- glmer(presence ~ PA + habitat + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m16 <- glmer(presence ~ PA + habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# # interaction model
m17 <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))



models <- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17)
AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17)
modellist <- lapply(models, formula)


calculate.AIC<-function(aictable,modellist) {
  modelnames<-modellist
  delta.aic<-aictable$AIC-min(aictable$AIC)
  lik.aic<-exp(-delta.aic/2)
  aic.w<-lik.aic/(sum(lik.aic))
  aic.table<-data.frame(modelnames,aictable,delta.aic,lik.aic,aic.w)
}

output <- calculate.AIC(AICoutput, as.character(modellist))
aic.ordered<-output[rev(order(output$aic.w)),]

# see AIC table
aic.ordered[,3:6] <- round(aic.ordered[,c(3:6)], digits = 3)
print(aic.ordered)
# 
# 
# 
# #
# 
# #### ---- TEST MODELS ---- ####
# 
# modelnames <- c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","m13","m14","m15","m16")
# 
# # NULL MODEL
# assign(modelnames[1], glmer(presence ~ 1 + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# # # protected area models
# # strategy + propPA
# # strategy + propPA + elevation
# # strategy + propPA + elevation + spei.Mar
# # strategy + propPA + elevation + spei.Aug
# # strategy + propPA + elevation + spei.Mar + spei.Aug
# 
# assign(modelnames[2], glmer(presence ~ PA.cont + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[3], glmer(presence ~ strat + PA.cont + rescale.elev + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[4], glmer(presence ~ strat + PA.cont + rescale.elev + spei.Mar + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[5], glmer(presence ~ strat + PA.cont + rescale.elev + spei.Aug + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[6], glmer(presence ~ strat + PA.cont + rescale.elev + spei.Mar + spei.Aug + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# # # habitat models
# # strategy + habitat
# # strategy + habitat + elevation
# # strategy + habitat + elevation + winterclimate
# # strategy + habitat + elevation + springclimate
# # strategy + habitat + elevation + winterclimate + springclimate
# 
# # using PCA for the habitat variables, pc1=arable vs natural, pc2=forest vs grassland
# assign(modelnames[7], glmer(presence ~ strat + pc1 + pc2 + hab.un + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[8], glmer(presence ~ strat + pc1 + pc2 + hab.un + hab.wet + rescale.elev + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[9], glmer(presence ~ strat + pc1 + pc2 + hab.un + hab.wet + rescale.elev + spei.Mar + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[10], glmer(presence ~ strat + pc1 + pc2 + hab.un + hab.wet + rescale.elev + spei.Aug + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[11], glmer(presence ~ strat + pc1 + pc2 + hab.un + hab.wet + rescale.elev + spei.Mar + spei.Aug + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# 
# # # protected area + habitat models
# # strategy + propPA + habitat
# # strategy + propPA + habitat + elevation
# # strategy + propPA + habitat + elevation + winterclimate
# # strategy + propPA + habitat + elevation + springclimate
# # strategy + propPA + habitat + elevation + winterclimate + springclimate
# 
# assign(modelnames[12], glmer(presence ~ strat + PA.cont + pc1 + pc2 + hab.un + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[13], glmer(presence ~ strat + PA.cont + pc1 + pc2 + hab.un + hab.wet + rescale.elev + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[14], glmer(presence ~ strat + PA.cont + pc1 + pc2 + hab.un + hab.wet + rescale.elev + spei.Mar + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[15], glmer(presence ~ strat + PA.cont + pc1 + pc2 + hab.un + hab.wet + rescale.elev + spei.Aug + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[16], glmer(presence ~ strat + PA.cont + pc1 + pc2 + hab.un + hab.wet + rescale.elev + spei.Mar + spei.Aug + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# models <- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
# AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
# modellist <- lapply(models, formula)
# 
# 
# assign(modelnames[1], glmer(presence ~ strat + PA.cont + pc1 + pc2 + hab.un + hab.wet + rescale.elev + spei.Mar + spei.Aug + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# 
# calculate.AIC<-function(aictable,modellist) {
#   modelnames<-modellist
#   delta.aic<-aictable$AIC-min(aictable$AIC)
#   lik.aic<-exp(-delta.aic/2)
#   aic.w<-lik.aic/(sum(lik.aic))
#   aic.table<-data.frame(modelnames,aictable,delta.aic,lik.aic,aic.w)
# }
# 
# output <- calculate.AIC(AICoutput, as.character(modellist))
# aic.ordered<-output[rev(order(output$aic.w)),]
# 
# # see AIC table
# aic.ordered[,3:6] <- round(aic.ordered[,c(3:6)], digits = 3)
# print(aic.ordered)
# 
# 
# 
# # # test HABITAT variables to include
# assign(modelnames[1], glmer(presence ~ strat + pc1 + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[2], glmer(presence ~ strat + pc2 + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[3], glmer(presence ~ strat + hab.un + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[4], glmer(presence ~ strat + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# assign(modelnames[5], glmer(presence ~ strat + pc1 + pc2 + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[6], glmer(presence ~ strat + pc1 + hab.un + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[7], glmer(presence ~ strat + pc1 + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[8], glmer(presence ~ strat + pc2 + hab.un + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[9], glmer(presence ~ strat + pc2 + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[10], glmer(presence ~ strat + hab.un + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# assign(modelnames[11], glmer(presence ~ strat + pc1 + pc2 + hab.un + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[12], glmer(presence ~ strat + pc1 + pc2 + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[13], glmer(presence ~ strat + pc1 + hab.un + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# assign(modelnames[14], glmer(presence ~ strat + pc2 + hab.un + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# assign(modelnames[15], glmer(presence ~ strat + pc1 + pc2 + hab.un + hab.wet + (1|name), data=newdat, family=binomial, control=glmerControl(optimizer="bobyqa")))
# 
# models <- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
# AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
# modellist <- lapply(models, formula)
# 
# 
# calculate.AIC<-function(aictable,modellist) {
#   modelnames<-modellist
#   delta.aic<-aictable$AIC-min(aictable$AIC)
#   lik.aic<-exp(-delta.aic/2)
#   aic.w<-lik.aic/(sum(lik.aic))
#   aic.table<-data.frame(modelnames,aictable,delta.aic,lik.aic,aic.w)
# }
# 
# output <- calculate.AIC(AICoutput, as.character(modellist))
# aic.ordered<-output[rev(order(output$aic.w)),]
# 
# # see AIC table
# aic.ordered[,3:6] <- round(aic.ordered[,c(3:6)], digits = 3)
# print(aic.ordered)
# 
# ######################################################
# ######################################################
# ######################################################
# ######################################################
# 
# m1 <- glm(presence ~ prop.overlap, data=alldata, family=binomial)
# cdplot(presence ~ prop.overlap, data=alldata)
# 
# m2 <- glm(presence ~ overlap, data=alldata)
# ### Number of stopovers which overlap with a protected area
# table(present$overlap)
# for (z in 1:length(absentsplit)){
#   print(names(absentsplit[z]))
#   print(table(absentsplit[[z]]$overlap))
# }
# 
# ### simple comparisions using melt/cast
# 
# melt.mcp <- melt(mcp.all[,4:10], id=c("laststop","strategy", "Sahara.success","random"))
# 
# #setwd(outputwd)
# #sink("MCP protected area summary.txt")
# cast(melt.mcp, random ~ variable, mean)
# cast(melt.mcp, random ~ laststop ~ variable, mean)
# #sink()
# 
# boxplot(mcp.all$prop.MCP.overlap~mcp.all$random, xlab=paste("random stopover at \n", rradius, " radius", sep=""), ylab="", cex=0.8, cex.lab=0.8)
# 
# 
# #####################################
# 
# mcp.pa <- read.csv("allbirds MCP overlap data.csv", header=T)
# mcp.random <- read.csv(paste("random_", rradius, "_allbirds MCP overlap data.csv", sep=""), header=T)
# 
# mcp.pa2 <- data.frame(mcp.pa, random=rep("N", nrow(mcp.pa)))
# mcp.random2 <- data.frame(mcp.random, random=rep("Y", nrow(mcp.random)))
# 
# mcp.all <- rbind(mcp.pa2, mcp.random2)
# 
# ### simple comparisions using melt/cast
# 
# melt.mcp <- melt(mcp.all[,4:10], id=c("laststop","strategy", "Sahara.success","random"))
# 
# #setwd(outputwd)
# #sink("MCP protected area summary.txt")
# cast(melt.mcp, random ~ variable, mean)
# cast(melt.mcp, random ~ laststop ~ variable, mean)
# #sink()
# 
# boxplot(mcp.all$prop.MCP.overlap~mcp.all$random, xlab=paste("random stopover at \n", rradius, " radius", sep=""), ylab="", cex=0.8, cex.lab=0.8)
# 
# #, ylab="proportion of stopover MCP \n overlapping with protected area")
# 
# #######
# 
# setwd(mapoutputdir)
# 
# dev.off()
# 
# ### proportional area of all MCPs in protected areas
# mean(mcp.pa$prop.MCP.overlap)
# mean(mcp.random$prop.MCP.overlap)
# 
# aggregate(mcp.all$prop.MCP.overlap, by=list(random=mcp.all$random), mean)
# 
# ### proportional area of all laststop MCPs in protected areas
# aggregate(mcp.pa$prop.MCP.overlap, by=list(laststop=mcp.pa$laststop), mean)
# aggregate(mcp.random$prop.MCP.overlap, by=list(laststop=mcp.random$laststop), mean)
# 
# aggregate(mcp.all$prop.MCP.overlap, by=list(laststop=mcp.all$laststop, ))
# 
# 
# 
# ### proportional area of MCPs in protected areas for Sahara successful vs unsuccessful birds
# aggregate(mcp.pa$prop.MCP.overlap, by=list(Sahara.success=mcp.pa$Sahara.success), mean)
# x <- aggregate(mcp.pa$prop.MCP.overlap, by=list(Sahara.success=mcp.pa$Sahara.success, laststop=mcp.pa$laststop), mean)
# 
# ### size of protected area
# aggregate(mcp.pa$PA.overlap.area.km2, by=list(laststop=mcp.pa$laststop), mean)
# aggregate(mcp.pa$PA.overlap.area.km2, by=list(Sahara.success=mcp.pa$Sahara.success), mean)
# 
# ### migration strategy
# aggregate(mcp.pa$PA.overlap.area.km2, by=list(strategy=mcp.pa$strategy), mean)
# aggregate(mcp.pa$prop.MCP.overlap, by=list(strategy=mcp.pa$strategy), mean)
# 
# 
#   
#   ###-------------------------------------------###
#   #----  BARPLOTS OF CORINE HABITAT USE  ----
#   ###-------------------------------------------###  
#   
#   ### --- FUNCTION TO PLOT CUCKOO HABITAT USE, MULTIPLE SITES POOLED --- ###
#   
#   habsummary <- function(dataset,landcover){ # datasets required are point land cover values, clclegend data (landcover), and par settings for specific # of stopovers
#     
#     dataset$mgroup <- as.factor(as.character(dataset$mgroup))
#     
#     # calculate proportion of each habitat type used
#     prophab <- prop.table(table(dataset$LABEL4))
#     
#     # produce colour labels for barplot according to RGB colour codes of CLC legend
#     tmp1 <- landcover[which(landcover$LABEL4 %in% levels(dataset$LABEL4)),]
#     tmp2 <- tmp1[order(tmp1$LABEL4),]
#     landcolours <- rgb(tmp2[,c("rgb1","rgb2","rgb3")], max=255)
#     
#     # plot the proportional habitat use
#     
#     x <- barplot(prophab, names.arg=c(""), ylim=c(0,max(prophab)+0.05), col=landcolours, cex.axis=0.8, las=2, ylab="Proportion of total stopover area")
#     text(x, par("usr")[3] - 0.01, srt = 45, adj = 1, xpd = TRUE, cex=0.7, labels = levels(dataset$LABEL4))
#     
#   }
#     
#   
#   Chance29 <- subset(newactual, mgroup==29)
#   Chance29 <- droplevels(Chance29)
#   
#   par(mar=c(9,6,3,1), oma=c(0,0,2,0)) # specific par settings for number of plots to produce
#   habsummary(Chance29, clclegend)
#   
