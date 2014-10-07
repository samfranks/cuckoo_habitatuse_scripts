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
#library(ggplot2)
#library(scales)


###----------------------------------------------###
####         SET WORKING DIRECTORIES		####
###----------------------------------------------###


HPCcluster <- TRUE

Mac <- FALSE

if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  
}

  
if (HPCcluster) parentwd <- c("cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

#extractionwd <- c("/corine PA elevation spei extracted values/")

####==== IMPORT POINT DATA ====####

  
setwd(paste(datawd, "/data for analysis/", sep=""))
present <- read.csv("presence data all variables points.csv", header=T)
# absent <- read.csv(paste("absence data all variables points ", randomradius, " km.csv", sep=""), header=T)

present <- rename(present, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))

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

### ----------------------------- ###
####        DATA PREPARATION     ####
### ------------------------------###

#### ---- CHOOSE LANDSCALE SCALE FOR RANDOM DATASET ---- ####

# newpresent <- rename(present, c("newlongs.epsg3035"="long.epsg3035", "newlats.epsg3035"="lat.epsg3035", "newlongs.epsg4326"="long.epsg4326", "newlats.epsg4326"="lat.epsg4326"))
# newabsent <- rename(absent, c("randomlong.epsg3035"="long.epsg3035", "randomlat.epsg3035"="lat.epsg3035", "randomlong.epsg4326"="long.epsg4326", "randomlat.epsg4326"="lat.epsg4326"))
# 
# #### ADD STOPOVER DURATION TO DATA ####
# 
# ### creates 2 files: present.LOS and absent.LOS with length of stay data added onto present and absent datasets
# # setwd(paste(parentwd, "/scripts/", sep=""))
# # source("source code to add stopover duration data for analysis.R")
# 
# ############# subset out only one round of absences
# #newabsent <- subset(newabsent, nabsence==1)
# newabsent <- newabsent[,-which(names(newabsent) %in% c("nabsence","random.scale"))]
# #absent.LOS <- absent.LOS[,-which(names(absent.LOS) %in% c("nabsence","random.scale"))]
# 
# #present.LOS <- present.LOS[,-which(names(present.LOS) %in% c("datetime"))]
# #present.LOS <- present.LOS[c(2,1,3:length(present.LOS))]
# 
# alldata <- rbind(newpresent, newabsent)
# 
# #alldata.LOS <- rbind(present.LOS, absent.LOS)


########################################################################
########################################################################
########################################################################


####===========  SAHARA CROSSING ==============####

#### ---- MODEL DATASET ---- ####


### convert Sahara success yes/no into integer 1/0

saharasuc <- subset(present)

saharasuc$Sahara.success <- as.numeric(as.character(revalue(saharasuc$Sahara.success, c("Y"=1, "N"=0))))

### add a variable for whether Sahara crossing was attempted or not (ie. remove birds that failed before attempting)
setwd(datawd)
othervar <- read.csv("cuckoo migratory strategy and Sahara crossing success.csv", header=TRUE)

attemptSahara <- with(othervar, data.frame(name, attemptSahara=attempted.Sahara.crossing))

attemptdat <- merge(saharasuc, attemptSahara)

saharadat <- subset(attemptdat, attemptSahara=="Y")
saharadat <- droplevels(saharadat)

saharadat <- with(saharadat, data.frame(success=Sahara.success, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel, country, laststop, attemptSahara))

#### RUN MODELS ####

saharamodels <- list()

saharamodels[[1]] <- glmer(success ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=saharadat, family=binomial, control=glmerControl(optimizer="bobyqa"))

saharamodels[[2]] <- glmer(success ~ PA + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=saharadat, family=binomial, control=glmerControl(optimizer="bobyqa"))

saharamodels[[3]] <- glmer(success ~ habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=saharadat, family=binomial, control=glmerControl(optimizer="bobyqa"))

saharamodels[[4]] <- glmer(success ~ elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=saharadat, family=binomial, control=glmerControl(optimizer="bobyqa"))

lapply(saharamodels, summary)

newdat <- data.frame(elevation=seq(min(saharadat$elevation), max(saharadat$elevation), by=0.001), spei.Mar=mean(saharadat$spei.Mar), spei.Aug=mean(saharadat$spei.Aug))
# 
# 
# tapply(saharasuc$elevation.m, list(as.factor(saharasuc$Sahara.success)), mean)
# ggplot(saharadat, aes(x=elevation, y=success))
# 
# p <- ggplot(byhab.country, aes(x=habitat, y=proportion, fill=presence)) + geom_bar(stat="identity", position="dodge")
# countrylabels <- labs(x="Country", y="Proportion of points per habitat", title="Proportion of points per habitat, by country", size=4)
# countrytheme <- theme(
#   axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
#   axis.text.x = element_text(size=12, angle=45, vjust=0.7), 
#   axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
#   axis.text.y = element_text(size=12, vjust=0.5),
#   legend.text=element_text(size=14), 
#   legend.title=element_text(size=14),
#   plot.title = element_text(face="bold", size=18)
# )
# p + countrylabels + countrytheme + scale_fill_manual(values=c("#0066CC", "#003366")) + facet_wrap(~ country, ncol=3)
# 
### Ben Bolker's population-level GLMM CI prediction function
easyPredCI <- function(model,newdata,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
                    newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link (logistic) function: could also use plogis()
  linkinv <- model@resp$family$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(lwr=pred0-crit*pred.se,
                upr=pred0+crit*pred.se)) 
}

pred <- predict(saharamodels[[4]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(saharamodels[[4]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)
# 
# ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of successful Sahara crossing", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(
#   axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
#   axis.text.x = element_text(size=12), 
#   axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
#   axis.text.y = element_text(size=12, vjust=0.5),
#   strip.text.x = element_text(size=14, face="bold"),
#   legend.text=element_text(size=14), 
#   legend.title=element_text(size=14)
# )
# 
# setwd(paste(outputwd, "/GLMM results/", sep=""))
# ggsave("final stopovers SW interaction plot 50km with error bars.jpg")