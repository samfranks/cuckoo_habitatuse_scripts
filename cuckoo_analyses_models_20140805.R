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
library(ggplot2)
library(scales)


###----------------------------------------------###
####         SET WORKING DIRECTORIES		####
###----------------------------------------------###


points <- TRUE

cluster <- FALSE

Mac <- FALSE

randomradius <- 50

if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  
}

  
if (cluster) parentwd <- c("/users1/samf/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

#extractionwd <- c("/corine PA elevation spei extracted values/")

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


####===========	PROTECTED AREA AND HABITAT USE EFFECT ON PRESENCE/ABSENCE ==============####

#### NOTES ####
# Check for co-linearity among continuous explanatory variables using pairplots pairs() (habitat features, climate, elevation), and cor(x, method="spearman") 
#   **** co-linearity is less of an issue now because habitat is coded as a categorical variable
# Have not yet included designation level as an explanatory variable since analysis would need to be conducted on a subset of data (only stopovers which overlap with PAs); otherwise NAs are introduced and models do not include the same number of observations)
# Initial analysis to include all sites, subsequent analysis may investigate the effect of designation level on use and interaction with habitat using only sites which intersect with Protected Areas

### ----------------------------- ###
####        DATA PREPARATION     ####
### ------------------------------###

#### ---- CHOOSE LANDSCALE SCALE FOR RANDOM DATASET ---- ####

newpresent <- rename(present, c("newlongs.epsg3035"="long.epsg3035", "newlats.epsg3035"="lat.epsg3035", "newlongs.epsg4326"="long.epsg4326", "newlats.epsg4326"="lat.epsg4326"))
newabsent <- rename(absent, c("randomlong.epsg3035"="long.epsg3035", "randomlat.epsg3035"="lat.epsg3035", "randomlong.epsg4326"="long.epsg4326", "randomlat.epsg4326"="lat.epsg4326"))

############# subset out only one round of absences
#newabsent <- subset(newabsent, nabsence==1)
newabsent <- newabsent[,-which(names(newabsent) %in% c("nabsence","random.scale"))]

alldata <- rbind(newpresent, newabsent)


#### ---- MODEL DATASET ---- ####

fulldat <- with(alldata, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS))

datSE <- subset(fulldat, strat=="SE")
datSE <- droplevels(datSE)

datSW <- subset(fulldat, strat=="SW")
datSW <- droplevels(datSW)

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


### ---------------------------------------- ###
####                 MODELS                 ####
### ------------------------ ----------------###
### change to fulldat, SE, or SW

SW <- TRUE

if (SW) dat <- datSW
if (!SW) dat <- datSE

### removing Reacher, who's single stopover is driving a strong interaction between PA and unsuitable habitat
# dat <- subset(dat, name!="Reacher")

# to change the habitat reference level from agriculture to unsuitable - a reference level of agriculture probably makes sense in Europe since that is the predominant land cover type
# dat$habitat <- factor(dat$habitat, levels=c("unsuitable","agriculture","scrub.grassland","wetland.water","forest"))

m <- list()

# NULL MODEL
m[[1]] <- glmer(presence ~ 1 + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

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
m[[2]] <- glmer(presence ~ PA + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[3]] <- glmer(presence ~ PA + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[4]] <- glmer(presence ~ PA + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[5]] <- glmer(presence ~ PA + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[6]] <- glmer(presence ~ PA + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

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
m[[7]] <- glmer(presence ~ habitat + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[8]] <- glmer(presence ~ habitat + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[9]] <- glmer(presence ~ habitat + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[10]] <- glmer(presence ~ habitat + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[11]] <- glmer(presence ~ habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# ### interaction between habitat x elevation
# m[[16]] <- glmer(presence ~ habitat*elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# m[[17]] <- glmer(presence ~ habitat*elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# m[[18]] <- glmer(presence ~ habitat*elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# m[[19]] <- glmer(presence ~ habitat*elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))


####------------------- PROTECTED AREA + HABITAT MODELS -------------------####
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
m[[12]] <- glmer(presence ~ PA + habitat + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[13]] <- glmer(presence ~ PA + habitat + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[14]] <- glmer(presence ~ PA + habitat + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[15]] <- glmer(presence ~ PA + habitat + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[16]] <- glmer(presence ~ PA + habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

### interaction between PA x habitat
m[[17]] <- glmer(presence ~ PA*habitat + elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[18]] <- glmer(presence ~ PA*habitat + elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[19]] <- glmer(presence ~ PA*habitat + elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[20]] <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[21]] <- glmer(presence ~ PA*habitat + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# SWmodels[[22]] <- glmer(presence ~ elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# SEmodels[[22]] <- glmer(presence ~ elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# ### interaction between PA x habitat x elevation
# m[[29]] <- glmer(presence ~ PA*habitat*elevation + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# m[[30]] <- glmer(presence ~ PA*habitat*elevation + spei.Mar + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# m[[31]] <- glmer(presence ~ PA*habitat*elevation + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# m[[32]] <- glmer(presence ~ PA*habitat*elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))


# # 
# saveobjects <- which(ls() == "SEmodels" | ls() == "SWmodels")

# #allPAobjects <- which(ls() == "PA.subset" | ls() == "r" | ls() == "Europeraster" | ls() == "corine.crs")

# toremove <- ls()[-saveobjects]

# rm(list=toremove)

####================= MODEL OUTPUT =================####

setwd(paste(parentwd, "/scripts/", sep=""))
source("source code cuckoo model selection.R")

### ---------------------------------------- ###
####              PLOT OUTPUT               ####
### ------------------------ ----------------###

### Ben Bolker's function for calculating CIs on predictions from an merMod object and plotting the results from his RPubs GLMM worked examples
# http://rpubs.com/bbolker/glmmchapter
# by specifying re.form=NA we're saying that we want the population-level prediction, i.e. setting the random effects to zero and getting a prediction for an average (or unknown) group
# Computing confidence intervals on the predicted values is relatively easy if we're willing to completely ignore the random effects, and the uncertainty of the random effects

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
  
### SW cuckoos

newdat <- data.frame(habitat=levels(datSW$habitat), PA=rep(levels(datSW$PA), 5), elevation=mean(datSW$elevation), spei.Mar=mean(datSW$spei.Mar), spei.Aug=mean(datSW$spei.Aug))

pred <- predict(SWmodels[[20]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(SWmodels[[20]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (SW) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(axis.title.x = element_text(face="bold", size=16, vjust=0.1), axis.text.x = element_text(size=12), axis.title.y = element_text(face="bold", size=16, vjust=0.9), axis.text.y = element_text(size=12, vjust=0.5), legend.text=element_text(size=14), legend.title=element_text(size=14))

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("top model interaction plot 50km SW with error bars.jpg")

# interaction.plot(newdat$habitat, newdat$PA, newdat$prob, xlab="Habitat", ylab=c"Predicted Probability", trace.label=c("PA"))
#title("Predicted probability of presence (SW) at mean elevation and mean climate values, 200km")

### SE cuckoos
newdat <- data.frame(habitat=levels(datSE$habitat), PA=rep(levels(datSE$PA), 5), elevation=mean(datSE$elevation), spei.Mar=mean(datSE$spei.Mar), spei.Aug=mean(datSE$spei.Aug))

pred <- predict(SEmodels[[20]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(SEmodels[[20]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (SE) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(axis.title.x = element_text(face="bold", size=16, vjust=0.1), axis.text.x = element_text(size=12), axis.title.y = element_text(face="bold", size=16, vjust=0.9), axis.text.y = element_text(size=12, vjust=0.5), legend.text=element_text(size=14), legend.title=element_text(size=14))

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("top model interaction plot 50km SE with error bars.jpg")

#interaction.plot(newdat$habitat, newdat$PA, newdat$prob, xlab=c("Habitat"), ylab=c("Predicted Probability"), trace.label=c("PA"))
#title("Predicted probability of presence (SE) at mean elevation and mean climate values, 200km")


########################################################################
########################################################################
########################################################################


####===========	DESIGNATION LEVEL ON PRESENCE/ABSENCE IN PROTECTED AREAS ==============####

#### ---- MODEL DATASET ---- ####

fulldat <- with(alldata, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel))


PAdat <- data.frame(fulldat, desig=alldata$desiglevel)
PAdat$desig <- factor(PAdat$desig, levels=c(levels(PAdat$desig), "N"))
PAdat$desig[which(is.na(PAdat$desig))] <- "N"

# PAdat <- subset(fulldat, PA=="Y")
# PAdat <- droplevels(PAdat)

datSE <- subset(PAdat, strat=="SE")
datSE <- droplevels(datSE)

datSW <- subset(PAdat, strat=="SW")
datSW <- droplevels(datSW)

dat <- datSW

desig.SWmodels[[1]] <- glmer(presence ~ desig*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

desig.SEmodels[[1]] <- glmer(presence ~ desig*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# dm2 <- glmer(presence ~ desig + habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# dm3 <- glmer(presence ~ desig + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# dm4 <- glmer(presence ~ habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# dm5 <- glmer(presence ~ 1 + (1|name/mgroup/id), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

#desig.SWmodels <- list()
#desig.SWmodels[[1]] <- dm1

#desig.SEmodels <- list()
#desig.SEmodels[[1]] <- dm1

relgrad <- with(desig.SEmodels[[1]]@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

setwd(paste(outputwd, "/GLMM results", sep=""))

sink(paste("SW GLMM 3-level PA designation models ", randomradius, " km.txt", sep=""))
cat("\n########==========  SW MODEL SUMMARIES ==========########\n", sep="\n")
print(summary(desig.SWmodels[[1]]))
sink()

sink(paste("SE GLMM 3-level PA designation models ", randomradius, " km.txt", sep=""))
cat("\n########==========  SW MODEL SUMMARIES ==========########\n", sep="\n")
print(summary(desig.SEmodels[[1]]))
sink()

### ---------------------------------------- ###
####              PLOT OUTPUT               ####
### ------------------------ ----------------###

### SE cuckoos
newdat <- data.frame(habitat=levels(datSE$habitat), desig=rep(levels(datSE$desig), 5), elevation=mean(datSE$elevation), spei.Mar=mean(datSE$spei.Mar), spei.Aug=mean(datSE$spei.Aug))

newdat$prob <- predict(dm1, newdata=newdat, type="response", re.form=NA)

interaction.plot(newdat$habitat, newdat$desig, newdat$prob, xlab=c("Habitat"), ylab=c("Predicted Probability"), trace.label=c("Designation"))
title("Predicted probability of presence (SE) at mean elevation and mean climate values, 50km")

### SW cuckoos
newdat <- data.frame(habitat=levels(datSW$habitat), desig=rep(levels(datSW$desig), 5), elevation=mean(datSW$elevation), spei.Mar=mean(datSW$spei.Mar), spei.Aug=mean(datSW$spei.Aug))

newdat$prob <- predict(dm1, newdata=newdat, type="response", re.form=NA)

interaction.plot(newdat$habitat, newdat$desig, newdat$prob, xlab=c("Habitat"), ylab=c("Predicted Probability"), trace.label=c("Designation"))
title("Predicted probability of presence (SW) at mean elevation and mean climate values, 50km")

