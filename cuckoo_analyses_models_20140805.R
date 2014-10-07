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
Args <- commandArgs()

route <- substr(Args[3],1,2)	# SW/SE
randomradius <- as.numeric(substr(Args[3],3,5)) # 50, 100, 200, or 500km

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


Mac <- TRUE

if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  
}

  
if (cluster) parentwd <- c("/users1/samf/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")


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

### ----------------------------- ###
####        DATA PREPARATION     ####
### ------------------------------###

#### ---- CHOOSE LANDSCALE SCALE FOR RANDOM DATASET ---- ####

newpresent <- rename(present, c("newlongs.epsg3035"="long.epsg3035", "newlats.epsg3035"="lat.epsg3035", "newlongs.epsg4326"="long.epsg4326", "newlats.epsg4326"="lat.epsg4326"))
newabsent <- rename(absent, c("randomlong.epsg3035"="long.epsg3035", "randomlat.epsg3035"="lat.epsg3035", "randomlong.epsg4326"="long.epsg4326", "randomlat.epsg4326"="lat.epsg4326"))

#### ADD STOPOVER DURATION TO DATA ####

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
# hist(plotLOSdat2$LOS)
####===========	PROTECTED AREA AND HABITAT USE EFFECT ON PRESENCE/ABSENCE ==============####

#### NOTES ####
# Check for co-linearity among continuous explanatory variables using pairplots pairs() (habitat features, climate, elevation), and cor(x, method="spearman") 
#   **** co-linearity is less of an issue now because habitat is coded as a categorical variable
# Have not yet included designation level as an explanatory variable since analysis would need to be conducted on a subset of data (only stopovers which overlap with PAs); otherwise NAs are introduced and models do not include the same number of observations)
# Initial analysis to include all sites, subsequent analysis may investigate the effect of designation level on use and interaction with habitat using only sites which intersect with Protected Areas

#### ---- MODEL DATASET - reclassifying SE/SW according to country rather than departure route to Africa ---- ####


fulldat <- with(alldata.LOS, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, country=usecountry))

dat.route <- subset(fulldat, strat==route)
dat.route <- droplevels(dat.route)

# datSE <- subset(fulldat, strat=="SE")
# datSE <- droplevels(datSE)

# datSW <- subset(fulldat, strat=="SW")
# datSW <- droplevels(datSW)



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
ggsave("top model elevation plot 50km SW.jpg")

# interaction.plot(newdat$habitat, newdat$PA, newdat$prob, xlab="Habitat", ylab=c"Predicted Probability", trace.label=c("PA"))
#title("Predicted probability of presence (SW) at mean elevation and mean climate values, 200km")

############## ELEVATION SW

elevationseq <- seq(min(datSW$elevation), max(datSW$elevation), 0.01)
newdat <- data.frame(elevation=elevationseq, habitat=rep(levels(datSW$habitat), each=length(elevationseq)), PA=rep(levels(datSW$PA), each=length(elevationseq)*length(levels(datSW$habitat))) , spei.Mar=mean(datSW$spei.Mar), spei.Aug=mean(datSW$spei.Aug))

pred <- predict(SWmodels[[20]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(SWmodels[[20]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

qplot(x=elevation, y=pred, data=newdat2, geom="point", colour=habitat, shape=PA, size=2)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("top model elevation plot 50km SW.jpg")

############## climate SW
climateseq <- seq(min(datSW$spei.Aug), max(datSW$spei.Aug), 0.01)
newdat <- data.frame(spei.Aug=climateseq, habitat=rep(levels(datSW$habitat), each=length(climateseq)), PA=rep(levels(datSW$PA), each=length(climateseq)*length(levels(datSW$habitat))) , elevation=mean(datSW$elevation), spei.Mar=mean(datSW$spei.Mar))

pred <- predict(SWmodels[[20]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(SWmodels[[20]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

qplot(x=spei.Aug, y=pred, data=newdat2, geom="point", colour=habitat, shape=PA, size=2)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("top model summer SPEI plot 50km SW.jpg")

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


############## ELEVATION SE

elevationseq <- seq(min(datSE$elevation), max(datSE$elevation), 0.01)
newdat <- data.frame(elevation=elevationseq, habitat=rep(levels(datSE$habitat), each=length(elevationseq)), PA=rep(levels(datSE$PA), each=length(elevationseq)*length(levels(datSE$habitat))) , spei.Mar=mean(datSE$spei.Mar), spei.Aug=mean(datSE$spei.Aug))

pred <- predict(SEmodels[[20]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(SEmodels[[20]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

qplot(x=elevation, y=pred, data=newdat2, geom="point", colour=habitat, shape=PA, size=2)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("top model elevation plot 50km SE.jpg")


############## climate SE
climateseq <- seq(min(datSE$spei.Mar), max(datSE$spei.Mar), 0.01)
newdat <- data.frame(spei.Mar=climateseq, habitat=rep(levels(datSE$habitat), each=length(climateseq)), PA=rep(levels(datSE$PA), each=length(climateseq)*length(levels(datSE$habitat))) , elevation=mean(datSE$elevation), spei.Aug=mean(datSE$spei.Aug))

pred <- predict(SEmodels[[20]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(SEmodels[[20]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

qplot(x=spei.Mar, y=pred, data=newdat2, geom="point", colour=habitat, shape=PA, size=2)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("top model winter SPEI plot 50km SE.jpg")
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

### SW cuckoos
newdat <- data.frame(habitat=levels(datSW$habitat), desig=rep(levels(datSW$desig), 5), elevation=mean(datSW$elevation), spei.Mar=mean(datSW$spei.Mar), spei.Aug=mean(datSW$spei.Aug))

pred <- predict(desig.SWmodels[[1]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(desig.SWmodels[[1]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=desig, linetype=desig, group=desig, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (SW) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(axis.title.x = element_text(face="bold", size=16, vjust=0.1), axis.text.x = element_text(size=12), axis.title.y = element_text(face="bold", size=16, vjust=0.9), axis.text.y = element_text(size=12, vjust=0.5), legend.text=element_text(size=14), legend.title=element_text(size=14))

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("3-level PA designation interaction plot 50km SW with error bars.jpg")


### SE cuckoos
newdat <- data.frame(habitat=levels(datSE$habitat), desig=rep(levels(datSE$desig), 5), elevation=mean(datSE$elevation), spei.Mar=mean(datSE$spei.Mar), spei.Aug=mean(datSE$spei.Aug))

pred <- predict(desig.SEmodels[[1]], newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(desig.SEmodels[[1]], newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=desig, linetype=desig, group=desig, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (SE) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(axis.title.x = element_text(face="bold", size=16, vjust=0.1), axis.text.x = element_text(size=12), axis.title.y = element_text(face="bold", size=16, vjust=0.9), axis.text.y = element_text(size=12, vjust=0.5), legend.text=element_text(size=14), legend.title=element_text(size=14))

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("3-level PA designation interaction plot 50km SE with error bars.jpg")

########################################################################
########################################################################
########################################################################


####===========  COUNTRY AND PRESENCE/ABSENCE IN PROTECTED AREAS ==============####

#### ---- MODEL DATASET ---- ####

SPA.FRA <- subset(datSW2, country=="Spain" | country=="France")
SPA.FRA <- droplevels(SPA.FRA)

bigcountries <- subset(alldata, country=="Spain" | country=="France" | country=="Germany" | country=="Italy")
bigcountries <- droplevels(bigcountries)

countrydat <- with(bigcountries, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel, country))

countrymodel <- glmer(presence ~ country*PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=countrydat, family=binomial, control=glmerControl(optimizer="bobyqa"))

#big4countrymodel <- countrymodel

summary(countrymodel)

newdat <- data.frame(PA=rep(levels(countrydat$PA), each=5), country=rep(levels(countrydat$country), each=10), habitat=rep(levels(countrydat$habitat), 4), elevation=mean(countrydat$elevation), spei.Mar=mean(countrydat$spei.Mar), spei.Aug=mean(countrydat$spei.Aug))

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

pred <- predict(countrymodel, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(countrymodel, newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + facet_wrap( ~ country, ncol=2) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (stopovers in large countries) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(
    axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
    axis.text.x = element_text(size=12), 
    axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
    axis.text.y = element_text(size=12, vjust=0.5),
    strip.text.x = element_text(size=14, face="bold"),
    legend.text=element_text(size=14), 
    legend.title=element_text(size=14)
    )

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("large countries interaction plot 50km with error bars.jpg")

########################################################################
########################################################################
########################################################################


####===========  STOPOVER DURATION - LONG STOPOVERS ONLY ==============####

#### ---- MODEL DATASET ---- ####

longstops <- subset(alldata.LOS, LOS >= 10)
longstops <- droplevels(longstops)

longstopdat <- with(longstops, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel, country))

longstopdatSW <- subset(longstopdat, strat=="SW")
longstopdatSW <- droplevels(longstopdatSW)

longstopdatSE <- subset(longstopdat, strat=="SE")
longstopdatSE <- droplevels(longstopdatSE)

longstopmodelSW <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=longstopdatSW, family=binomial, control=glmerControl(optimizer="bobyqa"))

longstopmodelSE <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=longstopdatSE, family=binomial, control=glmerControl(optimizer="bobyqa"))


setwd(paste(outputwd, "/GLMM results", sep=""))

sink(paste("SW GLMM models long stopovers ", randomradius, " km.txt", sep=""))
cat("\n########==========  SW MODELS ==========########\n", sep="\n")
print(longstopmodelSW)
cat("\n##### CONVERGENCE WARNING degenerate Hessian with 2 negative eigenvalues\n")
cat("\n########==========  SW MODEL SUMMARIES ==========########\n", sep="\n")
print(summary(longstopmodelSW))

sink()

sink(paste("SE GLMM models long stopovers ", randomradius, " km.txt", sep=""))
cat("\n########==========  SE MODELS ==========########\n", sep="\n")
print(longstopmodelSE)
cat("\n########==========  SE MODEL SUMMARIES ==========########\n", sep="\n")
print(summary(longstopmodelSE))


sink()


#### SW

newdat <- data.frame(PA=rep(levels(longstopdatSW$PA), each=5), habitat=rep(levels(longstopdatSW$habitat)), elevation=mean(longstopdatSW$elevation), spei.Mar=mean(longstopdatSW$spei.Mar), spei.Aug=mean(longstopdatSW$spei.Aug))

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

pred <- predict(longstopmodelSW, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(longstopmodelSW, newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (LONG stopovers, SW) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(
  axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
  axis.text.x = element_text(size=12), 
  axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
  axis.text.y = element_text(size=12, vjust=0.5),
  strip.text.x = element_text(size=14, face="bold"),
  legend.text=element_text(size=14), 
  legend.title=element_text(size=14)
)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("long stopovers SW interaction plot 50km with error bars.jpg")


#### SE

newdat <- data.frame(PA=rep(levels(longstopdatSE$PA), each=5), habitat=rep(levels(longstopdatSE$habitat)), elevation=mean(longstopdatSE$elevation), spei.Mar=mean(longstopdatSE$spei.Mar), spei.Aug=mean(longstopdatSE$spei.Aug))

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

pred <- predict(longstopmodelSE, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(longstopmodelSE, newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (LONG stopovers, SE) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(
  axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
  axis.text.x = element_text(size=12), 
  axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
  axis.text.y = element_text(size=12, vjust=0.5),
  strip.text.x = element_text(size=14, face="bold"),
  legend.text=element_text(size=14), 
  legend.title=element_text(size=14)
)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("long stopovers SE interaction plot 50km with error bars.jpg")

########################################################################
########################################################################
########################################################################


####===========  STOPOVER DURATION - SHORT STOPOVERS ONLY ==============####

#### ---- MODEL DATASET ---- ####

shortstops <- subset(alldata.LOS, LOS < 10)
shortstops <- droplevels(shortstops)
shortstops$LOSlength <- "short"

shortstopdat <- with(shortstops, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel, country))

shortstopdatSW <- subset(shortstopdat, strat=="SW")
shortstopdatSW <- droplevels(shortstopdatSW)

shortstopdatSE <- subset(shortstopdat, strat=="SE")
shortstopdatSE <- droplevels(shortstopdatSE)

shortstopmodelSW <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=shortstopdatSW, family=binomial, control=glmerControl(optimizer="bobyqa"))

shortstopmodelSE <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=shortstopdatSE, family=binomial, control=glmerControl(optimizer="bobyqa"))

setwd(paste(outputwd, "/GLMM results", sep=""))

sink(paste("SW GLMM models short stopovers ", randomradius, " km.txt", sep=""))
cat("\n########==========  SW MODELS ==========########\n", sep="\n")
print(shortstopmodelSW)
cat("\n##### CONVERGENCE WARNING degenerate Hessian with 2 negative eigenvalues\n")
cat("\n########==========  SW MODEL SUMMARIES ==========########\n", sep="\n")
print(summary(shortstopmodelSW))

sink()

sink(paste("SE GLMM models short stopovers ", randomradius, " km.txt", sep=""))
cat("\n########==========  SE MODELS ==========########\n", sep="\n")
print(shortstopmodelSE)
cat("\n########==========  SE MODEL SUMMARIES ==========########\n", sep="\n")
print(summary(shortstopmodelSE))


sink()


#### SW

newdat <- data.frame(PA=rep(levels(shortstopdatSW$PA), each=5), habitat=rep(levels(shortstopdatSW$habitat)), elevation=mean(shortstopdatSW$elevation), spei.Mar=mean(shortstopdatSW$spei.Mar), spei.Aug=mean(shortstopdatSW$spei.Aug))

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

pred <- predict(shortstopmodelSW, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(shortstopmodelSW, newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (SHORT stopovers, SW) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(
  axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
  axis.text.x = element_text(size=12), 
  axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
  axis.text.y = element_text(size=12, vjust=0.5),
  strip.text.x = element_text(size=14, face="bold"),
  legend.text=element_text(size=14), 
  legend.title=element_text(size=14)
)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("short stopovers SW interaction plot 50km with error bars.jpg")


#### SE

newdat <- data.frame(PA=rep(levels(shortstopdatSE$PA), each=5), habitat=rep(levels(shortstopdatSE$habitat)), elevation=mean(shortstopdatSE$elevation), spei.Mar=mean(shortstopdatSE$spei.Mar), spei.Aug=mean(shortstopdatSE$spei.Aug))

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

pred <- predict(shortstopmodelSE, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(shortstopmodelSE, newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (SHORT stopovers, SE) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(
  axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
  axis.text.x = element_text(size=12), 
  axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
  axis.text.y = element_text(size=12, vjust=0.5),
  strip.text.x = element_text(size=14, face="bold"),
  legend.text=element_text(size=14), 
  legend.title=element_text(size=14)
)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("short stopovers SE interaction plot 50km with error bars.jpg")

########################################################################
########################################################################
########################################################################


####===========  LAST STOPOVERS ONLY - only birds that attempted Sahara crossing ==============####

#### ---- MODEL DATASET ---- ####

### Convert Sahara success yes/no into integer 1/0

#saharasuc <- subset(alldata.LOS, presence==1 & laststop=="Y")
saharasuc <- subset(alldata.LOS, presence==1)


saharasuc$Sahara.success <- as.numeric(as.character(revalue(saharasuc$Sahara.success, c("Y"=1, "N"=0))))

## add a variable for whether Sahara crossing was attempted or not (ie. remove birds that failed before attempting)
setwd(datawd)
othervar <- read.csv("cuckoo migratory strategy and Sahara crossing success.csv", header=TRUE)

attemptSahara <- with(othervar, data.frame(name, attemptSahara=attempted.Sahara.crossing))

attemptdat <- merge(saharasuc, attemptSahara)

attemptdat <- subset(attemptdat, attemptSahara=="Y")
attemptdat <- droplevels(attemptdat)

saharadat <- with(attemptdat, data.frame(success=Sahara.success, name, year, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=elevation.m, spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel, country, laststop, attemptSahara, LOS))

######################################################

## add a variable for whether Sahara crossing was attempted or not (ie. remove birds that failed before attempting)
setwd(datawd)
othervar <- read.csv("cuckoo migratory strategy and Sahara crossing success.csv", header=TRUE)

attemptSahara <- with(othervar, data.frame(name, attemptSahara=attempted.Sahara.crossing))

attemptdat <- merge(alldata, attemptSahara)

finalstops <- subset(attemptdat, laststop=="Y")
finalstops <- droplevels(finalstops)

finalstopdat <- with(finalstops, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel, country, laststop))

finalstopdatSW <- subset(finalstopdat, strat=="SW")
finalstopdatSW <- droplevels(finalstopdatSW)

finalstopdatSE <- subset(finalstopdat, strat=="SE")
finalstopdatSE <- droplevels(finalstopdatSE)

finalstopmodelSW <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=finalstopdatSW, family=binomial, control=glmerControl(optimizer="bobyqa"))

finalstopmodelSE <- glmer(presence ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=finalstopdatSE, family=binomial, control=glmerControl(optimizer="bobyqa"))

setwd(paste(outputwd, "/GLMM results", sep=""))

sink(paste("SW GLMM models final stopovers Sahara attempts ", randomradius, " km.txt", sep=""))
cat("\n########==========  SW MODELS ==========########\n", sep="\n")
print(finalstopmodelSW)
cat("\n########==========  SW MODEL SUMMARIES ==========########\n", sep="\n")
print(summary(finalstopmodelSW))

sink()

sink(paste("SE GLMM models final stopovers Sahara attempts ", randomradius, " km.txt", sep=""))
cat("\n########==========  SE MODELS ==========########\n", sep="\n")
print(finalstopmodelSE)
cat("\n########==========  SE MODEL SUMMARIES ==========########\n", sep="\n")
print(summary(finalstopmodelSE))


sink()



#### SW

newdat <- data.frame(PA=rep(levels(finalstopdatSW$PA), each=5), habitat=rep(levels(finalstopdatSW$habitat)), elevation=mean(finalstopdatSW$elevation), spei.Mar=mean(finalstopdatSW$spei.Mar), spei.Aug=mean(finalstopdatSW$spei.Aug))

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

pred <- predict(finalstopmodelSW, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(finalstopmodelSW, newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (FINAL stopovers, Sahara attempts, SW) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(
  axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
  axis.text.x = element_text(size=12), 
  axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
  axis.text.y = element_text(size=12, vjust=0.5),
  strip.text.x = element_text(size=14, face="bold"),
  legend.text=element_text(size=14), 
  legend.title=element_text(size=14)
)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("final stopovers Sahara attempts SW interaction plot 50km with error bars.jpg")


#### SE

newdat <- data.frame(PA=rep(levels(finalstopdatSE$PA), each=5), habitat=rep(levels(finalstopdatSE$habitat)), elevation=mean(finalstopdatSE$elevation), spei.Mar=mean(finalstopdatSE$spei.Mar), spei.Aug=mean(finalstopdatSE$spei.Aug))

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

pred <- predict(finalstopmodelSE, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(finalstopmodelSE, newdat))

newdat2 <- data.frame(newdat, pred, pred.CI)

ggplot(newdat2, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr)) + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) + labs(x="Habitat", y="Predicted Probability", title="Predicted probability of presence (FINAL stopovers, Sahara attempts, SE) at mean elevation and mean climate values, 50km", size=4) + scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water")) + theme(
  axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
  axis.text.x = element_text(size=12), 
  axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
  axis.text.y = element_text(size=12, vjust=0.5),
  strip.text.x = element_text(size=14, face="bold"),
  legend.text=element_text(size=14), 
  legend.title=element_text(size=14)
)

setwd(paste(outputwd, "/GLMM results/", sep=""))
ggsave("final stopovers Sahara attempts SE interaction plot 50km with error bars.jpg")


########################################################################
########################################################################
########################################################################


####===========  SAHARA CROSSING ==============####

#### ---- MODEL DATASET ---- ####

### Convert Sahara success yes/no into integer 1/0

#saharasuc <- subset(alldata.LOS, presence==1 & laststop=="Y")
saharasuc <- subset(alldata.LOS, presence==1)


saharasuc$Sahara.success <- as.numeric(as.character(revalue(saharasuc$Sahara.success, c("Y"=1, "N"=0))))

## add a variable for whether Sahara crossing was attempted or not (ie. remove birds that failed before attempting)
setwd(datawd)
othervar <- read.csv("cuckoo migratory strategy and Sahara crossing success.csv", header=TRUE)

attemptSahara <- with(othervar, data.frame(name, attemptSahara=attempted.Sahara.crossing))

attemptdat <- merge(saharasuc, attemptSahara)

attemptdat <- subset(attemptdat, attemptSahara=="Y")
attemptdat <- droplevels(attemptdat)

saharadat <- with(attemptdat, data.frame(success=Sahara.success, name, year, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=elevation.m, spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel, country, laststop, attemptSahara, LOS))

### Convert point data to observations summarizing the final stopover

### habitat

bybird <- split(saharadat, list(saharadat$name))

bybird2 <- lapply(bybird, function(x) {x$mgroup <- droplevels(x$mgroup); return(x)})

habitatsum <- lapply(bybird2, function(x) prop.table(table(x$mgroup, x$habitat)))

### PA: proportion of points in mgroup overlapping with PA, if no overlap (100% of points are a N), then classify that mgroup as a N for PAoverlap
PAsum <- lapply(bybird2, function(x) prop.table(table(x$mgroup, x$PA), 1))
PA.YN <- lapply(PAsum, function(x) revalue(as.factor(x[,1] == 1), c("TRUE"="N", "FALSE"="Y"))) # identifies as TRUE mgroups which have 100% of points in the "No PA overlap" column

### elevation
meanelevation <- lapply(bybird2, function(x) tapply(x$elevation, x$mgroup, mean))

### climate

mean.spei.Mar <- lapply(bybird2, function(x) tapply(x$spei.Mar, x$mgroup, mean))
mean.spei.Aug <- lapply(bybird2, function(x) tapply(x$spei.Aug, x$mgroup, mean))

### success
success <- lapply(bybird2, function(x) prop.table(table(x$year, x$success), 1))

### LOS
LOS <- lapply(bybird2, function(x) tapply(x$LOS, x$mgroup, mean))


names <- do.call(rbind, lapply(bybird2, function(x) data.frame(name=levels(droplevels(x$name)), mgroup=levels(x$mgroup))))

fulldat <- data.frame(names, do.call(rbind, habitatsum), PA=unlist(PA.YN), elevation=unlist(meanelevation), spei.Mar=unlist(mean.spei.Mar), spei.Aug=unlist(mean.spei.Aug), year=unlist(lapply(success, rownames)), success=unlist(lapply(success, function(x) x[,2])), LOS=unlist(LOS))

m1 <- glmer(success ~ rescale(elevation) + (1|name), data=fulldat, family=binomial)


# saharamodels <- list()
# 
# saharamodels[[1]] <- glmer(success ~ PA*habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=saharadat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# saharamodels[[2]] <- glmer(success ~ PA + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=saharadat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# saharamodels[[3]] <- glmer(success ~ habitat + elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=saharadat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# saharamodels[[4]] <- glmer(success ~ elevation + spei.Mar + spei.Aug + (1|name/mgroup/id), data=saharadat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# 
# tapply(saharasuc$elevation.m, list(as.factor(saharasuc$Sahara.success)), mean)
# 
# newdat <- data.frame(elevation=seq(min(saharadat$elevation), max(saharadat$elevation), by=0.001), spei.Mar=mean(saharadat$spei.Mar), spei.Aug=mean(saharadat$spei.Aug))
# 
# ### Ben Bolker's population-level GLMM CI prediction function
# easyPredCI <- function(model,newdata,alpha=0.05) {
#   ## baseline prediction, on the linear predictor (logit) scale:
#   pred0 <- predict(model,re.form=NA,newdata=newdata)
#   ## fixed-effects model matrix for new data
#   X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
#                     newdata)
#   beta <- fixef(model) ## fixed-effects coefficients
#   V <- vcov(model)     ## variance-covariance matrix of beta
#   pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
#   ## inverse-link (logistic) function: could also use plogis()
#   linkinv <- model@resp$family$linkinv
#   ## construct 95% Normal CIs on the link scale and
#   ##  transform back to the response (probability) scale:
#   crit <- -qnorm(alpha/2)
#   linkinv(cbind(lwr=pred0-crit*pred.se,
#                 upr=pred0+crit*pred.se)) 
# }
# 
# pred <- predict(saharamodels[[4]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(saharamodels[[4]], newdat))
# 
# newdat2 <- data.frame(newdat, pred, pred.CI)
# 
# plot(success ~ elevation, saharadat)
# 
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


####===========  STOPOVER DURATION AS A FUNCTION OF HABITAT + PA ==============####

#### ---- MODEL DATASET ---- ####

stopdat <- subset(alldata.LOS, presence==1)

usedat <- with(stopdat, data.frame(LOS, name, year, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=elevation.m, spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, desig=desiglevel, country, laststop))

### Convert point data to observations summarizing stopovers

### habitat

bybird <- split(usedat, list(usedat$name))

bybird2 <- lapply(bybird, function(x) {x$mgroup <- droplevels(x$mgroup); return(x)})

habitatsum <- lapply(bybird2, function(x) prop.table(table(x$mgroup, x$habitat), 1))

### PA: proportion of points in mgroup overlapping with PA, if no overlap (100% of points are a N), then classify that mgroup as a N for PAoverlap
PAsum <- lapply(bybird2, function(x) prop.table(table(x$mgroup, x$PA), 1))
PA.YN <- lapply(PAsum, function(x) revalue(as.factor(x[,1] == 1), c("TRUE"="N", "FALSE"="Y"))) # identifies as TRUE mgroups which have 100% of points in the "No PA overlap" column

### elevation
meanelevation <- lapply(bybird2, function(x) tapply(x$elevation, x$mgroup, mean))

### climate

mean.spei.Mar <- lapply(bybird2, function(x) tapply(x$spei.Mar, x$mgroup, mean))
mean.spei.Aug <- lapply(bybird2, function(x) tapply(x$spei.Aug, x$mgroup, mean))

### LOS
LOS <- lapply(bybird2, function(x) tapply(x$LOS, x$mgroup, mean))

### year
year <- lapply(bybird2, function(x) tapply(x$year, x$mgroup, mean))

### put all data together
names <- do.call(rbind, lapply(bybird2, function(x) data.frame(name=levels(droplevels(x$name)), mgroup=levels(x$mgroup))))

fulldat <- data.frame(names, do.call(rbind, habitatsum), PA=unlist(PA.YN), elevation=unlist(meanelevation), spei.Mar=unlist(mean.spei.Mar), spei.Aug=unlist(mean.spei.Aug), year=unlist(year), LOS=unlist(LOS))

LOSdat <- data.frame(lapply(fulldat, function(y) if(is.numeric(y)) round(y, 3) else y))

LOSdat2 <- data.frame(LOSdat, sqrtLOS=sqrt(LOSdat$LOS))

###======== COLLINEARITY ========###

### explore collinearity of continuous explanatory variables
exp.cont <- LOSdat2[,c(3:7,9:11)]

#exp.cont <- alldata[,c("agriculture","forest","scrub_grassland","unsuitable","precip.mm.winter","temp.C.winter","precip.mm.spring","temp.C.spring","elevation.m")]

source(paste(parentwd,"/scripts/Zuur_functions.r", sep=""))
# for 6 habitat categories

pairs(exp.cont, lower.panel=panel.cor, panel=panel.smooth)

habpc <- LOSdat2[,3:7]

model <- prcomp(habpc,scale=TRUE)
biplot(model)
summary(model)
model
pc1<-predict(model)[,1] # predicted values of PC1 (agriculture component)
pc2 <- predict(model)[,2] # predicted valeus of PC2 (wetland + unsuitable component)
pc3 <- predict(model)[,3] # forest/scrub component

LOSdat2 <- data.frame(LOSdat2, pc1, pc2, pc3)

LOSmodels <- list()

LOSmodels[[1]] <- lmer(sqrtLOS ~ 1 + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[2]] <- lmer(sqrtLOS ~ PA + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[3]] <- lmer(sqrtLOS ~ elevation + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[4]] <- lmer(sqrtLOS ~ spei.Mar + spei.Aug + (1|name), data=LOSdat2, REML=FALSE)

LOSmodels[[5]] <- lmer(sqrtLOS ~ forest + scrub.grassland + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[6]] <- lmer(sqrtLOS ~ agriculture + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[7]] <- lmer(sqrtLOS ~ wetland.water + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[8]] <- lmer(sqrtLOS ~ forest + scrub.grassland + wetland.water + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[9]] <- lmer(sqrtLOS ~ agriculture + forest + scrub.grassland + wetland.water + unsuitable + (1|name), data=LOSdat2, REML=FALSE)


LOSmodels <- list()

LOSmodels[[1]] <- lmer(sqrtLOS ~ 1 + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[2]] <- lmer(sqrtLOS ~ PA + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[3]] <- lmer(sqrtLOS ~ elevation + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[4]] <- lmer(sqrtLOS ~ spei.Mar + spei.Aug + (1|name), data=LOSdat2, REML=FALSE)

LOSmodels[[5]] <- lmer(sqrtLOS ~ pc1 + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[6]] <- lmer(sqrtLOS ~ pc2 + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[7]] <- lmer(sqrtLOS ~ pc3 + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[8]] <- lmer(sqrtLOS ~ pc1 + pc2 + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[9]] <- lmer(sqrtLOS ~ pc1 + pc3 + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[10]] <- lmer(sqrtLOS ~ pc2 + pc3 + (1|name), data=LOSdat2, REML=FALSE)
LOSmodels[[11]] <- lmer(sqrtLOS ~ pc1 + pc2 + pc3 + (1|name), data=LOSdat2, REML=FALSE)



### AIC function
calculate.AIC<-function(aictable,modellist) {
  modelnames <- modellist
  delta.aic <- aictable$AIC-min(aictable$AIC)
  lik.aic <- exp(-delta.aic/2)
  aic.w <- lik.aic/(sum(lik.aic))
  aic.table <- data.frame(modelnames,aictable,delta.aic,lik.aic,aic.w)
}

m <- LOSmodels

# models <- m[-(which(sapply(m,is.null),arr.ind=TRUE))] # need this line if leaving out models so that list of models has NULL values

# assign models in list m names m1, m2, m3, etc (m1 = m[[1]], m2 = m[[2]])
for (i in 1:length(m)) {
  assign(paste("m", i, sep=""), m[[i]])
}

modelsummary <- lapply(m, summary)
names(modelsummary) <- rep(paste("m", 1:length(m), sep=""))
AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)
                 #,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22)
modellist <- lapply(m, formula)
output <- calculate.AIC(AICoutput, as.character(modellist))
aic.ordered<-output[rev(order(output$aic.w)),]
output <- aic.ordered
output[,3:6] <- round(output[,c(3:6)], digits = 3)
output

plot(sqrt(LOS) ~ wetland.water, data=LOSdat2)
