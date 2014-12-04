##########################################################
#
#   CUCKOO ANALYSES: SUCCESS/FAILURE DIFFERENCES
#
#  Samantha Franks
#  7 Oct 2014 - base code taken from cuckoo_analyses_habitatuse_models_20140915
#
##########################################################

Mac <- FALSE

if(.Platform$OS =='windows') cluster <- FALSE
if(.Platform$OS=='unix' & !Mac) cluster <- TRUE

# if (cluster) {
#   Args <- commandArgs()[3]
#   
#   route <- substr(Args,1,2)  # SW/SE
#   randomradius <- as.numeric(substr(Args,3,5)) # 50, 100, 200, or 500km
#   analysistype <- substr(Args,6,7)
# }
# 
# if (!cluster) {
#   route <- "DE"
#   randomradius <- 200
#   analysistype <- "DE"
# }

if (cluster) {
  library(plyr)
  library(reshape)
  #library(glmmML)
  library(lme4)
  library(arm)
}

if (!cluster) {
  library(plyr)
  library(reshape)
  #library(glmmML)
  library(lme4)
  library(arm)
  library(ggplot2)
  library(scales)
  library(sp)
  library(rgeos)
  library(rgdal)
  library(MuMIn)
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
outputwd <- paste(parentwd, "/output/autumn LOS and Sahara success analysis", sep="")
workspacewd <- paste(parentwd, "/workspaces", sep="")


####==== IMPORT POINT DATA - only need presences for Sahara success ====####

setwd(paste(datawd, "/data for analysis/", sep=""))
present <- read.csv("presence data all variables points all with country data.csv", header=T)

# revalue Sahara success for John (incorrect in main data sheet) => "Y" should be "N"
present[present$name=="John", "Sahara.success"] <- "N"

### ----------------------------- ###
####        DATA PREPARATION     ####
### ------------------------------###

#### ---- CORRECT COLUMN NAMES OF DATASET SO CAN MERGE ---- ####

newpresent <- rename(present, c("newlongs.epsg3035"="long.epsg3035", "newlats.epsg3035"="lat.epsg3035", "newlongs.epsg4326"="long.epsg4326", "newlats.epsg4326"="lat.epsg4326"))

#### ---- ADD STOPOVER DURATION TO DATA ---- ####

### creates 2 files: present.LOS and absent.LOS with length of stay data added onto present and absent datasets
setwd(paste(parentwd, "/scripts/", sep=""))
source("source code to add stopover duration data for analysis.R")

present.LOS <- present.LOS[,-which(names(present.LOS) %in% c("datetime"))]
present.LOS <- present.LOS[c(2,1,3:length(present.LOS))]

present.LOS$mgroup <- as.factor(present.LOS$mgroup)

present.LOS <- rename(present.LOS, c("LAND.CLASS"="habitat", "PAoverlap"="PA", "elevation.m"="elevation", "Sahara.success"="success"))

########################################################################
########################################################################
########################################################################

###########################################################################
#
####    LOGISTIC REGRESSION MIXED MODEL ANALYSIS & MODEL SELECTION     ####
#
###########################################################################


####========================   SAHARA SUCCESS ==========================####


#### ---- CONVERT POINT DATA TO STOPOVER DATA FOR **FINAL STOPOVERS ONLY** ---- ####

### remove birds that had unsuccessful crossing due to tag failure: Wallace
# only want to include birds that are presumed dead in Europe or on crossing or which likely died in Europe or on crossing (gleaned from text of BTO cuckoos page)

finalstops <- subset(present.LOS, laststop=="Y" & name!="Wallace")
finalstops <- droplevels(finalstops)

### habitat

bybird <- split(finalstops, list(finalstops$name))

bybird2 <- lapply(bybird, function(x) {x$mgroup <- droplevels(x$mgroup); x$name <- droplevels(x$name); return(x)})

habitatsum <- lapply(bybird2, function(x) prop.table(table(x$mgroup, x$habitat), 1))

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

# names <- do.call(rbind, lapply(bybird2, function(x) unique(x[c("name","mgroup","year")])))

names <- do.call(rbind, lapply(bybird2, function(x) data.frame(name=levels(x$name), mgroup=levels(x$mgroup))))

#### ---- MODEL DATASET ---- ####

dat <- data.frame(names, do.call(rbind, habitatsum), PA=unlist(PA.YN), elevation=rescale(unlist(meanelevation)), spei.Mar=unlist(mean.spei.Mar), spei.Aug=unlist(mean.spei.Aug), year=unlist(lapply(success, rownames)), success=unlist(lapply(success, function(x) x[,2])), LOS=unlist(LOS))

########################################################################
########################################################################
########################################################################

### ---------------------------------- ###
####        SAHARA SUCCESS MODELS     ####
### --------------------------------- -###

m <- list()

m[[1]] <- glmer(success ~ 1 + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

m[[2]] <- glmer(success ~ PA + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[3]] <- glmer(success ~ elevation + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[4]] <- glmer(success ~ spei.Mar + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[5]] <- glmer(success ~ spei.Aug + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

m[[6]] <- glmer(success ~ agriculture + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[7]] <- glmer(success ~ forest + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[8]] <- glmer(success ~ scrub.grassland + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[9]] <- glmer(success ~ unsuitable + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[10]] <- glmer(success ~ wetland.water + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

m[[11]] <- glmer(success ~ forest + scrub.grassland + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# drop model, max|grad| is too high (11.57)
# m[[12]] <- glmer(success ~ forest + scrub.grassland + wetland.water + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

m[[12]] <- glmer(success ~ PA + agriculture + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[13]] <- glmer(success ~ PA + forest + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[14]] <- glmer(success ~ PA + scrub.grassland + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[15]] <- glmer(success ~ PA + unsuitable + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[16]] <- glmer(success ~ PA + wetland.water + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

m[[17]] <- glmer(success ~ PA + forest + scrub.grassland + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[18]] <- glmer(success ~ PA + forest + scrub.grassland + wetland.water + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

m[[19]] <- glmer(success ~ agriculture*PA + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[20]] <- glmer(success ~ forest*PA + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
m[[21]] <- glmer(success ~ scrub.grassland*PA + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

m[[22]] <- glmer(success ~ PA + agriculture + forest + scrub.grassland + wetland.water + elevation + spei.Mar + spei.Aug + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))



# drop model max|grad| is too high (84.9, tol=0.001)
# m[[14]] <- glmer(success ~ unsuitable*PA + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

# drop model, rank deficient
# m[[15]] <- glmer(success ~ wetland.water*PA + (1|name), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))


###======== OUTPUT RESULTS ========###

output.mumin <- model.sel(m)
modavg.mumin <- model.avg(output.mumin)

setwd(outputwd)
#write.table(output.mumin, file=modseltabfile, row.names=FALSE, sep=",")

modseltabfile <- file("sahara success analysis model selection table.txt", "w")

cat("TERM CODES\n", file=modseltabfile, append=TRUE)
write.table(modavg.mumin$term.codes, file=modseltabfile, na="", sep=",")
cat("\n", file=modseltabfile, append=TRUE)

cat("MODEL SELECTION TABLE\n", file=modseltabfile, append=TRUE)
write.table(data.frame(models=row.names(modavg.mumin$summary), modavg.mumin$summary), file=modseltabfile, na="", row.names=FALSE, sep=",")
cat("\n", file=modseltabfile, append=TRUE)

cat("MODEL-AVERAGED PARAMETER ESTIMATES\n", file=modseltabfile, append=TRUE)
write.table(data.frame(Parameter=row.names(modavg.mumin$avg.model), modavg.mumin$avg.model), file=modseltabfile, row.names=FALSE, sep=",")
cat("\n", file=modseltabfile, append=TRUE)

cat("95% CONFIDENCE INTERVALS\n", file=modseltabfile, append=TRUE)
write.table(data.frame(Parameter=row.names(confint(modavg.mumin)), confint(modavg.mumin, level=0.95)), file=modseltabfile, row.names=FALSE, sep=",")
cat("\n", file=modseltabfile, append=TRUE)

cat("84% CONFIDENCE INTERVALS\n", file=modseltabfile, append=TRUE)
write.table(data.frame(Parameter=row.names(confint(modavg.mumin)), confint(modavg.mumin, level=0.84)), file=modseltabfile, row.names=FALSE, sep=",")
cat("\n\n", file=modseltabfile, append=TRUE)


close(modseltabfile)



easyPredCI <- function(model,newdata,alpha=0.16) {
  
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

# # plot the effects of scrub, ag, spei.Mar, PA, (elevation)
# p <- modavg.mumin$avg.model[,"Estimate"]
# 
# med.elev <- median(dat$elevation)
# med.scr <- median(dat$scrub.grassland)
# med.ag <- median(dat$agriculture)
# med.for <- median(dat$forest)
# med.spM <- median(dat$spei.Mar)
# med.spA <- median(dat$spei.Aug)
# med.uns <- median(dat$unsuitable)
# 
# median.dN <- median(dat$dN)
# int <- p[1]
# N.max <- max(dat$dN)
# N.min <- min(dat$dN)
# D.max <- max(dat$Distance)
# D.min <- min(dat$Distance)
# C.max <- max(dat$dC)
# C.min <- min(dat$dC)
# 
# #parameters p[5:13] are YearSite parameters
# #int #dist #dN #dC #site(change to appropriate parameter depending on site) #dist:dN #dist:dC #dN:dC #dist:dN:dC
# line.dN <- function(x)(int + (p[2]*median.dist) + (p[3]*x) + (p[4]*median.dC) + (p[5]*median.dist*x) + (p[6]*median.dist*median.dC) + (p[7]*x*median.dC) + (p[8]*median.dist*x*median.dC))
# 
# line.D <- function(x)(int + (p[2]*x) + (p[3]*median.dN) + (p[4]*median.dC) + (p[5]*x*median.dN) + (p[6]*x*median.dC) + (p[7]*median.dN*median.dC) + (p[8]*x*median.dN*median.dC))
# 
# line.dC <- function(x)(int + (p[2]*median.dist) + (p[3]*median.dN) + (p[4]*x) + (p[5]*median.dist*median.dN) + (p[6]*median.dist*x) + (p[7]*median.dN*x) + (p[8]*median.dist*median.dN*x))
# 
# plot(Response~Distance, data=dat, type="p", pch=16, ylab = "Standardized nest initiation date", xlab = "Migration distance (km)", las=1, cex.axis=1.2, cex.lab=1.2)
# lines(seq(D.min,D.max,100),line.D(seq(D.min,D.max,100)),col="black",lty=3,lwd=2)
# mtext(c("d) Females"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)
# 
# plot(Response~dC, data=dat, type="p", pch=16, ylab="", xlab = (expression(bold(paste(delta^13,"C"," (???)")))), las=1, cex.axis=1.2, cex.lab=1.2)
# lines(seq(C.min,C.max,0.1),line.dC(seq(C.min,C.max,0.1)),col="black",lty=3,lwd=2)
# mtext(c("e)"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)
# 
# plot(Response~dN, data=dat, type="p", pch=16, ylab="", xlab = (expression(bold(paste(delta^15,"N"," (???)")))), las=1, cex.axis=1.2, cex.lab=1.2)
# lines(seq(N.min,N.max,0.1),line.dN(seq(N.min,N.max,0.1)),col="black",lty=3,lwd=2)
# mtext(c("f)"),side=3,line=0.5,font=2,cex=0.8,adj=0,las=1)


# ### DATA for each habitat x PA combination
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(forest=habvar, PA=rep(levels(dat$PA), each=length(habvar)), agriculture=mean(dat$agriculture), scrub.grassland=mean(dat$scrub.grassland), wetland.water=mean(dat$wetland.water), elevation=mean(dat$elevation), spei.Mar=mean(dat$spei.Mar), spei.Aug=mean(dat$spei.Aug))
# 
# pred <- predict(m[[1]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[1]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=forest, y=pred, data=plotdat, colour=PA)


#######################################
# 
# # agriculture
# 
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(agriculture=habvar)
# 
# pred <- predict(m[[3]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[3]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=agriculture, y=pred, data=plotdat)
# 
# # forest
# 
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(forest=habvar)
# 
# pred <- predict(m[[4]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[4]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=forest, y=pred, data=plotdat)
# 
# # scrub
# 
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(scrub.grassland=habvar)
# 
# pred <- predict(m[[5]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[5]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=scrub.grassland, y=pred, data=plotdat)
# 
# # unsuitable
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(unsuitable=habvar)
# 
# pred <- predict(m[[6]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[6]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=unsuitable, y=pred, data=plotdat)
# 
# # wetland-water
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(wetland.water=habvar)
# 
# pred <- predict(m[[7]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[7]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=wetland.water, y=pred, data=plotdat)
# 
# # elevation
# habvar <- seq(min(dat$elevation), max(dat$elevation), 0.01)
# newdat <- data.frame(elevation=habvar)
# 
# pred <- predict(m[[8]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[8]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=elevation, y=pred, data=plotdat)
# 
# # spei Mar
# habvar <- seq(min(dat$spei.Mar), max(dat$spei.Mar), 0.01)
# newdat <- data.frame(spei.Mar=habvar)
# 
# pred <- predict(m[[9]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[9]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=spei.Mar, y=pred, data=plotdat)
# 
# # spei Aug
# habvar <- seq(min(dat$spei.Aug), max(dat$spei.Aug), 0.01)
# newdat <- data.frame(spei.Aug=habvar)
# 
# pred <- predict(m[[10]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[10]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=spei.Aug, y=pred, data=plotdat)
# 
# # PA + agriculture
# 
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(agriculture=habvar, PA=rep(levels(dat$PA), each=length(habvar)))
# 
# pred <- predict(m[[11]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[11]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=agriculture, y=pred, data=plotdat, colour=PA)
# 
# # PA + forest
# 
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(forest=habvar, PA=rep(levels(dat$PA), each=length(habvar)))
# 
# pred <- predict(m[[12]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[12]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=forest, y=pred, data=plotdat, colour=PA)
# 
# # PA + scrub
# 
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(scrub.grassland=habvar, PA=rep(levels(dat$PA), each=length(habvar)))
# 
# pred <- predict(m[[13]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[13]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=scrub.grassland, y=pred, data=plotdat, colour=PA)
# 
# # PA + unsuitable
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(unsuitable=habvar, PA=rep(levels(dat$PA), each=length(habvar)))
# 
# pred <- predict(m[[14]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[14]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=unsuitable, y=pred, data=plotdat, colour=PA)
# 
# # PA + wetland.water
# 
# habvar <- seq(0,1, 0.01)
# newdat <- data.frame(wetland.water=habvar, PA=rep(levels(dat$PA), each=length(habvar)))
# 
# pred <- predict(m[[15]], newdata=newdat, type="response", re.form=NA)
# pred.CI <- data.frame(easyPredCI(m[[15]], newdat))
# 
# plotdat <- data.frame(newdat, pred, pred.CI)
# 
# qplot(x=wetland.water, y=pred, data=plotdat, colour=PA)

########################################################################
########################################################################
########################################################################

####========================   STOPOVER DURATION    ==========================####

### ----------------------------- ###
####        DATA PREPARATION     ####
### ------------------------------###


#### ---- CONVERT POINT DATA TO STOPOVER DATA FOR **ALL STOPOVERS** ---- ####

### habitat

bybird <- split(present.LOS, list(present.LOS$name))

bybird2 <- lapply(bybird, function(x) {x$mgroup <- droplevels(x$mgroup); x$name <- droplevels(x$name); return(x)})

habitatsum <- lapply(bybird2, function(x) prop.table(table(x$mgroup, x$habitat), 1))

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

# names <- do.call(rbind, lapply(bybird2, function(x) unique(x[c("name","mgroup","year")])))

names <- do.call(rbind, lapply(bybird2, function(x) data.frame(name=levels(x$name), mgroup=levels(x$mgroup))))

#### ---- MODEL DATASET ---- ####

fulldat <- data.frame(names, do.call(rbind, habitatsum), PA=unlist(PA.YN), elevation=rescale(unlist(meanelevation)), spei.Mar=unlist(mean.spei.Mar), spei.Aug=unlist(mean.spei.Aug), LOS=unlist(LOS))

LOSdat <- data.frame(lapply(fulldat, function(y) if(is.numeric(y)) round(y, 3) else y))

LOSdat2 <- data.frame(LOSdat, sqrtLOS=sqrt(LOSdat$LOS))



###################################################################
###################################################################
###################################################################

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

LOSdat3 <- data.frame(LOSdat2, pc1, pc2, pc3)

# exclude Reacher's stopover (he died and tag may have continued transmitting for some time following death without death being detected)
# exclude Wallace's final stopover (mgroup 13) - tag failed during this stop, so stopover may have been longer than detected by tag

LOSdat4 <- LOSdat3[-which(LOSdat3$name=="Wallace" & LOSdat3$mgroup=="13"),]
LOSdat5 <- subset(LOSdat4, name!="Reacher")



###======== MODELS - using true habitat variables, not pc's ========###

output.mumin <- list()
modavg.mumin <- list()
formulas <- list()

dat <- LOSdat5

m <- list()

m[[1]] <- lmer(sqrtLOS ~ 1 + (1|name), data=dat, REML=FALSE)
m[[2]] <- lmer(sqrtLOS ~ PA + (1|name), data=dat, REML=FALSE)
m[[3]] <- lmer(sqrtLOS ~ elevation + (1|name), data=dat, REML=FALSE)
m[[4]] <- lmer(sqrtLOS ~ spei.Mar + spei.Aug + (1|name), data=dat, REML=FALSE)

m[[5]] <- lmer(sqrtLOS ~ forest + scrub.grassland + (1|name), data=dat, REML=FALSE)
m[[6]] <- lmer(sqrtLOS ~ agriculture + (1|name), data=dat, REML=FALSE)
m[[7]] <- lmer(sqrtLOS ~ wetland.water + (1|name), data=dat, REML=FALSE)
m[[8]] <- lmer(sqrtLOS ~ forest + scrub.grassland + wetland.water + (1|name), data=dat, REML=FALSE)
m[[9]] <- lmer(sqrtLOS ~ agriculture + forest + scrub.grassland + wetland.water + unsuitable + (1|name), data=dat, REML=FALSE)

m[[10]] <- lmer(sqrtLOS ~ PA + forest + scrub.grassland + (1|name), data=dat, REML=FALSE)
m[[11]] <- lmer(sqrtLOS ~ PA + agriculture + (1|name), data=dat, REML=FALSE)
m[[12]] <- lmer(sqrtLOS ~ PA + wetland.water + (1|name), data=dat, REML=FALSE)
m[[13]] <- lmer(sqrtLOS ~ PA + forest + scrub.grassland + wetland.water + (1|name), data=dat, REML=FALSE)
m[[14]] <- lmer(sqrtLOS ~ PA + agriculture + forest + scrub.grassland + wetland.water + unsuitable + (1|name), data=dat, REML=FALSE)

output.mumin[[1]] <- model.sel(m)
modavg.mumin[[1]] <- model.avg(output.mumin[[1]])
formulas[[1]] <- lapply(m, formula)

##################################################################
##################################################################
##################################################################

m <- list()

m[[1]] <- lmer(sqrtLOS ~ 1 + (1|name), data=dat, REML=FALSE)
m[[2]] <- lmer(sqrtLOS ~ PA + (1|name), data=dat, REML=FALSE)
m[[3]] <- lmer(sqrtLOS ~ elevation + (1|name), data=dat, REML=FALSE)
m[[4]] <- lmer(sqrtLOS ~ spei.Mar + spei.Aug + (1|name), data=dat, REML=FALSE)

m[[5]] <- lmer(sqrtLOS ~ forest + scrub.grassland + (1|name), data=dat, REML=FALSE)
m[[6]] <- lmer(sqrtLOS ~ agriculture + (1|name), data=dat, REML=FALSE)
m[[7]] <- lmer(sqrtLOS ~ wetland.water + (1|name), data=dat, REML=FALSE)
m[[8]] <- lmer(sqrtLOS ~ forest + scrub.grassland + wetland.water + (1|name), data=dat, REML=FALSE)

m[[9]] <- lmer(sqrtLOS ~ PA + forest + scrub.grassland + (1|name), data=dat, REML=FALSE)
m[[10]] <- lmer(sqrtLOS ~ PA + agriculture + (1|name), data=dat, REML=FALSE)
m[[11]] <- lmer(sqrtLOS ~ PA + wetland.water + (1|name), data=dat, REML=FALSE)
m[[12]] <- lmer(sqrtLOS ~ PA + forest + scrub.grassland + wetland.water + (1|name), data=dat, REML=FALSE)

output.mumin[[2]] <- model.sel(m)
modavg.mumin[[2]] <- model.avg(output.mumin[[2]])
formulas[[2]] <- lapply(m, formula)

##################################################################
##################################################################
##################################################################

###======== MODELS - same as set 2, but including additive climate/elevation to habitat ========###


m <- list()

m[[1]] <- lmer(sqrtLOS ~ 1 + (1|name), data=dat, REML=FALSE)
m[[2]] <- lmer(sqrtLOS ~ PA + (1|name), data=dat, REML=FALSE)
m[[3]] <- lmer(sqrtLOS ~ elevation + (1|name), data=dat, REML=FALSE)
m[[4]] <- lmer(sqrtLOS ~ spei.Mar + spei.Aug + (1|name), data=dat, REML=FALSE)

m[[5]] <- lmer(sqrtLOS ~ forest + scrub.grassland + (1|name), data=dat, REML=FALSE)
m[[6]] <- lmer(sqrtLOS ~ agriculture + (1|name), data=dat, REML=FALSE)
m[[7]] <- lmer(sqrtLOS ~ wetland.water + (1|name), data=dat, REML=FALSE)
m[[8]] <- lmer(sqrtLOS ~ forest + scrub.grassland + wetland.water + (1|name), data=dat, REML=FALSE)

m[[9]] <- lmer(sqrtLOS ~ PA + forest + scrub.grassland + (1|name), data=dat, REML=FALSE)
m[[10]] <- lmer(sqrtLOS ~ PA + agriculture + (1|name), data=dat, REML=FALSE)
m[[11]] <- lmer(sqrtLOS ~ PA + wetland.water + (1|name), data=dat, REML=FALSE)
m[[12]] <- lmer(sqrtLOS ~ PA + forest + scrub.grassland + wetland.water + (1|name), data=dat, REML=FALSE)

m[[13]] <- lmer(sqrtLOS ~ forest + scrub.grassland + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[14]] <- lmer(sqrtLOS ~ agriculture + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[15]] <- lmer(sqrtLOS ~ wetland.water + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[16]] <- lmer(sqrtLOS ~ forest + scrub.grassland + wetland.water + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)

m[[17]] <- lmer(sqrtLOS ~ PA + forest + scrub.grassland + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[18]] <- lmer(sqrtLOS ~ PA + agriculture + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[19]] <- lmer(sqrtLOS ~ PA + wetland.water + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[20]] <- lmer(sqrtLOS ~ PA + forest + scrub.grassland + wetland.water + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)

output.mumin[[4]] <- model.sel(m)
modavg.mumin[[4]] <- model.avg(output.mumin[[4]])
formulas[[4]] <- lapply(m, formula)


##################################################################
##################################################################
##################################################################


###======== MODELS - using PCs ========###

# 2 major outlier points in pc2 (the wetland/unsuitable component) - remove them
# one with large proportion of wetland (0.80 David 35); one with large proportion of unsuitable (0.68 Kasper 2)
LOSdat6 <- subset(LOSdat5, pc2 < 4)
LOSdat6 <- droplevels(LOSdat6)

dat <- with(LOSdat6, data.frame(sqrtLOS, name, pc1, pc2, pc3, PA, elevation, spei.Mar, spei.Aug))


m <- list()

m[[1]] <- lmer(sqrtLOS ~ 1 + (1|name), data=dat, REML=FALSE)
m[[2]] <- lmer(sqrtLOS ~ PA + (1|name), data=dat, REML=FALSE)
m[[3]] <- lmer(sqrtLOS ~ elevation + (1|name), data=dat, REML=FALSE)
m[[4]] <- lmer(sqrtLOS ~ spei.Mar + spei.Aug + (1|name), data=dat, REML=FALSE)

m[[5]] <- lmer(sqrtLOS ~ pc1 + (1|name), data=dat, REML=FALSE)
m[[6]] <- lmer(sqrtLOS ~ pc2 + (1|name), data=dat, REML=FALSE)
m[[7]] <- lmer(sqrtLOS ~ pc3 + (1|name), data=dat, REML=FALSE)
m[[8]] <- lmer(sqrtLOS ~ pc1 + pc2 + (1|name), data=dat, REML=FALSE)
m[[9]] <- lmer(sqrtLOS ~ pc1 + pc3 + (1|name), data=dat, REML=FALSE)
m[[10]] <- lmer(sqrtLOS ~ pc2 + pc3 + (1|name), data=dat, REML=FALSE)
m[[11]] <- lmer(sqrtLOS ~ pc1 + pc2 + pc3 + (1|name), data=dat, REML=FALSE)

output.mumin[[3]] <- model.sel(m)
modavg.mumin[[3]] <- model.avg(output.mumin[[3]])
formulas[[3]] <- lapply(m, formula)


##################################################################
##################################################################
##################################################################

###======== MODELS - using PCs as set 3, but including additive climate/elevation to habitat========###


m <- list()

m[[1]] <- lmer(sqrtLOS ~ 1 + (1|name), data=dat, REML=FALSE)
m[[2]] <- lmer(sqrtLOS ~ PA + (1|name), data=dat, REML=FALSE)
m[[3]] <- lmer(sqrtLOS ~ elevation + (1|name), data=dat, REML=FALSE)
m[[4]] <- lmer(sqrtLOS ~ spei.Mar + spei.Aug + (1|name), data=dat, REML=FALSE)

m[[5]] <- lmer(sqrtLOS ~ pc1 + (1|name), data=dat, REML=FALSE)
m[[6]] <- lmer(sqrtLOS ~ pc2 + (1|name), data=dat, REML=FALSE)
m[[7]] <- lmer(sqrtLOS ~ pc3 + (1|name), data=dat, REML=FALSE)
m[[8]] <- lmer(sqrtLOS ~ pc1 + pc2 + (1|name), data=dat, REML=FALSE)
m[[9]] <- lmer(sqrtLOS ~ pc1 + pc3 + (1|name), data=dat, REML=FALSE)
m[[10]] <- lmer(sqrtLOS ~ pc2 + pc3 + (1|name), data=dat, REML=FALSE)
m[[11]] <- lmer(sqrtLOS ~ pc1 + pc2 + pc3 + (1|name), data=dat, REML=FALSE)

m[[12]] <- lmer(sqrtLOS ~ PA + pc1 + (1|name), data=dat, REML=FALSE)
m[[13]] <- lmer(sqrtLOS ~ PA + pc2 + (1|name), data=dat, REML=FALSE)
m[[14]] <- lmer(sqrtLOS ~ PA + pc3 + (1|name), data=dat, REML=FALSE)
m[[15]] <- lmer(sqrtLOS ~ PA + pc1 + pc2 + (1|name), data=dat, REML=FALSE)
m[[16]] <- lmer(sqrtLOS ~ PA + pc1 + pc3 + (1|name), data=dat, REML=FALSE)
m[[17]] <- lmer(sqrtLOS ~ PA + pc2 + pc3 + (1|name), data=dat, REML=FALSE)
m[[18]] <- lmer(sqrtLOS ~ PA + pc1 + pc2 + pc3 + (1|name), data=dat, REML=FALSE)

m[[19]] <- lmer(sqrtLOS ~ PA + pc1 + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[20]] <- lmer(sqrtLOS ~ PA + pc2 + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[21]] <- lmer(sqrtLOS ~ PA + pc3 + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[22]] <- lmer(sqrtLOS ~ PA + pc1 + pc2 + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[23]] <- lmer(sqrtLOS ~ PA + pc1 + pc3 + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[24]] <- lmer(sqrtLOS ~ PA + pc2 + pc3 + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)
m[[25]] <- lmer(sqrtLOS ~ PA + pc1 + pc2 + pc3 + spei.Mar + spei.Aug + elevation + (1|name), data=dat, REML=FALSE)

output.mumin[[5]] <- model.sel(m)
modavg.mumin[[5]] <- model.avg(output.mumin[[5]])
formulas[[5]] <- lapply(m, formula)

##################################################################
##################################################################
##################################################################


###======== OUTPUT RESULTS ========###


setwd(outputwd)
#write.table(output.mumin, file=modseltabfile, row.names=FALSE, sep=",")

modseltabfile <- file("stopover duration analysis model selection table.csv", "w")
#modavgsummary <- file("stopover duration analysis model summary.txt", "w")

#data.frame(models=as.character(formulas[[1]]), modavg.mumin[[1]]$summary)

for (b in 1:length(output.mumin)) {
  cat(paste("MODEL SET ", b, "\n\n", sep=""), file=modseltabfile, append=TRUE)
  cat("TERM CODES\n", file=modseltabfile, append=TRUE)
  write.table(modavg.mumin[[b]]$term.codes, file=modseltabfile, na="", sep=",")
  cat("\n", file=modseltabfile, append=TRUE)
  
  cat("MODEL SELECTION TABLE\n", file=modseltabfile, append=TRUE)
  write.table(data.frame(models=row.names(modavg.mumin[[b]]$summary), modavg.mumin[[b]]$summary), file=modseltabfile, na="", row.names=FALSE, sep=",")
  cat("\n", file=modseltabfile, append=TRUE)
  
  cat("MODEL-AVERAGED PARAMETER ESTIMATES\n", file=modseltabfile, append=TRUE)
  write.table(data.frame(Parameter=row.names(modavg.mumin[[b]]$avg.model), modavg.mumin[[b]]$avg.model), file=modseltabfile, row.names=FALSE, sep=",")
  cat("\n", file=modseltabfile, append=TRUE)
  
  cat("95% CONFIDENCE INTERVALS\n", file=modseltabfile, append=TRUE)
  write.table(data.frame(Parameter=row.names(confint(modavg.mumin[[b]])), confint(modavg.mumin[[b]], level=0.95)), file=modseltabfile, row.names=FALSE, sep=",")
  cat("\n", file=modseltabfile, append=TRUE)
  
  cat("84% CONFIDENCE INTERVALS\n", file=modseltabfile, append=TRUE)
  write.table(data.frame(Parameter=row.names(confint(modavg.mumin[[b]])), confint(modavg.mumin[[b]], level=0.84)), file=modseltabfile, row.names=FALSE, sep=",")
  cat("\n\n", file=modseltabfile, append=TRUE)
}

close(modseltabfile)


##################################################################
##################################################################
##################################################################

# ### AIC function
# calculate.AIC<-function(aictable,modellist) {
#   modelnames <- modellist
#   delta.aic <- aictable$AIC-min(aictable$AIC)
#   lik.aic <- exp(-delta.aic/2)
#   aic.w <- lik.aic/(sum(lik.aic))
#   aic.table <- data.frame(modelnames,aictable,delta.aic,lik.aic,aic.w)
# }
# 
# # assign models in list m names m1, m2, m3, etc (m1 = m[[1]], m2 = m[[2]])
# for (i in 1:length(m)) {
#   assign(paste("m", i, sep=""), m[[i]])
# }
# 
# modelsummary <- lapply(m, summary)
# names(modelsummary) <- rep(paste("m", 1:length(m), sep=""))
# AICoutput <- AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14)
# #AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9)
# #AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)
# #,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22)
# modellist <- lapply(m, formula)
# output <- calculate.AIC(AICoutput, as.character(modellist))
# aic.ordered<-output[rev(order(output$aic.w)),]
# output <- aic.ordered
# output[,3:6] <- round(output[,c(3:6)], digits = 3)
# output

output.mumin <- model.sel(m)
modavg.mumin <- model.avg(output.mumin)
summary(modavg.mumin)


############################################### MODEL PREDICTIONS using AICcmodavg

modavg(m, "pc3", names(modelsummary), conf.level=0.86)

var <- seq(min(dat$pc2), max(dat$pc2), 0.01)
newdat <- data.frame(pc2=var, PA=rep(levels(dat$PA), each=length(var)), pc1=mean(dat$pc1), pc3=mean(dat$pc3), elevation=mean(dat$elevation), spei.Mar=mean(dat$spei.Mar), spei.Aug=mean(dat$spei.Aug))

x <- modavgpred(m, names(modelsummary), newdat, type="response")

plotdat <- data.frame(newdat, x)

plot(mod.avg.pred~pc2, plotdat)

qplot(x=pc2, y=mod.avg.pred, data=plotdat, geom="line")

p0 <- ggplot(plotdat, aes(x=pc2,y=mod.avg.pred, group=1))
p0 + geom_smooth()


modavgpred()

plot(sqrtLOS ~ pc2, data=dat)
