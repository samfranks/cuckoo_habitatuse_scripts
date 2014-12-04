##########################################################
#
#   CUCKOO ANALYSES: MODEL OUTPUT - AIC OUTPUT AND GRAPHS OF TOP MODELS
#
#   Samantha Franks
#   16 Sep 2014 (based on initial code from cuckoo_analyses_models_20140805.R)
#
##########################################################

######## NOTES ON ANALYSIS ########

# Model output derived from RData files generated from cuckoo_analyses_habitatuse_models_20140915.R
# Source file to output AIC table and table of parameter estimates

######################
#
####   FUNCTIONS  ####
#
######################

### Ben Bolker's function for calculating CIs on predictions from an merMod object and plotting the results from his RPubs GLMM worked examples
# http://rpubs.com/bbolker/glmmchapter
# by specifying re.form=NA we're saying that we want the population-level prediction, i.e. setting the random effects to zero and getting a prediction for an average (or unknown) group
# Computing confidence intervals on the predicted values is relatively easy if we're willing to completely ignore the random effects, and the uncertainty of the random effects

###########################################  alpha=0.16 is approximately equal to 84% CIs  ######################################


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

### FUNCTION TO LOAD REQUIRED PACKAGES
habitatusepackages <- function() {
  library(plyr)
  library(reshape)
  library(lme4)
  library(arm)
  library(ggplot2)
  library(grid)
  library(scales)
}

### Capture warnings and errors from models
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

#########################################################################################
#########################################################################################
#########################################################################################

# loads habitatusepackages (as a function so it can be saved for use later)
habitatusepackages()

Mac <- FALSE

###----------------------------------------------###
####         SET WORKING DIRECTORIES		####
###----------------------------------------------###

if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output/autumn stopover GLMM results", sep="")
workspacewd <- paste(parentwd, "/workspaces/rerun workspaces", sep="")

setwd(workspacewd)
n <- length(list.files())

for (a in 1:n) {
  
  setwd(workspacewd)
  load(list.files()[a])
  
####################################################
#
####   GENERATE AIC OUTPUT & MODEL SUMMARIES	####
#
####################################################

# m <- SEmodels
# route <- "SE"
# randomradius <- 200
# analysistype <- "FM"

setwd(paste(parentwd, "/scripts/", sep=""))
source("source code cuckoo model selection.R")

### CLEAN WORKSPACE OF LARGE OBJECTS, leave packages and certain variables
# loaded model objects are very large and slow R down
# clear workspace of all objects EXCEPT: libraries, Mac, parentwd, datawd, outputwd, workspacewd, n, and a
# # rm(list=(ls()[-which(ls()=="Mac" | ls()=="parentwd" | ls()=="datawd" | ls()=="outputwd" | ls()=="workspacewd" | ls()=="n" | ls()=="a" | ls()=="analysistype" | ls()=="dat" | ls()=="route" | ls()=="randomradius"]))

keepvars <- c("Mac", "parentwd", "datawd", "outputwd", "workspacewd", "n", "a", "analysistype", "dat", "route", "randomradius", "topmodel")

# removes everything but functions and variables to keep; function habitatusepackages will reload required packages
allvars.nofunctions <- setdiff(ls(), lsf.str())
rm(list=setdiff(allvars.nofunctions, keepvars))

######################
#
####   PLOTS  ####
#
######################

####==== HABITAT & PROTECTED AREA INTERACTION PLOT ====####

### DATA
newdat <- data.frame(habitat=levels(dat$habitat), PA=rep(levels(dat$PA), 5), elevation=mean(dat$elevation), spei.Mar=mean(dat$spei.Mar), spei.Aug=mean(dat$spei.Aug))

pred <- predict(topmodel, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(topmodel, newdat))

plotdat <- data.frame(newdat, pred, pred.CI)

### PLOT SPECIFICATIONS
plottheme <- theme(
  axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
  axis.text.x = element_text(size=14, colour="black"), 
  axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
  axis.text.y = element_text(size=14, vjust=0.5, colour="black"),
  plot.title = element_text(size=16, face="bold", vjust=0.5),
  legend.text=element_text(size=14), 
  legend.title=element_text(size=14),
  legend.key=element_rect(fill="white"),
  legend.key.width=unit(2,"cm"),
  panel.background=element_rect(colour="black", fill="white"),
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank()
)

if (analysistype=="FM") analysislabel <- "OVERALL"
if (analysistype=="DL") analysislabel <- "DESIGNATION LEVEL"
if (analysistype=="SS") analysislabel <- "SHORT STOPOVERS"
if (analysistype=="LS") analysislabel <- "LONG STOPOVERS"
if (analysistype=="FS") analysislabel <- "FINAL STOPOVERS"
if (analysistype=="DW") analysislabel <- "DESIGNATION LEVEL, LUMPING No and National for SW"
if (analysistype=="FR") analysislabel <- "FRANCE"
if (analysistype=="DE") analysislabel <- "GERMANY & LOW COUNTRIES"
if (analysistype=="IT") analysislabel <- "ITALY"
if (analysistype=="SP") analysislabel <- "SPAIN"

plotlabels <- labs(x="Habitat", y="Predicted Probability", title=paste("Predicted probability of presence at mean elevation and climate values", route, randomradius, "km", analysislabel, sep=" "))

xscale <- scale_x_discrete("Habitat", labels = c("scrub.grassland"="scrub/grassland", "wetland.water"="wetland/water"))

#habitatPA.yscale <- scale_y_continuous(limits=c(0,max(plotdat$upr) + 0.2))

### PLOT
p0 <- ggplot(plotdat, aes(x=habitat,y=pred, shape=PA, linetype=PA, group=PA, ymin=lwr, ymax=upr))
p0 + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) +  plottheme + plotlabels + xscale


### SAVE PLOT
setwd(outputwd)
ggsave(paste(analysistype, route, randomradius, "km top model habitat PA interaction plot.jpg", sep=" "), height=9.4, width=14, units="in")

####==== HABITAT PLOT (holding all else at their means) ====####

plotdat2 <- aggregate(plotdat[,c("pred","lwr","upr")], by=list(habitat=plotdat$habitat, elevation=plotdat$elevation, spei.Mar=plotdat$spei.Mar, spei.Aug=plotdat$spei.Aug), mean)
plotdat2 <- rename(plotdat2, c("x"="pred"))

p1 <- ggplot(plotdat2, aes(x=habitat, y=pred, ymin=lwr, ymax=upr, group=1))
p1 + geom_point(size=4) + geom_line() + geom_errorbar(width=0.1) +  plottheme + plotlabels + xscale

### SAVE PLOT
setwd(outputwd)
ggsave(paste(analysistype, route, randomradius, "km top model mean habitat plot.jpg", sep=" "), height=9.4, width=14, units="in")


####==== ELEVATION PLOT ====####

### DATA for each habitat x PA combination

elevationseq <- seq(min(dat$elevation), max(dat$elevation), 0.01)
newdat <- data.frame(elevation=elevationseq, habitat=rep(levels(dat$habitat), each=length(elevationseq)), PA=rep(levels(dat$PA), each=length(elevationseq)*length(levels(dat$habitat))) , spei.Mar=mean(dat$spei.Mar), spei.Aug=mean(dat$spei.Aug))

pred <- predict(topmodel, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(topmodel, newdat))

plotdat <- data.frame(newdat, pred, pred.CI)

### PLOT SPECIFICATIONS
elevationplotlabels <- labs(x="Elevation", y="Predicted Probability", title=paste("Predicted probability of presence vs elevation", route, randomradius, "km", analysislabel, sep=" "))

elevation.yscale <- scale_y_continuous(limits=c(0,max(plotdat$pred) + 0.1))

### PLOT
p0 <- ggplot(plotdat, aes(x=elevation, y=pred, colour=habitat, linetype=PA, ymin=lwr, ymax=upr))

# without confidence bands
p0 + geom_line(size=1.5) + plottheme + elevationplotlabels + elevation.yscale + scale_color_manual(values=c("#FFCC33", "#009900", "#FF6600", "#FF0000", "#0000FF"))

# # with confidence bands
# p0 + geom_line(size=1.5) + geom_ribbon(aes(fill=habitat, linetype=NA), alpha=0.1) + plottheme + elevationplotlabels +  elevation.yscale _ scale_color_manual(values=c("#FFCC33", "#009900", "#FF6600", "#FF0000", "#0000FF")) + scale_fill_manual(values=c("#FFCC33", "#009900", "#FF6600", "#FF0000", "#0000FF"))

### SAVE PLOT
setwd(outputwd)
ggsave(paste(analysistype, route, randomradius, "km top model elevation habitat + PA plot.jpg", sep=" "), height=9.4, width=14, units="in")

#-----------------------

### DATA for averaging for each habitat type across PA level

plotdat2 <- aggregate(plotdat$pred, by=list(elevation=plotdat$elevation, habitat=plotdat$habitat), mean)
plotdat2 <- rename(plotdat2, c("x"="pred"))

### PLOT
p1 <- ggplot(plotdat2, aes(x=elevation, y=pred, colour=habitat))
p1 + geom_line(size=1.5) + plottheme + elevationplotlabels +  elevation.yscale + scale_color_manual(values=c("#FFCC33", "#009900", "#FF6600", "#FF0000", "#0000FF"))

### SAVE PLOT
setwd(outputwd)
ggsave(paste(analysistype, route, randomradius, "km top model elevation mean habitat plot.jpg", sep=" "), height=9.4, width=14, units="in")


####==== SUMMER CLIMATE PLOT ====####

### DATA for each habitat x PA combination
climateseq <- seq(min(dat$spei.Aug), max(dat$spei.Aug), 0.01)
newdat <- data.frame(spei.Aug=climateseq, habitat=rep(levels(dat$habitat), each=length(climateseq)), PA=rep(levels(dat$PA), each=length(climateseq)*length(levels(dat$habitat))) , elevation=mean(dat$elevation), spei.Mar=mean(dat$spei.Mar))

pred <- predict(topmodel, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(topmodel, newdat))

plotdat <- data.frame(newdat, pred, pred.CI)

### PLOT SPECIFICATIONS
climateplotlabels <- labs(x="Climate", y="Predicted Probability", title=paste("Predicted probability of presence vs growing season climate", route, randomradius, "km", analysislabel, sep=" "))

climate.yscale <- scale_y_continuous(limits=c(0,max(plotdat$pred) + 0.1))

### PLOT
p0 <- ggplot(plotdat, aes(x=spei.Aug, y=pred, colour=habitat, linetype=PA, ymin=lwr, ymax=upr))

# without confidence bands
p0 + geom_line(size=1.5) + plottheme + climateplotlabels + climate.yscale + scale_color_manual(values=c("#FFCC33", "#009900", "#FF6600", "#FF0000", "#0000FF"))

### SAVE PLOT
setwd(outputwd)
ggsave(paste(analysistype, route, randomradius, "km top model spring climate habitat + PA plot.jpg", sep=" "), height=9.4, width=14, units="in")

#-----------------------

### DATA for averaging for each habitat type across PA level

plotdat2 <- aggregate(plotdat$pred, by=list(spei.Aug=plotdat$spei.Aug, habitat=plotdat$habitat), mean)
plotdat2 <- rename(plotdat2, c("x"="pred"))

### PLOT
p1 <- ggplot(plotdat2, aes(x=spei.Aug, y=pred, colour=habitat))
p1 + geom_line(size=1.5) + plottheme + climateplotlabels +  climate.yscale + scale_color_manual(values=c("#FFCC33", "#009900", "#FF6600", "#FF0000", "#0000FF"))

### SAVE PLOT
setwd(outputwd)
ggsave(paste(analysistype, route, randomradius, "km top model spring climate mean habitat plot.jpg", sep=" "), height=9.4, width=14, units="in")


####==== WINTER CLIMATE PLOT ====####

climateseq <- seq(min(dat$spei.Mar), max(dat$spei.Mar), 0.01)
newdat <- data.frame(spei.Mar=climateseq, habitat=rep(levels(dat$habitat), each=length(climateseq)), PA=rep(levels(dat$PA), each=length(climateseq)*length(levels(dat$habitat))) , elevation=mean(dat$elevation), spei.Aug=mean(dat$spei.Aug))

pred <- predict(topmodel, newdata=newdat, type="response", re.form=NA)
pred.CI <- data.frame(easyPredCI(topmodel, newdat))

plotdat <- data.frame(newdat, pred, pred.CI)

### PLOT SPECIFICATIONS
climateplotlabels <- labs(x="Climate", y="Predicted Probability", title=paste("Predicted probability of presence vs winter climate", route, randomradius, "km", analysislabel, sep=" "))

climate.yscale <- scale_y_continuous(limits=c(0,max(plotdat$pred) + 0.1))

### PLOT
p0 <- ggplot(plotdat, aes(x=spei.Mar, y=pred, colour=habitat, linetype=PA, ymin=lwr, ymax=upr))

# without confidence bands
p0 + geom_line(size=1.5) + plottheme + climateplotlabels + climate.yscale + scale_color_manual(values=c("#FFCC33", "#009900", "#FF6600", "#FF0000", "#0000FF"))

### SAVE PLOT
setwd(outputwd)
ggsave(paste(analysistype, route, randomradius, "km top model winter climate habitat + PA plot.jpg", sep=" "), height=9.4, width=14, units="in")

#-----------------------

### DATA for averaging for each habitat type across PA level

plotdat2 <- aggregate(plotdat$pred, by=list(spei.Mar=plotdat$spei.Mar, habitat=plotdat$habitat), mean)
plotdat2 <- rename(plotdat2, c("x"="pred"))

### PLOT
p1 <- ggplot(plotdat2, aes(x=spei.Mar, y=pred, colour=habitat))
p1 + geom_line(size=1.5) + plottheme + climateplotlabels +  climate.yscale + scale_color_manual(values=c("#FFCC33", "#009900", "#FF6600", "#FF0000", "#0000FF"))

### SAVE PLOT
setwd(outputwd)
ggsave(paste(analysistype, route, randomradius, "km top model winter climate mean habitat plot.jpg", sep=" "), height=9.4, width=14, units="in")

########################################################################
########################################################################
########################################################################

}

