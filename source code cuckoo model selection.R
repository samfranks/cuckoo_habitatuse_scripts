########################################################################################
#
#  SOURCE CODE TO OUTPUT RESULTS OF CUCKOO MODEL SELECTION
#
#  Samantha Franks
#  21 August 2014
#
########################################################################################

### change to FALSE to load SE model results, change randomradius to load 200km scale results
# SW <- FALSE
# randomradius <- 50

if (SW) load(paste("~/Git/cuckoos/SW ", randomradius, "km models 20140821.RData", sep=""))
if (!SW) load(paste("~/Git/cuckoos/SE ", randomradius, "km models 20140821.RData", sep=""))

### AIC function
calculate.AIC<-function(aictable,modellist) {
  modelnames <- modellist
  delta.aic <- aictable$AIC-min(aictable$AIC)
  lik.aic <- exp(-delta.aic/2)
  aic.w <- lik.aic/(sum(lik.aic))
  aic.table <- data.frame(modelnames,aictable,delta.aic,lik.aic,aic.w)
}

if (SW) m <- SWmodels
if (!SW) m <- SEmodels

# models <- m[-(which(sapply(m,is.null),arr.ind=TRUE))] # need this line if leaving out models so that list of models has NULL values

# assign models in list m names m1, m2, m3, etc (m1 = m[[1]], m2 = m[[2]])
for (i in 1:length(m)) {
  assign(paste("m", i, sep=""), m[[i]])
}

modelsummary <- lapply(m, summary)
names(modelsummary) <- rep(paste("m", 1:21, sep=""))
AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21)
modellist <- lapply(m, formula)
output <- calculate.AIC(AICoutput, as.character(modellist))
aic.ordered<-output[rev(order(output$aic.w)),]
output <- aic.ordered
output[,3:6] <- round(output[,c(3:6)], digits = 3)

setwd(paste(outputwd, "/GLMM results", sep=""))
if (SW) write.csv(output, paste("SW models AIC output ", randomradius, "km.csv"))
if (!SW) write.csv(output, paste("SE models AIC output ", randomradius, "km.csv"))


if (SW){
  sink(paste("SW GLMM models ", randomradius, " km.txt", sep=""))
  cat("########==========  SW MODELS AIC OUTPUT ==========########\n", sep="\n")
  print(output)
  cat("\n########==========  SW MODELS ==========########\n", sep="\n")
  print(m)
  cat("\n########==========  SW MODEL SUMMARIES ==========########\n", sep="\n")
  print(modelsummary)
}

if (!SW) {
  sink(paste("SE GLMM models ", randomradius, " km.txt", sep=""))
  cat("\n########==========  SE MODELS AIC OUTPUT ==========########\n", sep="\n")
  print(output)
  cat("\n########==========  SE MODELS ==========########\n", sep="\n")
  print(m)
  cat("\n########==========  SE MODEL SUMMARIES ==========########\n", sep="\n")
  print(modelsummary)
}

sink()

#dat$habitat <- factor(dat$habitat, levels=c("unsuitable","agriculture","scrub.grassland","wetland.water","forest"))