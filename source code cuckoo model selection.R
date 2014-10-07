########################################################################################
#
#  SOURCE CODE TO OUTPUT RESULTS OF CUCKOO MODEL SELECTION
#
#  Samantha Franks
#  21 August 2014
#
########################################################################################
  
### AIC function
calculate.AIC<-function(aictable,modellist) {
  modelnames <- modellist
  delta.aic <- aictable$AIC-min(aictable$AIC)
  lik.aic <- exp(-delta.aic/2)
  aic.w <- lik.aic/(sum(lik.aic))
  aic.table <- data.frame(modelnames,aictable,delta.aic,lik.aic,aic.w)
}

### Capture warnings and errors from models
# myTryCatch <- function(expr) {
#   warn <- err <- NULL
#   value <- withCallingHandlers(
#     tryCatch(expr, error=function(e) {
#       err <<- e
#       NULL
#     }), warning=function(w) {
#       warn <<- w
#       invokeRestart("muffleWarning")
#     })
#   list(value=value, warning=warn, error=err)
# }


# models <- m[-(which(sapply(m,is.null),arr.ind=TRUE))] # need this line if leaving out models so that list of models has NULL values

# assign models in list m names m1, m2, m3, etc (m1 = m[[1]], m2 = m[[2]])
for (i in 1:length(m)) {
  assign(paste("m", i, sep=""), m[[i]])
}

modelsummary <- lapply(m, summary)
names(modelsummary) <- rep(paste("m", 1:length(m), sep=""))
names(m) <- rep(paste("m", 1:length(m), sep=""))

### capture warning messages
# modelsummary.withwarn <- lapply(m, function(x) myTryCatch(summary(x))) # uses tryCatch() function
warningmessages <- lapply(m, function(x) slot(x, "optinfo")$conv$lme4$messages)
                                           
# AICoutput <- AIC(m1,m2)
# AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21)

AICoutput <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28)
modellist <- lapply(m, formula)
output <- calculate.AIC(AICoutput, as.character(modellist))
aic.ordered<-output[rev(order(output$aic.w)),]
output <- aic.ordered
output[,3:6] <- round(output[,c(3:6)], digits = 3)

topmodel <- m[[rownames(output)[1]]]

setwd(outputwd)

write.csv(output, paste(analysistype, route, randomradius, "km models AIC table output.csv", sep=" "))

sink(paste(analysistype, route, randomradius, "km models.txt", sep=" "))
cat("########==========  ", route, " ", analysistype, " ", randomradius, " km MODELS AIC OUTPUT ==========########\n\n", sep="")
print(output)
cat("\n########==========  TOP MODEL PARAMETER TABLE ==========########\n", sep="\n")
print(formula(topmodel))
print(summary(topmodel)$coef)
cat("\n########==========  MODEL SUMMARIES ==========########\n", sep="\n")
print(modelsummary)
cat("\n########==========  WARNING MESSAGES ==========########\n", sep="\n")
print(warningmessages)

sink()

