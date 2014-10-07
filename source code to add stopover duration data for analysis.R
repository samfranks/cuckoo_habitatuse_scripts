##########################################################
#
# Source code to add stopover duration to original data
#
# Samantha Franks
# 29 August 2014
#
##########################################################

library(chron)

# parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
# 
# datawd <- paste(parentwd, "/data", sep="")

####==== IMPORT POINT DATA ====####

#### Add LOS to PRESENT data

#setwd(paste(datawd, "/data for analysis/", sep=""))

#present <- read.csv("presence data all variables points.csv", header=T)

presentordered <- newpresent[order(newpresent$name, newpresent$location_date),]
x <- do.call(rbind,(strsplit(as.character(presentordered$timestamp), "  ")))

dat.datetime <- chron(x[,1], x[,2], format=c(dates="y/m/d", times="h:m:s"))

present.date <- data.frame(presentordered, datetime=dat.datetime)

newdat <- split(present.date, list(present.date$name)) # split dataset by individual bird

calculate.LOS <- function(fulldata) {
  x <- tapply(fulldata$datetime, fulldata$mgroup, diff)
  LOS <- do.call(rbind, lapply(x, sum))
  LOStable <- data.frame(mgroup=levels(as.factor(fulldata$mgroup)), LOS)
  newdat.LOS <- merge(fulldata, LOStable, by.x="mgroup", by.y="mgroup")
  return(newdat.LOS)
}

makeLOStable <- function(fulldata) {
  x <- tapply(fulldata$datetime, fulldata$mgroup, diff)
  LOS <- do.call(rbind, lapply(x, sum))
  LOStable <- data.frame(name=fulldata$name[1], mgroup=levels(as.factor(fulldata$mgroup)), LOS)
  return(LOStable)
}

present.LOS <- do.call(rbind, lapply(newdat, calculate.LOS))
LOStable.ref <- do.call(rbind, lapply(newdat, makeLOStable))

#### Add LOS to ABSENT data

# absent <- read.csv(paste("absence data all variables points ", 50, " km.csv", sep=""), header=T)
absent.LOS <- merge(newabsent, LOStable.ref, by.x=c("name","mgroup"), by.y=c("name","mgroup"))