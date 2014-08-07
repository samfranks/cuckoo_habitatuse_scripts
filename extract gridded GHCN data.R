########################################################
#
# Extract gridded GHCN data
#
########################################################

climatedata <- c("precipitation")

if (climatedata == "precipitation") {
  
  setwd(paste(datawd,"/GHCN grid/", sep=""))
  dd <- read.csv("precipitation anomaly grid 09-2010 to 02-2014.csv", header=T)
  
  by.yearmonth <- dd$yearmonth
  
  newdat <- split(dd, list(by.yearmonth))
  
  precip.yearmonth <- list()
  
  long.step <- seq(-180,175, by=5)
  lat.step <- seq(90, -85, by=-5)
  
  for (a in 1:length(newdat)) {
    
    newdat[[a]] <- droplevels(newdat[[a]])
    
    long <- numeric()
    lat <- numeric()
    precip <- numeric()
    
    usedat <- newdat[[a]][,4:15]    
    
    latbandcount <- 1
    
    for (i in 1:36) {
      # i counts latitude bands
      
      for (j in 1:6) {
        # j counts rows of single latitude band
        
        for (k in 1:12) {
          # k counts longitude columns
          precip <- c(precip,usedat[latbandcount,k])
        }
        
        latbandcount <- latbandcount + 1
        
      }
    }
    
    # create matrix to check that new variable precip creates correct and same dataset as usedat
    # comparison <- matrix(data=precip, nrow=nrow(usedat), ncol=ncol(usedat), byrow=TRUE)
    
    ###--- Create latitude and longitude variables ---###
    
    lat <- rep(lat.step,72)
    lat <- sort(lat, decreasing=TRUE)      
    long <- rep(long.step,36)
    
    ###--- Replace -32768 with NA ---###
    
    precip2 <- precip
    
    precip2[which(precip2 == -32768)] <- NA
    
    precip3 <- precip2/100
    
    ###-- Create new dataframe ---###
    
    precip.yearmonth[[a]] <- data.frame(
      year=rep(newdat[[a]]$year, 12),
      month=rep(newdat[[a]]$month, 12),
      yearmonth=rep(newdat[[a]]$yearmonth, 12),
      lat,
      long,
      precip.mm=precip3
    )
    
    precip.yearmonth[[a]]
    
  }
  
  allprecip <- do.call(rbind,precip.yearmonth)
  write.csv(allprecip, "GHCN grid precip data 201009-201402.csv", row.names=F)
  
}

if (climatedata == "temperature") {
  
  setwd(paste(datawd,"/GHCN grid/", sep=""))
  dd <- read.csv("temperature anomaly grid 09-2010 to 02-2014.csv", header=T)
  
  by.yearmonth <- dd$yearmonth
  
  newdat <- split(dd, list(by.yearmonth))
  
  temp.yearmonth <- list()
  
  long.step <- seq(-180,175, by=5)
  lat.step <- seq(90, -85, by=-5)
  
  for (a in 1:length(newdat)) {
    
    newdat[[a]] <- droplevels(newdat[[a]])
    
    long <- numeric()
    lat <- numeric()
    temp <- numeric()
    
    usedat <- newdat[[a]][,4:75]    
    
    latbandcount <- 1
    
    for (i in 1:36) {
      # i counts latitude bands
      
      for (k in 1:72) {
        # k counts longitude columns
        temp <- c(temp,usedat[latbandcount,k])
      }
      
      latbandcount <- latbandcount + 1
      
    }
    

    ###--- Create latitude and longitude variables ---###
    
    lat <- rep(lat.step,72)
    lat <- sort(lat, decreasing=TRUE)      
    long <- rep(long.step,36)
    
    ###--- Replace -9999 with NA ---###
    
    temp2 <- temp
    
    temp2[which(temp2 == -9999)] <- NA
    
    temp3 <- temp2/100
    
    ###-- Create new dataframe ---###
    
    temp.yearmonth[[a]] <- data.frame(
      year=rep(newdat[[a]]$year, 72),
      month=rep(newdat[[a]]$month, 72),
      yearmonth=rep(newdat[[a]]$yearmonth, 72),
      lat,
      long,
      temp.C=temp3
    )
        
  }
  
  alltemp <- do.call(rbind, temp.yearmonth)
  write.csv(alltemp, "GHCN grid temp data 201009-201402.csv", row.names=F)
  
}