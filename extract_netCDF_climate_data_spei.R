######################################################
#
#     EXTRACT ENSEMBLES EUROPEAN CLIMATE DATA - TO CALCULATE SPEI
#
######################################################

# uses ENSEMBLES European climate data downloaded as NETCDF files from http://www.ecad.eu/download/ensembles/download.php
# uses the 15-year chunk datafiles on the 0.25 degree grid for daily mean temperature (TG) and daily precipitation sum (RR) for 1995-2013

#### TO DO ####

# convert 3 dimensional data array into a 2 dimensional file with variables date, lat, long, and value (climvar) - extract only needed time series and lat/long cross-section



require(ncdf)
require(chron)
require(RColorBrewer)
require(reshape)
require(SPEI)

###-----------------------------------------------------------###
#         SET WORKING DIRECTORIES
###-----------------------------------------------------------###

if (!cluster){
  parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
}

if (cluster){
  parentwd <- c("/users1/samf/cuckoos")
}


datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

setwd(paste(datawd, "/E-OBS climate dataset", sep=""))

if (climvarname == "precipitation") {
  filename1 <- "rr_0.25deg_reg_1980-1994_v10.0.nc"
  filename2 <- "rr_0.25deg_reg_1995-2013_v10.0.nc"
}

if (climvarname == "temperature") {
  filename1 <- "tg_0.25deg_reg_1980-1994_v10.0.nc"
  filename2 <- "tg_0.25deg_reg_1995-2013_v10.0.nc"
}

###-----------------------------------------------------------###
#         OPEN NETCDF AND EXTRACT DATA
###-----------------------------------------------------------###

for (a in 1:2) {
  
  if (a==1) {filename <- filename1} else {filename <- filename2}
  
  ncin <- open.ncdf(filename) 
  print(ncin)
  
  # summary function for ncdf objects
  summary.ncdf <-
    function( control){
      for (i in 1:control$nvars ){
        
        vname = control$var[[i]]$name    # variable name
        ndims = control$var[[i]]$ndims   # number of dimensions for this variable
        
        dimstring = paste(vname,'( variable',as.character(i),') has shape')
        for (j in 1:ndims) {
          dimstring <- paste(dimstring, as.character(control$var[[i]]$dim[[j]]$len))
        }
        
        cat(dimstring, fill=TRUE)
      }
    }
  
  summary.ncdf(ncin)
  
  # time units for extracting whole time series for calculating SPI
  t <- get.var.ncdf(ncin, "time") # all t; if specific time series, then use start= and count=; e.g. get.var.ncdf(ncin, "time", start=5723, count=1218)
  tunits <- att.get.ncdf(ncin, "time", "units")
  nt <- dim(t)  # get the number of values and list them
  
  
  lon <- get.var.ncdf(ncin, "longitude", start=90, count=194)
  #lon <- get.var.ncdf(ncin, "longitude") # all lon
  nlon <- dim(lon)
  # 0.25 degree grid lon units -18.125 lon[90] to 30.125 lon[283]
  # need lon units -18.25 lon[45] to 30.25 lon[142] (0.5 degree grid)
  
  lat <- get.var.ncdf(ncin, "latitude", start=39, count=102)
  #lat <- get.var.ncdf(ncin, "latitude") # all lat
  nlat <- dim(lat)
  # 0.25 degree grid lat units 34.875 lat[39] to lat[140]
  # need lat units 34.75 lat[20] to 60.25 lat[71]
  
  print(c(nlon, nlat, nt))  # confirms the dimensions of the data
  
  ###--- EXTRACT CLIMATE VARIABLE FROM NETCDF ---###
  
  if (a==1) { # 1980-1994
    if (climvarname == "precipitation") {
      climvar <- get.var.ncdf(ncin, "rr", start=c(90,39,1), count=c(194,102,5479)) # for start and count, 1st dimension=lat, 2nd dimension=long, 3rd dimension=time
      dlname <- att.get.ncdf(ncin, "rr", "long_name")
      dunits <- att.get.ncdf(ncin, "rr", "units")
      fillvalue <- att.get.ncdf(ncin, "rr", "_FillValue")
      #fillvalue$value <- -999.9
    }
    
    
    if (climvarname == "temperature") {
      climvar <- get.var.ncdf(ncin, "tg", start=c(90,39,1), count=c(194,102,5479))
      dlname <- att.get.ncdf(ncin, "tg", "long_name")
      dunits <- att.get.ncdf(ncin, "tg", "units")
      fillvalue <- att.get.ncdf(ncin, "tg", "_FillValue")
      #fillvalue$value <- -999.9
    } 
  }
  
  
  if (a==2) { # 1995-2013
    if (climvarname == "precipitation") {
      climvar <- get.var.ncdf(ncin, "rr", start=c(90,39,1), count=c(194,102,6940)) # for start and count, 1st dimension=lat, 2nd dimension=long, 3rd dimension=time
      dlname <- att.get.ncdf(ncin, "rr", "long_name")
      dunits <- att.get.ncdf(ncin, "rr", "units")
      fillvalue <- att.get.ncdf(ncin, "rr", "_FillValue")
      #fillvalue$value <- -999.9
    }
    
    
    if (climvarname == "temperature") {
      climvar <- get.var.ncdf(ncin, "tg", start=c(90,39,1), count=c(194,102,6940))
      dlname <- att.get.ncdf(ncin, "tg", "long_name")
      dunits <- att.get.ncdf(ncin, "tg", "units")
      fillvalue <- att.get.ncdf(ncin, "tg", "_FillValue")
      #fillvalue$value <- -999.9
    }
  }
  
  dim(climvar)
  
  close.ncdf(ncin)
  
  climvar[climvar < -99] <- NA # convert missing values (-999 to NA's)
  
  ###-----------------------------------------------------------###
  #       WRITE WINTER & SPRING DAILY DATA TO RASTERBRICK
  ###-----------------------------------------------------------###
  
  
  # split the time units string into fields, and check that dates are correct
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth = as.integer(unlist(tdstr)[2])
  tday = as.integer(unlist(tdstr)[3])
  tyear = as.integer(unlist(tdstr)[1])
  dates <- chron(t, origin = c(tday, tmonth, tyear), format="d m year")
  #dates
  
  # create a date lookup variable for the array so that appropriate date's raster layer can be pulled from the raster brick
  # datelookup <- data.frame(dates, year=years(dates), month=as.numeric(months(dates)), day=days(dates), year.month=paste(years(dates),as.numeric(months(dates)), sep="_"))
  # datelookup$month <- as.ordered(datelookup$month)
  
  datelookup <- data.frame(dates, year=years(dates), month=months(dates), day=days(dates), year.month=paste(years(dates),months(dates), sep="_"))
  
  if (a==1) {datelookup1980.1994 <- datelookup} else {datelookup1995.2013 <- datelookup}
  
  ### create a list of numeric indices for the time dimension for each year and month combination so that daily climate values can be summed (precipitation) or averaged (temperature) across months
  yearmonthindices <- list()
  
  for (i in 1:length(levels(datelookup$year))) {
    chooseyear <- levels(datelookup$year)[i]
    yearsubset <- subset(datelookup, year==chooseyear)
    
    for (j in 1:length(levels(yearsubset$month))) {
      choosemonth <- levels(yearsubset$month)[j]
      monthsubset <- subset(yearsubset, month==choosemonth)
      #yearmonthindices[a] <- row.names(monthsubset)
      yearmonthindices <- c(yearmonthindices,list(row.names(monthsubset)))
      names(yearmonthindices)[length(yearmonthindices)] <- paste(chooseyear, choosemonth, sep="")
    }
  }
  
  ####---- CREATE MONTHLY CLIMATE VALUES (precip=sum, temp=mean) ----####
  
  if (climvarname == "precipitation") {
    
    precip.monthtotal <- list()
    
    for (i in 1:length(yearmonthindices)) {
      
      monthtoextract <- as.numeric(yearmonthindices[[i]])
      clim.month <- climvar[, ,monthtoextract]
      
      climbrick.month <- brick(clim.month, xmn=-15, xmx=30, ymn=35, ymx=60, crs="+init=epsg:4326", transpose=TRUE)
      flipclimbrick.month <- flip(climbrick.month, direction="y")
      
      precip.monthtotal[[i]] <- stackApply(flipclimbrick.month, indices=rep(1,length(monthtoextract)), fun=sum, na.rm=FALSE) # ground-truthed against UK MetOffice totals (wettest spot in Scotland Dec 2013) and values seem to make sense
      
    }
    
    # turn by month/by year rasterbricks of total monthly precip into a single raster object, layer.1 = Jan in the first year of time series
    
    allmonthsyears.precip <- do.call(stack,precip.monthtotal)
    names(allmonthsyears.precip) <- names(yearmonthindices)
    
    if (a==1) {precip1980.1994 <- allmonthsyears.precip} else {precip1995.2013 <- allmonthsyears.precip}
    
  }
  
  if (climvarname == "temperature") {
    
    temp.monthmean <- list()
    
    for (i in 1:length(yearmonthindices)) {
      
      monthtoextract <- as.numeric(yearmonthindices[[i]])
      clim.month <- climvar[, ,monthtoextract]
      
      climbrick.month <- brick(clim.month, xmn=-15, xmx=30, ymn=35, ymx=60, crs="+init=epsg:4326", transpose=TRUE)
      flipclimbrick.month <- flip(climbrick.month, direction="y")
      
      temp.monthmean[[i]] <- stackApply(flipclimbrick.month, indices=rep(1,length(monthtoextract)), fun=mean)
      
    }
    
    # turn by month/by year rasterbricks of mean monthly temp into a single raster object, layer.1 = Jan in the first year of time series
    
    allmonthsyears.temp <- do.call(stack,temp.monthmean)
    names(allmonthsyears.temp) <- names(yearmonthindices)
    
    if (a==1) {temp1980.1994 <- allmonthsyears.temp} else {temp1995.2013 <- allmonthsyears.temp}
    
  }
  
}

if (climvarname=="precipitation") {
  precip1980.2013 <- stack(precip1980.1994,precip1995.2013)
}

if (climvarname=="temperature") {
  temp1980.2013 <- stack(temp1980.1994,temp1995.2013)
}

########################################
########################################
########################################
