######################################################
#
#     EXTRACT ENSEMBLES EUROPEAN CLIMATE DATA
#
######################################################

# uses ENSEMBLES European climate data downloaded as NETCDF files from http://www.ecad.eu/download/ensembles/download.php
# uses the 15-year chunk datafiles on the 0.25 degree grid for daily mean temperature (TG) and daily precipitation sum (RR) for 1995-2013

#### TO DO ####

# convert 3 dimensional data array into a 2 dimensional file with variables date, lat, long, and value (climvar) - extract only needed time series and lat/long cross-section

climvarname <- "precipitation"

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
  filename <- "rr_0.25deg_reg_1995-2013_v10.0.nc"
}

if (climvarname == "temperature") {
  filename <- "tg_0.25deg_reg_1995-2013_v10.0.nc"
}

###-----------------------------------------------------------###
#         OPEN NETCDF AND EXTRACT DATA
###-----------------------------------------------------------###

ncin <- open.ncdf(filename) 
print(ncin)

# summary function for ncdf objects
summary.ncdf <-
  function( control){
    for ( i in 1:control$nvars ){
      
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

# need time units 5723-6940 (1 Sep 2010 to 31 Dec 2013)
t <- get.var.ncdf(ncin, "time", start=5723, count=1218)
# t <- get.var.ncdf(ncin, "time") # all t
tunits <- att.get.ncdf(ncin, "time", "units")
nt <- dim(t)  # get the number of values and list them
nt


lon <- get.var.ncdf(ncin, "longitude", start=90, count=194)
#lon <- get.var.ncdf(ncin, "longitude") # all lon
nlon <- dim(lon)
head(lon)
# 0.25 degree grid lon units -18.125 lon[90] to 30.125 lon[283]
# need lon units -18.25 lon[45] to 30.25 lon[142] (0.5 degree grid)

lat <- get.var.ncdf(ncin, "latitude", start=39, count=102)
#lat <- get.var.ncdf(ncin, "latitude") # all lat
nlat <- dim(lat)
head(lat)
# 0.25 degree grid lat units 34.875 lat[39] to lat[140]
# need lat units 34.75 lat[20] to 60.25 lat[71]

print(c(nlon, nlat, nt))  # confirms the dimensions of the data

###--- EXTRACT CLIMATE VARIABLE FROM NETCDF ---###

if (climvarname == "precipitation") {
  if (fullseries) {climvar <- get.var.ncdf(ncin, "rr", start=c(90,39,1), count=c(194,102,6940))}
  #if (fullseries) {climvar <- get.var.ncdf(ncin, "rr", start=c(90,39,5723), count=c(194,102,1218))}
  dlname <- att.get.ncdf(ncin, "rr", "long_name")
  dunits <- att.get.ncdf(ncin, "rr", "units")
  fillvalue <- att.get.ncdf(ncin, "rr", "_FillValue")
  #fillvalue$value <- -999.9
}

if (climvarname == "temperature") {
  climvar <- get.var.ncdf(ncin, "tg", start=c(90,39,5723), count=c(194,102,1218))
  dlname <- att.get.ncdf(ncin, "tg", "long_name")
  dunits <- att.get.ncdf(ncin, "tg", "units")
  fillvalue <- att.get.ncdf(ncin, "tg", "_FillValue")
  #fillvalue$value <- -999.9
}

dim(climvar)

close.ncdf(ncin)

climvar[climvar < -99] <- NA # convert missing values (-999 to NA's)

###-----------------------------------------------------------###
#       WRITE WINTER & SPRING DAILY DATA TO RASTERBRICK
###-----------------------------------------------------------###


# split the time units string into fields, and check that dates are correct (check that correct date series is being subsetted - 01/09/2010 to 31/12/2013)
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
dates <- chron(t, origin = c(tday, tmonth, tyear), format="d m year")
dates

# number of layers = 1218, where layer1 = 1 Sep 2010 and layer1218 = 31 Dec 2013
# create a date lookup variable so that appropriate date's raster layer can be pulled from the raster brick
# datelookup <- data.frame(dates, year=years(dates), month=as.numeric(months(dates)), day=days(dates), year.month=paste(years(dates),as.numeric(months(dates)), sep="_"))
# datelookup$month <- as.ordered(datelookup$month)

datelookup <- data.frame(dates, year=years(dates), month=months(dates), day=days(dates), year.month=paste(years(dates),as.numeric(months(dates)), sep="_"))

# extract winter daily climate data (Nov-Feb) and average into new rasterbrick
winterdates <- which(datelookup$month == 11 | datelookup$month == 12 | datelookup$month == 1 | datelookup$month == 2)

winter <- list()
winter[[1]] <- winterdates[1:120]
winter[[2]] <- winterdates[121:241]
winter[[3]] <- winterdates[242:361]

# extract spring daily climate data (Mar-May) and average into new rasterbrick
springdates <- which(datelookup$month == 3 | datelookup$month == 4 | datelookup$month == 5)

spring <- list()
spring[[1]] <- springdates[1:92]
spring[[2]] <- springdates[93:184]
spring[[3]] <- springdates[185:276]

# replace netCDF fill values with NAs
climvar[climvar < -99] <- NA

# subset only the relevant dates to a new array
climvarwinter <- climvar[, ,winterdates] 
climvarspring <- climvar[, ,springdates]

# convert 3-dimensional array to a RasterBrick, but transpose the array (otherwise ends up with lat and long on the wrong axes) - then the resulting raster layers end up upside down (!), so need to flip the rasterbrick around the y axis as well to end up with the correct dimensions in the correct orientation

climbrickwinter <- brick(climvarwinter, xmn=-15, xmx=30, ymn=35, ymx=60, crs="+init=epsg:4326", transpose=TRUE)
flipclimbrickwinter <- flip(climbrickwinter, direction="y")

climbrickspring <- brick(climvarspring, xmn=-15, xmx=30, ymn=35, ymx=60, crs="+init=epsg:4326", transpose=TRUE)
flipclimbrickspring <- flip(climbrickspring, direction="y")

###-----------------------------------------------------------###
#       WRITE WINTER & SPRING AVERAGED DATA TO RASTERBRICK
###-----------------------------------------------------------###

# create a new RasterBrick with separate layers (using the indices argument in stackApply) representing the average winter and spring climate variables in each winter (2010-11, 2011-12, 2012-13) and spring (2011, 2012, 2013) of the analysis

winterindices <- c(rep(1,length(winter[[1]])), rep(2,length(winter[[2]])), rep(3,length(winter[[3]])))

springindices <- c(rep(1,length(spring[[1]])), rep(2,length(spring[[2]])), rep(3, length(spring[[3]])))

if (climvarname == "precipitation") {
  precipbrickwinter <- flipclimbrickwinter
  precipbrickspring <- flipclimbrickspring
  
  mean.winter.precip <- stackApply(precipbrickwinter, indices=winterindices, fun=mean)  
  mean.spring.precip <- stackApply(precipbrickspring, indices=springindices, fun=mean)
  
#   plot(mean.winter.precip)
#   plot(mean.spring.precip)
  
}

if (climvarname == "temperature") {
  tempbrickwinter <- flipclimbrickwinter
  tempbrickspring <- flipclimbrickspring
  
  mean.winter.temp <- stackApply(tempbrickwinter, indices=winterindices, fun=mean)  
  mean.spring.temp <- stackApply(tempbrickspring, indices=springindices, fun=mean)
  
#   plot(mean.winter.temp)
#   plot(mean.spring.temp)
  
}

# group climate data RasterBricks into a list of RasterBricks once it has all been run

if (("mean.winter.precip" %in% ls()) & ("mean.winter.temp" %in% ls())){
  climate.data <- list(mean.winter.precip, mean.spring.precip, mean.winter.temp, mean.spring.temp)
  names(climate.data) <- c("mean.winter.precip", "mean.spring.precip", "mean.winter.temp", "mean.spring.temp")
}



# ######################################################
# ######################################################
# 
# #### the following code converts the climate netCDF data to a long dataframe with all combinations of lat/long, date, and associated climate variable
# 
# length(na.omit(as.vector(climvar[, , 1])))
# 
# m <- 1
# climvar.slice <- climvar[, , m]
# image(lon, lat, climvar.slice, xlim=c(-18,30), ylim=c(35,60), col = rev(brewer.pal(10, "RdBu")))
# 
# # Convert the nlon by nlat by nt array into a nlon by nlat by nt matrix. (This will work if the NetCDF data set was written as a CF-compliant data set, with arrays dimensioned as in Fortran, as nlon x nlat x nt arrays) First, create a long vector climvar.vec.long using the as.vector() reshaping function, and verify its length, which should be 98 * 52 * 1218 = 6206928.
# climvar.vec.long <- as.vector(climvar)
# length(climvar.vec.long)
# 
# # Then reshape that vector into a 5096 by 1218 matrix using the matrix() function, and verify its dimensions, which should be 5096 by 1218.
# 
# climvar.mat <- matrix(climvar.vec.long, nrow = nlon * nlat, ncol = nt)
# dim(climvar.mat)
# 
# # Create the second data frame from the climvar.mat matrix.
# 
# lonlat <- expand.grid(lon, lat)
# climvar2 <- data.frame(cbind(lonlat, climvar.mat))
# names(climvar2) <- c("lon", "lat", t)
# 
# climvar3 <- melt(climvar2, id=c("lon", "lat"))
# climvar4 <- na.omit(climvar3)
# dates2 <- chron(as.numeric(as.character(climvar4$variable)), origin = c(tmonth, tday, tyear))
# 
# if (climvarname == "precipitation") {
#   climvar.final <- data.frame(long=climvar4$lon, lat=climvar4$lat, date=dates2, precip.mm=climvar4$value)
#   write.csv(climvar.final, "precipitation 20100901-20131231.csv", row.names=FALSE)
# }
# 
# if (climvarname == "temperature") {
#   climvar.final <- data.frame(long=climvar4$lon, lat=climvar4$lat, date=dates2, temp.degreesC=climvar4$value)
#   write.csv(climvar.final, "temperature 20100901-20131231.csv", row.names=FALSE)
#   
# }
