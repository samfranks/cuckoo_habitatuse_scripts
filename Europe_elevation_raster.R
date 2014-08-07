########################################################
#
#  CREATE ELEVATION RASTER FOR EUROPE
#
# Samantha Franks
# 2 May 2014
#
########################################################

library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
library(rgeos)
library(raster)

load("/users1/samf/cuckoos/PA corine extraction.RData")

### ------- Elevation data ------- ###

# corine.crs <- CRS("+init=epsg:3035")

countries <- levels(PA.subset@data$country)[which(levels(PA.subset@data$country) %in% getData('ISO3'))]

countries2 <- countries[-which(countries=="MLT")]

elev <- list()

for (i in 1:length(countries2)){
  elev[[i]] <- getData("alt", country=countries2[i], mask=TRUE)
}

elev[[27]] <- merge(elev[[27]][[1]], elev[[27]][[2]])
fullelev <- do.call(merge, elev)

# plot(newelev, col=terrain.colors(20))

# newelev <- projectRaster(fullelev, crs=corine.crs)

# writeRaster(newelev, filename="Europe elevation raster epsg 3035.tif", format="GTiff", overwrite=TRUE)
writeRaster(fullelev, filename="Europe elevation raster.tif", format="GTiff", overwrite=TRUE)