#########################################################################
#
# Rasterize vector data for European protected areas
#
# Samantha Franks
# 24 Dec 2013
#
#########################################################################

rm(list=ls())

### LOAD PACKAGES
library(sp)
library(raster)
library(rgdal)
library(rgeos)

### Load shapefile

# set working directory (alter if working on cluster or not)

cluster <- FALSE

if (!cluster)
GISwd <- c("D:/Sam Franks/GIS/cuckoos")

if (cluster)
GISwd <- c("/users1/samf/cuckoos")

PAs <- readOGR(GISwd, "terrestrial PAs mainland W Europe corine countries only EPSG 3035")

### Clip shapefile to extent of corine map
# crop with SPDF (rather than raster) does I think the same as gIntersection in library(rgeos)
#clipPA <- drawExtent()
clipPA <- extent(2483500,5890400,1276600,4286800) # same extent as Europeraster (corine clipped raster layer)
EuropePA <- crop(PAs,clipPA)

### Rasterize the SpatialPolygons shapefile - 50m x 50m raster output
## Set up a raster "template" to use in rasterize()
extpoly <- extent(EuropePA)
col50 <- ceiling((extpoly@xmax-extpoly@xmin)/50) # create number of columns for raster (aim for approx 50m blocks)
row50 <- ceiling((extpoly@ymax-extpoly@ymin)/50) # create number of rows
newraster50 <- raster(extpoly, ncol=col50, nrow=row50)
  
## Rasterize the shapefile polygons

rastertime50 <- system.time({

PArasterized50 <-rasterize(EuropePA, newraster50)
  
})

rastertime50

# colors <- c("white",rep("blue",9407))
# rpoly50@legend@colortable <- colors

#newrpoly <- crop(rpoly,clipUKr)
# default of rpoly colortable is logical(0)

setwd(GISwd)
writeRaster(PArasterized50, filename="Europe PA raster 50m x 50m.tif", format="GTiff", overwrite=TRUE)

######################################################################

################################### TEST  CODE ###############################

######################################################################

# ### Load shapefile
# 
# GISwd <- c("D:/Sam Franks/GIS/cuckoos")
# 
# testshp <- readOGR(GISwd, "terrestrial PAs UK EPSG 4326")
# 
# #proj4string(testshp) <- CRS("+init=epsg:4326")
# 
# PAtrans <- spTransform(testshp,CRS("+init=epsg:3035"))
# 
# ## crop polygon - crop with SPDF does I think the same as gIntersection in library(rgeos)
# #clipPA <- drawExtent()
# clipPA <- extent(3039340,3879513,3008405,4308673)
# UK.PAs <- crop(PAtrans,clipPA)
# 
# ## create SpatialLinesDataFrame
# UKPAlines <- as(UK.PAs, "SpatialLinesDataFrame")
# 
# ### Create raster using SpatialLines
# ## Set up a raster "template" to use in rasterize()
# extlines <- extent(UKPAlines)
# col <- ceiling((extlines@xmax-extlines@xmin)/100) # create number of columns for raster (aim for approx 100m blocks)
# row <- ceiling((extlines@ymax-extlines@ymin)/100) # create number of rows
# newrasterlines <- raster(extlines, ncol=col, nrow=row)
#   
# ## Rasterize the shapefile lines
# 
# rastertimelines <- system.time({
# 
# rlines <-rasterize(UKPAlines, newrasterlines)
#   
# })
# 
# ### Create raster using SpatialPolygons - 100m x 100m
# ## Set up a raster "template" to use in rasterize()
# extpoly <- extent(UK.PAs)
# col <- ceiling((extpoly@xmax-extpoly@xmin)/100) # create number of columns for raster (aim for approx 100m blocks)
# row <- ceiling((extpoly@ymax-extpoly@ymin)/100) # create number of rows
# newrasterpoly <- raster(extpoly, ncol=col, nrow=row)
#   
# ## Rasterize the shapefile polygons
# 
# rastertimepoly <- system.time({
# 
# rpoly <-rasterize(UK.PAs, newrasterpoly)
#   
# })
# 
# colors <- c("white",rep("blue",9407))
# rpoly@legend@colortable <- colors
# 
# #newrpoly <- crop(rpoly,clipUKr)
# # default of rpoly colortable is logical(0)
# 
# setwd(GISwd)
# writeRaster(rpoly, filename="UK PAs test raster 100m x 100m.tif", format="GTiff", overwrite=TRUE)
# 
# ### Create raster using SpatialPolygons - 50m x 50m
# ## Set up a raster "template" to use in rasterize()
# extpoly <- extent(UK.PAs)
# col50 <- ceiling((extpoly@xmax-extpoly@xmin)/50) # create number of columns for raster (aim for approx 100m blocks)
# row50 <- ceiling((extpoly@ymax-extpoly@ymin)/50) # create number of rows
# newrasterpoly50 <- raster(extpoly, ncol=col50, nrow=row50)
#   
# ## Rasterize the shapefile polygons
# 
# rastertimepoly50 <- system.time({
# 
# rpoly50 <-rasterize(UK.PAs, newrasterpoly50)
#   
# })
# 
# colors <- c("white",rep("blue",9407))
# rpoly50@legend@colortable <- colors
# 
# #newrpoly <- crop(rpoly,clipUKr)
# # default of rpoly colortable is logical(0)
# 
# setwd(GISwd)
# writeRaster(rpoly50, filename="UK PAs test raster 50m x 50m.tif", format="GTiff", overwrite=TRUE)
# 
# 
# ### Create raster using SpatialPolygons - 200m x 200m
# ## Set up a raster "template" to use in rasterize()
# extpoly <- extent(UK.PAs)
# col200 <- ceiling((extpoly@xmax-extpoly@xmin)/200) # create number of columns for raster (aim for approx 100m blocks)
# row200 <- ceiling((extpoly@ymax-extpoly@ymin)/200) # create number of rows
# newrasterpoly200 <- raster(extpoly, ncol=col200, nrow=row200)
#   
# ## Rasterize the shapefile polygons
# 
# rastertimepoly200 <- system.time({
# 
# rpoly200 <-rasterize(UK.PAs, newrasterpoly200)
#   
# })
# 
# ### Create raster using SpatialPolygons - 500m x 500m
# ## Set up a raster "template" to use in rasterize()
# extpoly <- extent(UK.PAs)
# col500 <- ceiling((extpoly@xmax-extpoly@xmin)/500) # create number of columns for raster (aim for approx 100m blocks)
# row500 <- ceiling((extpoly@ymax-extpoly@ymin)/500) # create number of rows
# newrasterpoly500 <- raster(extpoly, ncol=col500, nrow=row500)
#   
# ## Rasterize the shapefile polygons
# 
# rastertimepoly500 <- system.time({
# 
# rpoly500 <-rasterize(UK.PAs, newrasterpoly200)
#   
# })
# 
# 
# ###################################
# ###################################
# 
# proj4string(shp) <- CRS("+init=epsg:4326")
# PAtrans <- spTransform(shp,CRS("+init=epsg:3035"))
# 
# 
# ### create color table for raster
# 
# colors <- c("white",rep("blue",9407))
# 
# 
# rpoly@legend@colortable <- colors
# 
# setwd(GISwd)
# tiff("UK SPA test plot 50m x 50m .tiff", width=3000, height=3000, units="px", res=300)
# plot(rpoly)
# dev.off()
# 
# rr@legend@colortable <- logical(0)
# 
# tiff("UK SPA test plot 2.tiff", width=3000, height=3000, units="px", res=300)
# plot(rr, col="blue")
# dev.off()
# 
# tiff("UK SPA test plot vector.tiff", width=3000, height=3000, units="px", res=300)
# plot(UK.PAs, col="blue")
# dev.off()
# 
# #################
# 
# PAs <- readOGR(GISwd, "terrestrial PAs mainland W Europe corine countries only EPSG 3035")
# 
# 
# rastertime <- system.time({
#   
# ## Set up a raster "template" to use in rasterize()
# ext <- extent(1500000,7400000,,5500000)
# newraster <- raster(ext, ncol=46000, nrow=59000)
#   
# ## Rasterize the shapefile
# rr <-rasterize(testshp, newraster)
#   
# })
# 
# ## A couple of outputs
# writeRaster(rr, "teow.asc")
# plot(rr)
# 
# testshp <- readOGR("C:/Users/samf/Documents/GIS/cuckoos/official_teow","wwf_terr_ecos")
# 
# ## Set up a raster "template" to use in rasterize()
# ext <-  extent (-95, -50, 24, 63)
# xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
# n <- 5
# ras <- raster(ext, ncol=xy[1]*5, nrow=xy[2]*5)
# 
# ## Rasterize the shapefile
# rr <-rasterize(testshp, ras)
# 
# ## A couple of outputs
# writeRaster(rr, "teow.asc")
# plot(rr)

