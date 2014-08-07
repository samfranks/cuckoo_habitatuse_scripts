##########################################################
#
#  Cuckoo resampled points and land cover value extractions
#
#	Samantha Franks
#	28 Nov 2013
#
#
##########################################################

library(adehabitat) # 
library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
#library(shapefiles) # allows you to read and write ESRI shapefiles
library(rgeos)
library(raster)
library(geosphere)
library(adehabitatHR)
library(rworldmap)
library(rworldxtra)
library(NCStats)

### -------- Working directories --------- ###

# GISwd <- c("C:/Users/samf/Documents/R/projects/cuckoos/ArcGIS/all original layers WGS84 export")
corinewd <- c("C:/Users/samf/Documents/R/projects/cuckoos/data/corine raster 100x100")

###-----------------------------------------------------------###
#         Load base layers
###-----------------------------------------------------------###

### CORINE layer ###

# projection and datum information for all maps
# epsg: 3395 - +proj=merc +ellps=wGS84 +datum=WGS84
# epsg: 4326 - +proj=latlong +ellps=wGS84 +datum=WGS84
### epsg: 3035 - LAEA CRS, used for the Corine Land Cover 2006 data

corine.crs <- CRS("+init=epsg:3035")

### ------- CORINE land cover raster file ------- ###
setwd(corinewd)
r <- raster("g100_06.tif")

Europebox <- extent(2483500,5890375,1276625,4286750) # bounding box extent for relevant parts of Europe (no Scandanavia or Turkey)
Europeraster <- crop(r,Europebox)

### ------- Europe/Africa national boundaries layer ------- ###

dat.world <- getMap(resolution="high")

Eur.Afr <- dat.world[which(dat.world$REGION=="Europe" | dat.world$REGION=="Africa"),]
Eur.Afr <- Eur.Afr[which(Eur.Afr$ADMIN!="Russia" & Eur.Afr$ADMIN!= "Azerbaijan" & Eur.Afr$ADMIN!="Gaza" & Eur.Afr$ADMIN!="Georgia" & Eur.Afr$ADMIN!="Israel" & Eur.Afr$ADMIN!="West Bank" & Eur.Afr$ADMIN!="Armenia"),]
Eur.Afr <- Eur.Afr[which(is.na(Eur.Afr$TERR_)),]

Eur.Afr <- spTransform(Eur.Afr,CRS("+proj=longlat +datum=WGS84")) # transform map to exactly the same projection as cuckoo newlongs/newlats data

Eur.only <- Eur.Afr[which(Eur.Afr$REGION=="Europe"),]
Afr.only <- Eur.Afr[which(Eur.Afr$REGION=="Africa"),]

### -------- Protected area network -------- ###

###------------------------------------------------------------###
#         Clip cuckoo stopovers to Europe map
###------------------------------------------------------------###

# create list to hold centroids of stopover points for all cuckoos

corine.values <- list()

# source file creates an object called "original", which are the tidied original datapoints for all individuals
setwd("C:/Users/samf/Documents/R/projects/cuckoos/scripts/")
source("cuckoo original data tidying source.R")

# for loop for each cuckoo file starts here
for (a in 1:31) {
  
  ### --- RESAMPLED cuckoo data: Load and add geoinfo --- ###
  
  setwd("C:/Users/samf/Documents/R/projects/cuckoos/data/resampled data with uncertainty 100 bootstraps + distance, movement groups, etc")
  cuckoofiles <- list.files()
  
  cuckoo <- read.csv(cuckoofiles[a], header=T)
  
  # long/lat are in CRS("+proj=longlat +datum=WGS84")
  # newlongs/newlats are in CRS("+init=epsg:3395")
  # centroidlong/lat are in CRS("+init=epsg:3395")
  
  coordinates(cuckoo) <- c("newlongs","newlats")
  proj4string(cuckoo) <- CRS("+init=epsg:3395")
  cuckoo <- spTransform(cuckoo,CRS("+proj=longlat +datum=WGS84")) # transform newlongs/newlats from World Mercator/WGS84 to longlat WGS84
  
  ### over() retrieves values of Eur.Afr at points described by cuckoo
  # extracts continent and administrative country names from Eur.Afr
  geoinfo <- over(cuckoo, Eur.Afr)[,c("REGION","ADMIN")]
  colnames(geoinfo) <- c("continent", "country")
  geoinfo <- droplevels(geoinfo)
  
  ### add a new factor for continent/country NAs, which should be points at sea
  atsea.pts <- which(is.na(geoinfo$continent))
  
  continent.factors <- c(levels(geoinfo$continent), "at sea")
  levels(geoinfo$continent) <- continent.factors
  geoinfo$continent[atsea.pts] <- "at sea"
  
  country.factors <- c(levels(geoinfo$country),"at sea")
  levels(geoinfo$country) <- country.factors
  geoinfo$country[atsea.pts] <- "at sea"
  
  ### add geographic info to cuckoo dataset and redefine projection of newlongs/newlats
  fullcuckoo <- data.frame(cuckoo, geoinfo)
  coordinates(fullcuckoo) <- c("newlongs","newlats")
  proj4string(fullcuckoo) <- CRS("+proj=longlat +datum=WGS84")
  
  ### create Europe only (+ at sea points) and Africa only data sets
  Eurcuckoo <- fullcuckoo[which(fullcuckoo$continent=="Europe" | fullcuckoo$continent=="at sea"),]
  Afrcuckoo <- fullcuckoo[which(fullcuckoo$continent=="Africa"),]
  
  ###------------------------------------------------------------###
  #         Trim European data - autumn stopovers only
  ###------------------------------------------------------------###
  
  # subset out only stopovers in Europe
  # do not include the first "stopover", which we assume is the capture and therefore breeding location
  # ???????????? other early stopovers (e.g. #2, or #3) may also need to be removed if they constitute movements around the breeding area (unlikely though, since distance threshold for a "stopover" is < 20km, and anything > 20km would likely consitute a migration movement away from the breeding area)
  # include only southward migration stopovers (month >= June 1)
  
  ### RESAMPLED CUCKOO DATA
  
  Eurcuckoo$mgroup <- as.factor(Eurcuckoo$mgroup)
  stopovers <- Eurcuckoo[which(Eurcuckoo$mtype=="C" | Eurcuckoo$mtype=="S" | is.na(Eurcuckoo$mtype)),] # includes capture (breeding) location, stopovers, as well as last known location (if in Europe)
  #stopovers <- stopovers[which(stopovers$mgroup!=1),] # this line removes the capture occasion
  autumn <- stopovers[which(stopovers$julian >= 131),] # all julian dates after May 10 (removes spring stopovers)
  autumn@data <- droplevels(autumn@data)
  
  autumn$mgroup <- as.factor(autumn$mgroup)
  
  # can add other subsetting variables here, and change whatever "dat" is called
  
  dat <- autumn
  
  # Transform cuckoo coordinates to same CRS used for Corine layers using spTransform() in library(rgdal)
  
  newdat <- spTransform(dat, CRS=corine.crs)
  
  ###-------------------------------------------------------------###
  #   Land-cover point value extractions
  ###-------------------------------------------------------------###
  
  corine.values[[a]] <- data.frame(newdat, corine.values=extract(Europeraster, newdat@coords))
    
}

setwd("C:/Users/samf/Documents/R/projects/cuckoos/data/resampled data corine extracted values")

for (i in 1:31){
  write.csv(corine.values[[i]], paste(corine.values[[i]]$name[1],".csv", sep=""), row.names=FALSE)
}


###---------------------------------------------------------###
#
#  Barplots of land-cover types used
#
###---------------------------------------------------------###

### --- FUNCTION TO PLOT CUCKOO HABITAT USE, MULTIPLE SITES POOLED --- ###

habsummary <- function(dataset,landcover){ # datasets required are point land cover values, clclegend data (landcover), and par settings for specific # of stopovers
  
  dataset <- droplevels(dataset)
  dataset$mgroup <- as.factor(as.character(dataset$mgroup))
  
  # calculate proportion of each habitat type used
  prophab <- prop.table(table(dataset$corine.name))
  
  # produce colour labels for barplot according to RGB colour codes of CLC legend
  landcolours <- rgb(landcover[which(landcover$LABEL4 %in% levels(dataset$corine.name)),c("rgb1","rgb2","rgb3")], max=255)
  
  # plot the proportional habitat use
  
  par(mar=c(9,4,2,1), oma=c(0,0,2,0)) # specific par settings for number of plots to produce
  
  x <- barplot(prophab, names.arg=c(""), ylim=c(0,max(prophab)+0.05), col=landcolours, cex.axis=0.8, las=2)
  text(x, par("usr")[3] - 0.01, srt = 45, adj = 1, xpd = TRUE, cex=0.7, labels = levels(dataset$corine.name))
  
  # if dataset has a single stopover, then label is for that specific stopover; otherwise, label is for Europe generally
  if (length(levels(dataset$mgroup)) == 1) {
    mtext(paste(dataset$name[1], ", habitat use in ", dataset$country[1], " ",dataset$year[1],"\nLast stopover? = ", dataset$laststop[1], sep=""), cex=0.8,font=2,line=0)
  } else {
    mtext(paste(dataset$name[1], ", overall habitat use in Europe", sep=""), cex=0.8,font=2,line=1)
  }
  
}

### --- FUNCTION TO PLOT CUCKOO HABITAT USE, SINGLE SITES, MULTIPLE PLOTS PER PAGE --- ###

habsummary.stopovers <- function(dataset,landcover){ # datasets required are point land cover values, clclegend data (landcover), and par settings for specific # of stopovers
  
  dataset <- droplevels(dataset)
  dataset$mgroup <- as.factor(as.character(dataset$mgroup))
  
  # calculate proportion of each habitat type used
  prophab <- prop.table(table(dataset$corine.name))
  
  # produce colour labels for barplot according to RGB colour codes of CLC legend
  landcolours <- rgb(landcover[which(landcover$LABEL4 %in% levels(dataset$corine.name)),c("rgb1","rgb2","rgb3")], max=255)
  
  # plot the proportional habitat use
  
  x <- barplot(prophab, names.arg=c(""), ylim=c(0,max(prophab)+0.05), col=landcolours, cex.axis=0.8, las=2)
  text(x, par("usr")[3] - 0.03, srt = 45, adj = 1, xpd = TRUE, cex=0.7, labels = levels(dataset$corine.name))
  
  # if dataset has a single stopover, then label is for that specific stopover; otherwise, label is for Europe generally
  mtext(paste(dataset$name[1], ", habitat use in ", dataset$country[1], " ",dataset$year[1],"\nLast stopover? = ", dataset$laststop[1], sep=""), cex=0.5,font=2,line=1)
  
}

### --- SET-UP DATA --- ###

### import corine legend
setwd("C:/Users/samf/Documents/R/projects/cuckoos/data/")
clclegend <- read.csv("clc_legend.csv", header=T)

for (a in 1:31) {
  ### select cuckoo dataset
  dat <- corine.values[[a]]
  
  ### identify last stopovers in Europe
  # identify last stopovers before Sahara crossing
  # for each level of year, identify max(mgroup)
  
  # convert mgroup to integer, via character
  dat$mgroup <- as.integer(as.character(dat$mgroup))
  
  yearlevels <- levels(as.factor(dat$year))
  maxmgroup.byyear <- rep(NA, length(yearlevels))
  
  for (i in 1:length(yearlevels)){
    maxmgroup.byyear[i] <- max(dat[dat$year==yearlevels[i], "mgroup"]) # extract max mgroup number for each level of year
  }
  
  laststop <- factor(levels=c("N","Y"))
  laststop[which(dat$mgroup %in% maxmgroup.byyear)] <- "Y"
  laststop[which(!dat$mgroup %in% maxmgroup.byyear)] <- "N"
  dat <- data.frame(dat,laststop)
  
  ### add names of corine land types to cuckoo dataset
  corine.name <- rep(NA,nrow(dat))
  
  for (i in 1:nrow(dat)){
    corine.name[i] <- as.character(clclegend[which(clclegend$GRID_CODE %in% dat[i,"corine.values"]), "LABEL4"])
  }
  
  newland <- data.frame(dat, corine.name)
  
  # order land cover levels in dataset so they appear in order of grid codes in barplot
  orderedcorine <- clclegend$LABEL4[which(clclegend$LABEL4 %in% levels(newland$corine.name))] 
  newland$corine.name <- factor(newland$corine.name, orderedcorine, ordered=TRUE)
  
  ### --- Land-use datasets --- ###
  
  # create and set wd for output habitat plots
  # create new directory for specific bird to put output maps in, and setwd() as this new directory
  
  #dir.create(paste("C:/Users/samf/Documents/R/projects/cuckoos/output/Europe habitat use/", levels(newland$name), sep=""))
  #outputdir <- paste("C:/Users/samf/Documents/R/projects/cuckoos/output/Europe habitat use/", levels(newland$name), sep="")
  setwd("C:/Users/samf/Documents/R/projects/cuckoos/output/Europe habitat use/")
  
  # all sites in Europe, breeding + stopovers
  tiff(paste(levels(newland$name)," - Europe.tiff"),res=300,height=2000,width=2400,units="px")
  
  habsummary(newland,clclegend)
  
  dev.off()
  
  # all stopovers, individually
  # convert mgroup to level and split data by stopover (level of mgroup)
  stopovers <- newland
  stopovers$mgroup <- as.factor(as.character(stopovers$mgroup))
  stopovers$mgroup <- reorder(stopovers$mgroup)
  bymgroup <- stopovers$mgroup
  splitstopovers <- split(stopovers, list(bymgroup))
  
  
  tiff(paste(levels(newland$name)," - individual stopovers.tiff"),res=300,height=2000,width=3000,units="px")
  
  if (length(splitstopovers) < 8){
  par(mfrow=c(2, ceiling(length(splitstopovers)/2)), mar=c(9,4,3,2), oma=c(0,0,2,1))
  } else {
    par(mfrow=c(3, ceiling(length(splitstopovers)/3)), mar=c(9,4,3,2), oma=c(0,0,2,1))

  }

  for (i in 1:length(splitstopovers)){
    
    habsummary.stopovers(splitstopovers[[i]],clclegend)
  }
  
  dev.off()
  
  # function to plot habitat use at all stopovers - use lapply to apply across list elements
  
}
