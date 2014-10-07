#########################################################
#
#     PROTECTED AREA DETAILS - FURTHER ANALYSIS
#
#########################################################

library(sp) # converts R dataframes into spatially explicit dta
library(maptools) # allows import of shapefiles into R
library(rgdal)
library(rgeos)
library(raster)
library(geosphere)
library(shape)
library(ggplot2)
library(grid)

Mac <- FALSE

if(.Platform$OS =='windows') cluster <- FALSE
if(.Platform$OS=='unix' & !Mac) cluster <- TRUE

if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  
}

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")

if (!cluster) {
  GISwd <- c("C:/Users/samf/Documents/GIS/cuckoos/protected areas")
}

if (cluster) {
  GISwd <- c("/users1/samf/cuckoos")
}

### load PA shapefile, subset out international and national PAs
epsg4326 <- CRS("+init=epsg:4326")
epsg3035 <- CRS("+init=epsg:3035")

PA <- readOGR(GISwd, "terrestrial_PAs_Europe_GCS_WGS_1984_final")
PA.epsg4326 <- spTransform(PA, epsg4326)

natPA <- subset(PA.epsg4326, !grepl("Ramsar | Directive", desig_eng))
internatPA <- subset(PA.epsg4326, grepl("Ramsar | Directive", desig_eng))

### load ESRI shapefile layer of countries that is used for ArcGIS layer in map of presences/absences (from G:/Common Themes/ESRIData)
# country code for Montenegro is MON in country GIS shapefile, but is MNE elsewhere - convert to "MNE"
country <- readOGR(dsn="C:/Users/samf/Documents/GIS/cuckoos", layer="ESRI_country_shapefile")
country <- spTransform(country, epsg4326)
country@data$code <- factor(country@data$GMI_CNTRY, levels=c(levels(country@data$GMI_CNTRY), "MNE"))
country@data[which(country@data$CNTRY_NAME=="Montenegro"), "code"] <- "MNE"

### load presences
setwd(paste(datawd, "/data for analysis/", sep=""))
present <- read.csv("presence data all variables points all with country data.csv", header=T)


### country names and codes
countries <- levels(present$usecountry)
countrycode <- country@data[country@data$CNTRY_NAME %in% countries, c("CNTRY_NAME","code")]
colnames(countrycode) <- c("countryname","country")

### load PA details for PAs that overlapped with stopover points
allPAdetails <- list()

setwd(paste(datawd, "/corine PA elevation spei extracted values points/protected areas details/presences", sep=""))

for (a in 1:length(list.files())) {
  allPAdetails[[a]] <- read.csv(list.files()[a], header=T)
}

alldat <- do.call(rbind, allPAdetails)

SWdat <- subset(alldat, strategy=="SW")
SWdat <- droplevels(SWdat)
SEdat <- subset(alldat, strategy=="SE")
SEdat <- droplevels(SEdat)

# SW
internatPA.SW <- subset(internatPA, (internatPA$wdpaid %in% SWdat$wdpaid))
internatPA.SW@data <- droplevels(internatPA.SW@data)
internatPA.SW <- spTransform(internatPA.SW, epsg3035)

natPA.SW <- subset(natPA, (natPA$wdpaid %in% SWdat$wdpaid))
natPA.SW@data <- droplevels(natPA.SW@data)
natPA.SW <- spTransform(natPA.SW, epsg3035)

# SE
internatPA.SE <- subset(internatPA, (internatPA$wdpaid %in% SEdat$wdpaid))
internatPA.SE@data <- droplevels(internatPA.SE@data)
internatPA.SE <- spTransform(internatPA.SE, epsg3035)


natPA.SE <- subset(natPA, (natPA$wdpaid %in% SEdat$wdpaid))
natPA.SE@data <- droplevels(natPA.SE@data)
natPA.SE <- spTransform(natPA.SE, epsg3035)

####==== AREA OF USED PROTECTED AREAS ====####

# calculate area (in square kms) of PAs on SW and SE routes that birds USE
internatPA.SW@data$area <- gArea(internatPA.SW, byid=TRUE)/(1000*1000)
natPA.SW@data$area <- gArea(natPA.SW, byid=TRUE)/(1000*1000)
internatPA.SE@data$area <- gArea(internatPA.SE, byid=TRUE)/(1000*1000)
natPA.SE@data$area <- gArea(natPA.SE, byid=TRUE)/(1000*1000)

area <- list()
area[[1]] <- data.frame(area=gArea(internatPA.SW, byid=TRUE)/(1000*1000), route="SW", desig="International", country=internatPA.SW$country)
area[[2]] <- data.frame(area=gArea(natPA.SW, byid=TRUE)/(1000*1000), route="SW", desig="National", country=natPA.SW$country)
area[[3]] <- data.frame(area=gArea(internatPA.SE, byid=TRUE)/(1000*1000), route="SE", desig="International", country=internatPA.SE$country)
area[[4]] <- data.frame(area=gArea(natPA.SE, byid=TRUE)/(1000*1000), route="SE", desig="National", country=natPA.SE$country)

usedPAs <- do.call(rbind, area)
usedPAs <- merge(usedPAs, countrycode, by="country")

### BOXPLOT OF AREA OF USED PROTECTED AREAS, BY COUNTRY
qplot(y=area, x=country, data=usedPAs, geom="boxplot")
setwd(paste(outputwd, "/autumn stopover GLMM results/", sep=""))
ggsave("boxplot area of used protected areas.jpg", width=15, height=12, units="in")

### HISTOGRAM OF AREA OF USED PAs, by COUNTRY
p <- ggplot(usedPAs, aes(x=area)) + geom_density() # + geom_histogram(aes(y=..density..), colour="black", fill="white")
p + facet_wrap(~ country)
setwd(paste(outputwd, "/autumn stopover GLMM results/", sep=""))
ggsave("histogram areas of used protected areas.jpg", width=15, height=12, units="in")

####==== AREA OF ALL PROTECTED AREAS (IN COUNTRIES USED BY CUCKOOS) ====####)


internatPA.sub <- subset(internatPA, (internatPA$country %in% countrycode$country))
internatPA.sub@data <- droplevels(internatPA.sub@data)
natPA.sub <- subset(natPA, (natPA$country %in% countrycode$country))
natPA.sub@data <- droplevels(natPA.sub@data)

PA.sub <- rbind(internatPA.sub, natPA.sub)
allPAs <- merge(PA.sub@data, countrycode, by="country")

# qplot(y=rep_area, x=countryname, geom="boxplot", data=allPAs)
# 
# p <- ggplot(allPAs, aes(x=rep_area)) + geom_histogram(binwidth=10)
# p + facet_wrap(~ countryname)
# 
# verylargePAs <- allPAs[allPAs$rep_area > 2000,]
# hist(verylargePAs$rep_area)
# smallPAs <- allPAs[allPAs$rep_area <= 200,]
# hist(smallPAs$rep_area, breaks=50)

### AREA OF PAs by designation level, by country
totalareas <- rbind(data.frame(area=tapply(internatPA.sub$rep_area, internatPA.sub$country, sum), country=levels(internatPA.sub$country), desig="International"), data.frame(area=tapply(natPA.sub$rep_area, natPA.sub$country, sum), country=levels(natPA.sub$country), desig="National"))

### merge with country data to extra area of countries that have PAs
countrydata <- with(country@data, data.frame(country=GMI_CNTRY, name=CNTRY_NAME, sqkm=SQKM_CNTRY))
MNEdata <- data.frame(country="MNE", name="Montenegro", sqkm=countrydata[countrydata$name=="Montenegro","sqkm"])
countrydata <- rbind(countrydata, MNEdata)
areas <- merge(totalareas, countrydata)
areas <- droplevels(areas)

### proportion of country area that is PA
propPA.country <- data.frame(areas, propPA=(areas$area/areas$sqkm))
propPA.country <- propPA.country[order(propPA.country$desig, propPA.country$country),]


### PLOT proportional area PA by country
plottheme <- theme(
  axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
  axis.text.x = element_text(size=14, colour="black", angle=45, vjust=0.5), 
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
  panel.grid.minor = element_blank(),
  strip.text.x=element_text(face="bold", size=14)
)

labels <- labs(x="Country", y="Proportion of total area covered by protected area", title="Proportion of total country area covered by protected area")


### PLOT

p <- ggplot(propPA.country, aes(x=name, y=propPA, fill=desig)) + geom_bar(stat="identity", position="dodge")
p + labels + plottheme + scale_fill_manual(values=c("#003300", "#339966"))
setwd(paste(outputwd, "/autumn stopover GLMM results/", sep=""))
ggsave("proportion of total country area covered by PA.jpg", width=15, height=12, units="in")

p <- ggplot(propPA.country, aes(x=name, y=propPA, fill=desig)) + geom_bar(stat="identity")
p + labels + plottheme + scale_fill_manual(values=c("#003300", "#339966"))
setwd(paste(outputwd, "/autumn stopover GLMM results/", sep=""))
ggsave("proportion of total country area covered by PA stacked.jpg", width=15, height=12, units="in")
