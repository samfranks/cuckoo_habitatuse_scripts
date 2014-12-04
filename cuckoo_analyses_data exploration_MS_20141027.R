##########################################################
#
#   CUCKOO HABITAT USE: DATA EXPLORATION & VISUALIZATION
#
#  Samantha Franks
#  11 April 2014
#  24 May 2014
#  17 Jul 2014: new random polygons generating 10 pseudoabsences/presence
#  5 Aug 2014: analysis using randomly generated points instead of polygons for absences, still 10 absent stopovers per presence
#  9 Sep 2014: using correct country data for absences, and filling in "at sea" points for presences
#
##########################################################

cluster <- FALSE

if (cluster) {
  library(plyr)
  library(reshape)
  #library(glmmML)
  library(lme4)
  #library(AICcmodavg)
  library(arm)
}

if (!cluster) {
  library(plyr)
  library(reshape)
  #library(glmmML)
  library(lme4)
  #library(AICcmodavg)
  library(arm)
  library(ggplot2)
  library(grid)
  library(scales)
  library(sp)
  library(rgeos)
  library(rgdal)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



###----------------------------------------------###
####         SET WORKING DIRECTORIES  	####
###----------------------------------------------###


Mac <- FALSE

randomradius <- 200

if (!cluster) {
  
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/cuckoos")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/cuckoos")
  
}


if (cluster) parentwd <- c("/users1/samf/cuckoos")

datawd <- paste(parentwd, "/data", sep="")
outputwd <- paste(parentwd, "/output", sep="")


####==== IMPORT POINT DATA ====####

setwd(paste(datawd, "/data for analysis/", sep=""))
present <- read.csv("presence data all variables points all with country data.csv", header=T)
absent <- read.csv(paste("absence data all variables points all with country data ", randomradius, " km.csv", sep=""), header=T)


### these data are without at sea country data and also use the incorrect country points for absences

# setwd(paste(datawd, "/data for analysis/", sep=""))
# present2 <- read.csv("presence data all variables points.csv", header=T)
# absent2 <- read.csv(paste("absence data all variables points ", randomradius, " km.csv", sep=""), header=T)

# present <- rename(present, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))
# 
# absent <- rename(absent, c("LAND.CLASS_agriculture"="agriculture", "LAND.CLASS_forest"="forest", "LAND.CLASS_scrub.grassland"="scrub.grassland","LAND.CLASS_unsuitable"="unsuitable", "LAND.CLASS_wetland.water"="wetland.water"))



###-----------------------------------------------------###
####   HABITAT USE - plot habitat types used vs random ####
###-----------------------------------------------------###

### ----------------------------- ###
####        DATA PREPARATION     ####
### ------------------------------###

newpresent <- rename(present, c("newlongs.epsg3035"="long.epsg3035", "newlats.epsg3035"="lat.epsg3035", "newlongs.epsg4326"="long.epsg4326", "newlats.epsg4326"="lat.epsg4326"))
newabsent <- rename(absent, c("randomlong.epsg3035"="long.epsg3035", "randomlat.epsg3035"="lat.epsg3035", "randomlong.epsg4326"="long.epsg4326", "randomlat.epsg4326"="lat.epsg4326"))

#### ADD STOPOVER DURATION TO DATA ####

### creates 2 files: present.LOS and absent.LOS with length of stay data added onto present and absent datasets
setwd(paste(parentwd, "/scripts/", sep=""))
source("source code to add stopover duration data for analysis.R")

############# subset out only one round of absences
#newabsent <- subset(newabsent, nabsence==1)
newabsent <- newabsent[,-which(names(newabsent) %in% c("nabsence","random.scale"))]
absent.LOS <- absent.LOS[,-which(names(absent.LOS) %in% c("nabsence","random.scale"))]

present.LOS <- present.LOS[,-which(names(present.LOS) %in% c("datetime"))]
present.LOS <- present.LOS[c(2,1,3:length(present.LOS))]

alldata <- rbind(newpresent, newabsent)

alldata.LOS <- rbind(present.LOS, absent.LOS)

### add 3-level designation level variable

alldata.LOS$desig <- factor(alldata.LOS$desiglevel, levels=c(levels(alldata.LOS$desiglevel), "N"))
alldata.LOS$desig[which(is.na(alldata.LOS$desiglevel))] <- "N"

################# TESTING PLOTS BY LOS ###############

plotLOSdat <- alldata.LOS[alldata.LOS$presence==1,c("name","mgroup","usecountry","LOS","laststop","strategy","Sahara.success","year")]

plotLOSdat2 <- unique(plotLOSdat)
plotLOSdat2 <- droplevels(plotLOSdat2)

meanLOSbycountry <- data.frame(mean=tapply(plotLOSdat2$LOS, plotLOSdat2$usecountry, mean), sd=tapply(plotLOSdat2$LOS, plotLOSdat2$usecountry, sd))
#boxplot(LOS~usecountry, data=plotLOSdat2, las=2)
qplot(factor(usecountry), LOS, data = plotLOSdat2, geom = "boxplot")
setwd(paste(outputwd, "/autumn stopover GLMM results/", sep=""))
ggsave("mean LOS by country.jpg", width=15, height=12, units="in")


hist(plotLOSdat2$LOS, breaks=35)

# 
# plotdistdat <- alldata.LOS[alldata.LOS$presence==1, c("name","mgroup","id","distance")]
# plotdistdat <- unique(plotdistdat)
# 
# hist(plotdistdat$distance[plotdistdat$distance<50], breaks=50)

###############################################
###############################################
###############################################

#### ---- PLOTTING DATASET ---- ####

analysistype2 <- c("FM")


SW <- TRUE #### ADDED FOR MS PLOTS ####


for (i in 1:length(analysistype2)) {
  
  analysistype <- analysistype2[i]
  
  fulldat <- with(alldata.LOS, data.frame(presence, name, mgroup=as.factor(mgroup), id=as.factor(id), strat=strategy, elevation=rescale(elevation.m), spei.Mar, spei.Aug, PA=PAoverlap, habitat=LAND.CLASS, country=usecountry, desig, laststop, Sahara.success, LOS))
  
  # if full dataset model
  if (analysistype=="FM") {
    
    ################ ADDED FOR MS PLOTS ###################
    if (SW) {
      dat <- subset(fulldat, strat=="SW")
      dat <- droplevels(dat)
    }
    
    if (!SW) {
      dat <- subset(fulldat, strat=="SE")
      dat <- droplevels(dat)
    }
    
    #######################################################
    # dat <- fulldat
  }
  
  
}
  
  ### PLOT SPECIFICATIONS
  
  # plot label
  if (analysistype=="FM") analysislabel <- "OVERALL"
  if (analysistype=="DL") analysislabel <- "DESIGNATION LEVEL"
  if (analysistype=="SS") analysislabel <- "SHORT STOPOVERS"
  if (analysistype=="LS") analysislabel <- "LONG STOPOVERS"
  if (analysistype=="FS") analysislabel <- "FINAL STOPOVERS"
  if (analysistype=="DW") analysislabel <- "DESIGNATION LEVEL, LUMPING No and National for SW"
  
  plottheme <- theme(
    axis.title.x = element_text(face="bold", size=16, vjust=0.1), 
    axis.text.x = element_text(size=14, colour="black", angle=0, vjust=0.5), 
    axis.title.y = element_text(face="bold", size=16, vjust=0.9), 
    axis.text.y = element_text(size=14, vjust=0.5, colour="black"),
    plot.title = element_text(size=16, face="bold", vjust=0.5, hjust=0),
    legend.text=element_text(size=14), 
    legend.title=element_text(size=14),
    legend.key=element_rect(fill="white"),
    legend.key.width=unit(1,"cm"),
    panel.background=element_rect(colour="black", fill="white"),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.text.x=element_text(face="bold", size=14)
    )
  

  
  ###-----------------------------------------###
  ####            PROTECTED AREAS            ####
  ###-----------------------------------------###
  
  #### PROPORTION OF ALL POINTS THAT OVERLAP WITH PROTECTED AREA, BY COUNTRY - POINTS ####
  
  overlap.present <- prop.table(table(dat$PA, dat$country, as.factor(dat$presence))[,,"1"], 2)
  overlap.absent <- prop.table(table(dat$PA, dat$country, as.factor(dat$presence))[,,"0"], 2)
  
  ### sample sizes
  
  p.tot <- margin.table(table(dat$PA, dat$country, as.factor(dat$presence))[,,"1"],2)
  a.tot <- margin.table(table(dat$PA, dat$country, as.factor(dat$presence))[,,"0"],2)
  zerocountries <- names(which(p.tot==0))
  
  npresent <- data.frame(country=names(p.tot), presence="present", n=p.tot)
  npresent <- subset(npresent, !(npresent$country %in% zerocountries))
  npresent2 <- data.frame(npresent, x=seq(1.2, (nrow(npresent)+0.2), by=1), y=1.05)
  
  nabsent <- data.frame(country=names(a.tot), presence="absent", n=a.tot)
  nabsent <- subset(nabsent, !(nabsent$country %in% zerocountries))
  nabsent2 <- data.frame(nabsent, x=seq(0.8, (nrow(nabsent)-0.2), by=1), y=1.05)
  
  # summary table, absences and presences
  overlap.summary <- rbind(overlap.present, overlap.absent)
  overlap.summary <- overlap.summary[, apply(overlap.summary, 2, function(x) !any(is.na(x)))]
  
overlap.summary <- data.frame(overlap.summary, presence=rep(c("present","absent"), each=2), PA=row.names(overlap.summary))
  
overlapcountries <- subset(overlap.summary, PA=="Y")
  
  #names(habitat.summary)[-which(names(habitat.summary) %in% c("presence","habitat"))]
  
  overlapcountries2 <- reshape(overlapcountries, times=c(names(overlap.summary)[-which(names(overlap.summary) %in% c("presence","PA"))]), timevar="country", varying=list(c(names(overlap.summary)[-which(names(overlap.summary) %in% c("presence","PA"))])), direction="long")
  
  inPA.country <- data.frame(presence=overlapcountries2$presence, country=overlapcountries2$country, proportion=overlapcountries2[,4])
  
  inPA.country$country <- revalue(inPA.country$country, c("United.Kingdom"="United Kingdom"))
  
  ### PLOT SPECIFICATION


################## PLOTS FOR MS #######################

if (SW) {
  countrylabels <- labs(x="Country", y="Proportion of points falling in a PA", title=paste("A)", levels(dat$strat), randomradius, "km", sep=" "))
}

if (!SW) {
  countrylabels <- labs(x="Country", y="Proportion of points falling in a PA", title=paste("B)", levels(dat$strat), randomradius, "km", sep=" "))
}


########################################################
  
  ### PLOT
  
  p <- ggplot(inPA.country, aes(x=country, y=proportion, fill=presence)) + geom_bar(stat="identity", position="dodge")
  
if (SW) {
  p.SW <- p + countrylabels + plottheme + scale_fill_manual(values=c("#0066CC", "#003366"), breaks=c("present", "absent"), labels=c("presences","absences")) + guides(fill=guide_legend(title=NULL)) + geom_text(aes(x,y,label=n, size=1, angle=0), data=nabsent2, show_guide=FALSE) + geom_text(aes(x,y,label=n, size=1, angle=0), data=npresent2, show_guide=FALSE) + scale_y_continuous(limits=c(0,1.1))
}

if (!SW) {
  p.SE <- p + countrylabels + plottheme + scale_fill_manual(values=c("#0066CC", "#003366"), breaks=c("present", "absent"), labels=c("presences","absences")) + guides(fill=guide_legend(title=NULL)) + geom_text(aes(x,y,label=n, size=1, angle=0), data=nabsent2, show_guide=FALSE) + geom_text(aes(x,y,label=n, size=1, angle=0), data=npresent2, show_guide=FALSE) + scale_y_continuous(limits=c(0,1.1))
}

########################################################
########################################################
########################################################

setwd(paste(outputwd, "/autumn stopover GLMM results/", sep=""))
tiff(paste("MS_proportion of points per country in a PA SE-SW", randomradius, "km.tiff", sep=" "), res=150, height=2000, width=2000, units="px")
multiplot(p.SW, p.SE, cols=1)
dev.off()
  


# 
# 
# ####=========== HABITAT USE: POLYGONS ============####
# 
# if (!points) {
#   
#   ################# MULTIPLE RANDOM SCALES ###################
#   
#   meancorine.present <- aggregate(present[,c(habitatvarnames)], list(presence=present$presence), mean)
#   meancorine.present <- meancorine.present[,-1]
#   rownames(meancorine.present) <- "present"
#   
#   corinesummary.absent <- list()
#   for (i in 1:length(absentsplit)) {
#     absentsubset <- absentsplit[[i]]
#     corinesummary.absent[[i]] <- aggregate(absentsubset[,habitatvarnames], list(presence=absentsubset$presence), mean)
#   }
#   
#   names(corinesummary.absent) <- paste("absent.",names(absentsplit), sep="")
#   
#   meancorine.absent <- do.call(rbind,corinesummary.absent)
#   meancorine.absent <- meancorine.absent[,-1]
#   
#   meancorine <- rbind(meancorine.present, meancorine.absent)
#   
#   # meancorine <- reshape(corinesummary, times=c("agriculture","artificial","bare_land","forest","scrub_grassland","wetland_water"), timevar="landcover", varying=list(c("agriculture","artificial","bare_land","forest","scrub_grassland","wetland_water")), direction="long")
#   
#   ### barplot showing all random scales
#   par(mar=c(5,5,2,1))
#   
#   barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present","random 50 km", "random 100km", "random 200km", "random 500km"), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
#   
#   ### pooling across random scales ###
#   absentpooled <- apply(meancorine.absent, 2, mean)
#   meancorine <- rbind(meancorine.present, absentpooled)
#   rownames(meancorine) <- c("present", "absent")
#   
#   ### barplot showing pooled random scales
#   par(mar=c(5,5,2,1))
#   
#   barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present","random"), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
#   
#   ################# SINGLE RANDOM SCALE ###################
#   
#   meancorine.absent <- aggregate(absent[,c(habitatvarnames)], list(presence=absent$presence), mean)
#   meancorine.absent <- meancorine.absent[,-1]
#   rownames(meancorine.absent) <- "absent"
#   
#   meancorine <- rbind(meancorine.present, meancorine.absent)
#   
#   ### barplot showing single random scale
#   par(mar=c(5,5,2,1))
#   
#   barplot(as.matrix(meancorine), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("present", paste("random ", randomradius, " km", sep="")), args.legend=list(bty="n", cex=0.8), xlab="land cover type", ylab="mean proportion of \n land cover type used")
#   
#   #########################################################
#   
#   ###!!!!!!!!!!!!!!!!! CAUTIONARY NOTE !!!!!!!!!!!!!!!!!!!!!!###
#   
#   # mtext("random wetland-water not necessarily representative due to sampling methodology", cex=0.8, side=3, line=1)
#   
#   # more wetland/water than random, but random habitat extraction specifically tries to avoid stopover polygons with too much water, so possibly an artifact of random stopover sampling procedure
#   
#   # may need to re-run with a lower avoidance rate of wetland-water sites (sea only?), or increase the proportion to 0.6-0.7
#   
#   # latest revision of random stopover generation code to generate 10 pseudoabsences per absence allowed more wetland.water habitat (particularly inland marshes, peat bogs, and water courses) - 0.5 threshold for choosing absence stopovers without too much other water categories (marine and water bodies)
#   
# }
# 
# 
# ###-----------------------------------------###
# ####            PROTECTED AREAS            ####
# ###-----------------------------------------###
# 
# #### PROPORTION OF ALL STOPOVER POLYGONS THAT OVERLAP WITH PROTECTED AREA ####
# 
# if (!points) {
#   
#   overlap.present <- prop.table(table(present$overlap))
#   overlap.absent <- prop.table(table(absent$overlap))
#   overlap.summary <- rbind(overlap.present, overlap.absent)
#   
#   par(mar=c(5,5,2,1))
#   barplot(as.matrix(overlap.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("N", "Y"), cex.names=0.8, xpd=TRUE, legend=c("present", paste("random ", randomradius, " km", sep="")), args.legend=list(bty="n", cex=0.8), xlab="In a protected area?", ylab="proportion of all points \n in protected areas")
#   
#   ####==== Proportion by area of stopover that overlaps with PA overlap ====####
#   
#   meanpropPA.present <- mean(present$prop.overlap)
#   
#   summarypropPA.absent <- list()
#   for (i in 1:length(absentsplit)) {
#     absentsubset <- absentsplit[[i]]
#     summarypropPA.absent[[i]] <- mean(absentsubset$prop.overlap)
#   }
#   
#   names(summarypropPA.absent) <- paste("absent.",names(absentsplit), sep="")
#   
#   meanpropPA.absent <- do.call(rbind,summarypropPA.absent)
#   
#   meanpropPA <- rbind(meanpropPA.present, meanpropPA.absent)
#   rownames(meanpropPA) <- c("present",rownames(meanpropPA.absent))
#   
#   par(mar=c(5,5,2,1))
#   barplot(meanpropPA, beside=TRUE, space=0.1, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), names.arg=c("present","random 50 km", "random 100km", "random 200km", "random 500km"), cex.names=0.7, xpd=TRUE, ylab="mean proportion of stopover \n overlapping protected areas", ylim=c(0,0.4))
#   
#   ### pooled random stopovers ###
#   meanpropPA <- rbind(meanpropPA.present, mean(meanpropPA.absent))
#   
#   ### pooled barplot ###
#   par(mar=c(5,5,2,1))
#   barplot(meanpropPA, beside=TRUE, space=0.1, col=c("darkblue","lightblue"), names.arg=c("present","random"), cex.names=0.7, xpd=TRUE, ylab="mean proportion of stopover \n overlapping protected areas", ylim=c(0,0.4))
#   
# }
# 
# 
# ####==== Proportion of international vs national PA use ====####
# 
# #### POINTS ####
# 
# if (points) {
#   
#   x <- subset(present, PAoverlap=="Y")
#   y <- subset(absent, PAoverlap=="Y")
#   desig.P <- prop.table(table(x$desiglevel))
#   desig.A <- prop.table(table(y$desiglevel))
#   desig.summary <- rbind(desig.P, desig.A)
#   
#   par(mar=c(5,5,2,1), mfrow=c(1,1))
#   barplot(as.matrix(desig.summary), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("International", "National"), cex.names=0.8, xpd=TRUE, legend=c("present", paste("random ", randomradius, " km", sep="")), args.legend=list(bty="n", cex=0.8), xlab="Protected Area designation level", ylab="proportion of all points \n in protected areas")
#   
#   
# }
# 
# #### POLYGONS ####
# 
# if (!points) {
#   
#   #### MULTIPLE SCALES ###
#   
#   desigprop.P <- prop.table(table(present$overlap, present$desig_type))[2,]
#   
#   desigprop.A <- list()
#   for (i in 1:4) {
#     absentsubset <- absentsplit[[i]]
#     overlapsubset <- subset(absentsubset, overlap=="Y")
#     desigprop.A[[i]] <- prop.table(table(overlapsubset$desig_type))
#   }
#   
#   designation <- rbind(desigprop.P, do.call(rbind, desigprop.A))
#   rownames(designation) <- c("present","absent.50","absent.100","absent.200","absent.500")
#   
#   par(mar=c(5,5,2,8), xpd=TRUE)
#   barplot(designation, beside=TRUE, col=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), xlab="Protected area designation level", ylab="Proportion of overlapping protected \n area at each designation level")
#   legend("right", xpd=TRUE, c("present","random 50 km", "random 100km", "random 200km", "random 500km"), bty="n", cex=0.8, fill=c("darkblue","lightblue","dodgerblue","cyan","lightgreen"), inset=c(-0.3,-0.2))
#   
#   ### pooled random stopovers ###
#   designation <- rbind(desigprop.P, apply(do.call(rbind, desigprop.A), 2, mean))
#   
#   ### pooled random stopovers barplot ###
#   par(mar=c(5,5,2,5), xpd=TRUE)
#   barplot(designation, beside=TRUE, col=c("darkblue","lightblue"), xlab="Protected area designation level", ylab="Proportion of overlapping protected \n area at each designation level", ylim=c(0,0.7))
#   legend("right", xpd=TRUE, c("present","random"), bty="n", cex=0.8, fill=c("darkblue","lightblue"), inset=c(-0.15,-0.2))
#   
# }
# 
# 
# 
# ###--------------------------------###
# ####   PROTECTED AREAS & HABITAT ####
# ###--------------------------------###
# 
# if (points) {
#   
#   PAhab.present <- prop.table(table(present$PAoverlap, present$LAND.CLASS))
#   PAhab.absent <- prop.table(table(absent$PAoverlap, absent$LAND.CLASS))
#   overlap.summary <- rbind(PAhab.present, PAhab.absent)
#   
#   par(mar=c(5,5,2,1), mfrow=c(1,2))
#   barplot(as.matrix(PAhab.present), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c("Present, Outside PA","Present, Inside PA"), args.legend=list(bty="n", cex=0.8), xlab="Land cover class", ylab="proportion of all points")
#   
#   barplot(as.matrix(PAhab.absent), beside=TRUE, col=c("darkblue","lightblue"), names.arg=c("agriculture", "forest", "scrub \n grassland", "unsuitable", "wetland \n water"), cex.names=0.8, xpd=TRUE, legend=c(paste("Absent ", randomradius, " km, Outside PA", sep=""), paste("Absent ", randomradius , " km, Inside PA", sep="")), args.legend=list(bty="n", cex=0.8), xlab="Land cover class", ylab="proportion of all points")
#   
#   
# }