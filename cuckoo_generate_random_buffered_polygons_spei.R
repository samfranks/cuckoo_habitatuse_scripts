##########################################################
#
#  SOURCE CODE - GENERATE RANDOM STOPOVER POLYGONS
#
#  Samantha Franks
#  11 April 2014
#  7 May 2014 - modified to use buffered point polygons (allstopSPDFs) instead of MCP
#  2 July 2014 - modified to calculate SPEI from precip and temp data from 1980-2013
#  31 July 2014 - modified to extract corine data from points
#
##########################################################


library(shape)
library(splancs)
  
############################ SOURCE STARTS HERE - CLIMATE ############################
  
###----------------------------------------------------------------###
#-----              RANDOM POLYGON GENERATION FUNCTION          ------
###----------------------------------------------------------------###

# draw a circle with centroid coordinates that are the centroid of the MCP using library(shape), then use csr() from library(splancs) to generate a random point somewhere within that circle polygon that will represent the new centroid coordinates of the randomized polygon

# calculate the difference (in metres) between the new randomized centroid coordinates and the original stopover polygon's centroid coordinates, and elide the original polygon by the difference

# for each Polygons object in MCP@polygons, create function that will randomly shift the original stopover polygon somewhere within a x km radius of the original MCP midpoint

#### FUNCTION TO GENERATE RANDOM POLYGONS IN DISTANCE BANDS ####

# function chooses a random point within distance bands (0-50km, 50-100km, 100-200km, 200-500km) of original stopover centroid

getnewMCP.distancebands <- function(spdf) {
  
  if (randomradius[z]==50000) {
    circleMCPouter <- getellipse(rx=randomradius[z], ry=randomradius[z], mid=spdf@labpt)
    outerradius <- Polygon(circleMCPouter)
    boundary <- Polygons(list(outerradius), ID="1")
    boundarySP <- SpatialPolygons(list(boundary))
    randpts <- spsample(boundarySP, n=1, type="random")
  } else {
    circleMCPouter <- getellipse(rx=randomradius[z], ry=randomradius[z], mid=spdf@labpt)
    outerradius <- Polygon(circleMCPouter)
    circleMCPinner <- getellipse(rx=randomradius[z-1], ry=randomradius[z-1], mid=spdf@labpt)
    innerradius <- Polygon(circleMCPinner, hole=TRUE)
    boundary <- Polygons(list(outerradius, innerradius), ID="1")
    boundarySP <- SpatialPolygons(list(boundary))
    randpts <- spsample(boundarySP, n=1, type="random")
  }
  
  centdiff <- spdf@labpt-coordinates(randpts) # difference in metres between original MCP centroid and new MCP centroid
  
  origMCP <- SpatialPolygons(list(spdf))
  proj4string(origMCP) <- corine.crs
  
  # use elide() and argument shift=c() in package(maptools) to shift coordinates of polygon by c(x,y) amount (centdiff) - c(x,y) is the difference between the centroid of original MCP and random pt in 50km radius circle around centroid of original MCP
  newMCP <- elide(origMCP, shift=c(centdiff)) # this outputs a SpatialPolygons object, I need a list of Polygons objects so I can stick all newMCP polygons together into a SpatialPolygons object  
  
  return(newMCP@polygons) # returns a list of Polygons objects
  
}

#### FUNCTION TO GENERATE RANDOM POLYGONS AS DISTANCE FROM ORIGINAL CENTROID ####

# function chooses a random point within x km of original stopover polygon centroid
getnewMCP <- function(spdf) {
  
  circleMCP <- getellipse(rx=randomradius[z], ry=randomradius[z], mid=spdf@labpt)
  
  randpts <- csr(circleMCP,1)
  centdiff <- spdf@labpt-randpts # difference in metres between original MCP centroid and new MCP centroid
  
  origMCP <- SpatialPolygons(list(spdf))
  proj4string(origMCP) <- corine.crs
  
  # use elide() and argument shift=c() in package(maptools) to shift coordinates of polygon by c(x,y) amount (centdiff) - c(x,y) is the difference between the centroid of original MCP and random pt in 50km radius circle around centroid of original MCP
  newMCP <- elide(origMCP, shift=c(centdiff)) # this outputs a SpatialPolygons object, I need a list of Polygons objects so I can stick all newMCP polygons together into a SpatialPolygons object  
 
  return(newMCP@polygons) # returns a list of Polygons objects
  
}


##############################################################

# function chooses a random point within x km of original stopover polygon centroid
chooserandompoint <- function(spdf) {
  
  #spdf <- allstopSPDFs@polygons[[1]]
  #poly.centroid <- sapply(slot(spdf, "polygons"), function(x) slot(x, "labpt"))
  
  #circleMCP <- getellipse(rx=randomradius[z], ry=randomradius[z], mid=spdf@labpt)
  circleMCP <- getellipse(rx=10000, ry=10000, mid=spdf@labpt)
  
  
  randpts <- csr(circleMCP,1)
  centdiff <- spdf@labpt-randpts # difference in metres between original MCP centroid and new MCP centroid
  
  return(centdiff)
  
}

movepointby <- t(sapply(allstopSPDFs@polygons, chooserandompoint))

split.newdat <- split(newdat, newdat$mgroup)
split.randomdat <- list()
for (i in 1:length(split.newdat)) {
  split.randomdat[[i]] <- elide(split.newdat[[i]], shift=movepointby[i,])
}


test <- subset(newdat, mgroup=="3")
randomtest <- elide(newdat, shift=x[1,])

test.r <- crop(r, extend(extent(split.randomdat[[2]]), 10000))
plot(test.r)
plot(split.randomdat[[2]], pch=16, add=TRUE)
plot(split.newdat[[2]], pch=16, col="blue", add=TRUE)
       
x <- data.frame(t(sapply(allstopSPDFs@polygons, chooserandompoint)), mgroup=allstopSPDFs$mgroup)
y <- merge(x, newdat)
randomlongs <- y$newlongs + y$xc
randomlats <- y$newlats + y$yc

randomdat <- data.frame(newdat, randomlongs, randomlats)
coordinates(randomdat) <- c("randomlongs","randomlats")
proj4string(randomdat) <- corine.crs

test <- subset(newdat, mgroup=="3")
randomtest <- elide(newdat, shift=c(7546.902, -3645.131))




randomtest <- subset(randomdat, mgroup=="3")



tapply(newdat, newdat$mgroup, )

movepoints <- function(sppts) {
  y <- merge(newdat, x)
}
  
  origMCP <- spdf@polygons
  proj4string(origMCP) <- corine.crs
  
  # use elide() and argument shift=c() in package(maptools) to shift coordinates of polygon by c(x,y) amount (centdiff) - c(x,y) is the difference between the centroid of original MCP and random pt in 50km radius circle around centroid of original MCP
  newMCP <- elide(spdf, shift=c(centdiff)) # this outputs a SpatialPolygons object, I need a list of Polygons objects so I can stick all newMCP polygons together into a SpatialPolygons object  
  
  return(newMCP@polygons)) # # returns a list of Polygons objects
  
}
  

sapply(slot(allstopSPDFs, "polygons"), function(i) slot(i, "labpt"))



###-----------------------------------------------------------###
# EXTRACT CLIMATE VALUES FOR RANDOM STOPOVER POLYGON - CHECK FOR NA'S
###-----------------------------------------------------------###

buffer.metres <- c(25000,50000,75000,100000,125000,150000)

countouter <- 1 # counts how many times the buffer limit was exceeded and the number of random polygons needing to be generated

repeat { # OUTER REPEAT LOOP - regenerates a new random polygon if buffer limit is reached
  
  ### CREATE RANDOMIZED POLYGONS FOR PSEUDO-ABSENCES
  newMCPs <- sapply(allstopSPDFs@polygons,getnewMCP) # don't use lapply, because that effectively makes it a nested list, use sapply so that the returned object is a single list of Polygons objects from function
  
  randomstop <- SpatialPolygons(newMCPs)
  proj4string(randomstop) <- corine.crs
  randomstopSPDFs <- SpatialPolygonsDataFrame(randomstop, data.frame(mgroupyear.combn[[a]], row.names=mgroupid))
  
  # transform from EPSG 3035 to same CRS as climate rasters (EPSG 4326)
  randomstop4326 <- spTransform(randomstopSPDFs, CRS=epsg4326)
  
  countinner <- 1 # counts the buffer levels
  
  extractprecip <- extractwithargs(precip1980.2013, randomstop4326)
  names(extractprecip) <- randomstop4326@data$mgroup
  extracttemp <- extractwithargs(temp1980.2013, randomstop4326)
  names(extracttemp) <- randomstop4326@data$mgroup      
  
  problems <- do.call(rbind,lapply(extractprecip, function(x) {all(is.na(x))}))
  probmgroupids <- randomstop4326@data[which(problems==TRUE), "mgroup"] # returns mgroupid of stopover which produces ALL NAs at current buffer level
  
  # if all mean extracted climate values for all stopovers an individual makes have values (are not NA), then the climate data extraction process is complete, breaks the outer repeat loop
  # if some stopovers are NA, then enter repeat loop, first identifying the problem stopovers (the ones with NAs)
  
  if (all(problems==FALSE)) { # if there are no problem stopovers that extract all NAs
    
    useprecip <- extractprecip
    usetemp <- extracttemp
    break # breaks the outer loop
    
  } else {
    
    okmgroups <- randomstop4326@data[which(problems==FALSE), "mgroup"]
    useprecip <- extractprecip[which(names(extractprecip) %in% as.character(okmgroups))]
    usetemp <- extracttemp[which(names(extractprecip) %in% as.character(okmgroups))]
    
    repeat { # INNER REPEAT LOOP - cycles through buffer levels
      
      # for mgroups with NAs, subset the allstopover SPDF to the problem mgroups, increase the buffer level by 1, and re-extract the climate raster
      probstops <- subset(randomstop4326, mgroup %in% probmgroupids)      
      countinner <- countinner + 1
      
      # checks if the count has gone beyond the max buffer level, if so, breaks the inner loop which cycles through buffer levels, but does not break the outer loop and will repeat the stopover randomization by generating a new random stopover, and starting from the first buffer level
      if (countinner > length(buffer.metres)) {
        breakouter <- FALSE
        break # breaks the inner loop, but not the outer loop
      }
      
      prob.extractprecip <- extractwithargs(precip1980.2013, probstops)
      names(prob.extractprecip) <- probstops@data$mgroup
      prob.extracttemp <- extractwithargs(temp1980.2013, probstops)
      names(prob.extracttemp) <- probstops@data$mgroup
      
      problems <- do.call(rbind,lapply(prob.extractprecip, function(x) {all(is.na(x))}))
      
      # if all these problem mgroups now have values (ie. all are non-NA), then break out of the repeat loop (which will break the outer repeat as well) and append extracted data to the list of data from the ok mgroups
      if (all(problems==FALSE)) {
        useprecip <- c(useprecip,prob.extractprecip)
        usetemp <- c(usetemp,prob.extracttemp)
        breakouter <- TRUE
        break # breaks the inner loop AND the outer loop
        
      } else {
        
        # if there is still > 1 mgroup with a problem, then need to continue loop; first, check if any mgroups that were a problem are now ok (ie. are non-NA) at the new buffer level. For ones that are ok, append data to the list of data from ok mgroups (created at the top of the "else" statement), and repeat the loop with the mgroups that are still a problem
        okmgroups <- probstops@data[which(problems==FALSE), "mgroup"]
        rerunmgroups.precip <- prob.extractprecip[which(names(prob.extractprecip) %in% as.character(okmgroups))]
        rerunmgroups.temp <- prob.extracttemp[which(names(prob.extractprecip) %in% as.character(okmgroups))]
        
        useprecip <- c(useprecip, rerunmgroups.precip)
        usetemp <- c(usetemp, rerunmgroups.temp)
        
        probmgroupids <- probstops@data[which(problems==TRUE), "mgroup"]
        
      } # end inner IF statement
      
    } # end INNER REPEAT loop
    
    if (breakouter) break # break outer repeat loop    
    
  } # end outer IF statement (if some mgroups are NA)
  
  countouter <- countouter + 1
  
} # end OUTER REPEAT loop

# print(paste(newdat$name[1], " random stopovers generated ", countouter, " times", sep=""))
    
#     #check whether ALL are NaN
#     
#     check <- logical()
#     
#     for (i in 1:length(climate.values.test)) {
#       checkNAN <- lapply(lapply(climate.values.test[[i]], is.nan),all)
#       checktemp <- do.call(rbind, checkNAN)
#       check[i] <- any(checktemp)
#     }
#     
#     if (any(check)) {
#       countinner <- countinner + 1
#     } else {
#       breakouterloop <- TRUE
#       #       print(buffer.metres[countinner])
#       #       print("breakinner")
#       break # breaks inner and outer loops if all stopovers manage to capture a raster cell with non-NaN values
#     }
#     
#     if (countinner > length(buffer.metres)) {
#       breakouterloop <- FALSE
#       break # breaks inner loop and repeats the random polygon process to try and find some that will capture climate raster data at the 50km buffer level or less (ie. will capture climate data within a 50km buffer zone)
#     }
#     
#   } # CLOSE inner repeat loop (extracts climate data)
#   
#   if (breakouterloop) break
#   if (!breakouterloop) countouter <- countouter + 1
#   
# } # CLOSE outer repeat loop (randomizes stopovers)
# 
# climate.values <- lapply(climate.data, extract.climate, randomstop4326)

############################ SOURCE ENDS HERE - CLIMATE ############################

################## SOURCE CODE BASED ON CORINE ##############
  
#   ###-----------------------------------------------------------###
#   # EXTRACT CORINE VALUES FOR RANDOM STOPOVER POLYGON - CHECK FOR NA's, WATER
#   ###-----------------------------------------------------------###
#   
#   # NOTE: for MCPs close to the coast or with points over water, inevitably I will get a randomized polygon (or several) with no corine values (NAs), or with the polygon mainly over sea grid cells
#   # if any of the randomstops are comprised of > 50% NAs or values 39-44 (intertidal, coastal, sea and ocean, etc), then repeat polygon randomization function and raster value extraction until all randomstops sample predominantly from land
# 
#   
#   repeat {
#     
#     ### CREATE RANDOMIZED MCPS FOR PSEUDO-ABSENCES
#     newMCPs <- sapply(allstopSPDFs@polygons,getnewMCP) # don't use lapply, because that effectively makes it a nested list, use sapply so that the returned object is a single list of Polygons objects from function
#     
#     
#     #newMCPs2 <- SpatialPolygons(newMCPs)
#     #randomstop <- newMCPs2
#     #randomstop <- SpatialPolygonsDataFrame(newMCPs2, allstopSPDFs@data)
#     randomstop <- SpatialPolygons(newMCPs)
#     proj4string(randomstop) <- corine.crs
#     
#     # plot randomized stopover polygons with original stopover polygons
#     #   windows(12,12)
#     #   plot(MCP[i,], col='black',xlim=c(MCPbox[[i]][1]-50000, MCPbox[[i]][2]+50000), ylim=c(MCPbox[[i]][3]-50000, MCPbox[[i]][4]+50000))
#     #   plot(randomstop[i,], col='blue',add=TRUE,xlim=c(MCPbox[[i]][1]-50000, MCPbox[[i]][2]+50000), ylim=c(MCPbox[[i]][3]-50000, MCPbox[[i]][4]+50000))
#     
#     # NOT ADDING A BUFFER ANYMORE bufferMCP <- gBuffer(MCP, width=50, byid=TRUE) # buffer each MCP by 50m and create new SPDF object
#     
#     corine.values[[a]] <- extract(r, randomstop)
#     
#     # check that > 50% of corine values in each random stopover polygon are not from sea
#     check <- logical()
#     
#     for (i in 1:length(corine.values[[a]])) {
#       if ( (length(which(corine.values[[a]][[i]]==44 | corine.values[[a]][[i]]==39 | corine.values[[a]][[i]]==40 | corine.values[[a]][[i]]==41 | corine.values[[a]][[i]]==42 | corine.values[[a]][[i]]==43 | is.na(corine.values[[a]][[i]])))/length(corine.values[[a]][[i]]) < 0.5 ) ) {
#         check[i] <- TRUE
#       } else {
#         check[i] <- FALSE
#       }
#     }
#     
#     if (all(check==TRUE)) break
#     
#   }
#   
#   ############################ SOURCE ENDS HERE ############################