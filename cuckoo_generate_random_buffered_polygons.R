##########################################################
#
#  SOURCE CODE - GENERATE RANDOM STOPOVER POLYGONS
#
#  Samantha Franks
#  11 April 2014
#  7 May 2014 - modified to use buffered point polygons (allstopSPDFs) instead of MCP
#
##########################################################


library(shape)
library(splancs)
  
  ############################ SOURCE STARTS HERE - CLIMATE ############################
  
  ###----------------------------------------------------------------###
  # FUNCTION TO GENERATE RANDOM POLYGON WITHIN x km's of ORIGINAL
  ###----------------------------------------------------------------###
  
  # draw a circle with centroid coordinates that are the centroid of the MCP using library(shape), then use csr() from library(splancs) to generate a random point somewhere within that circle polygon that will represent the new centroid coordinates of the randomized polygon
  
  # calculate the difference (in metres) between the new randomized centroid coordinates and the original stopover polygon's centroid coordinates, and elide the original polygon by the difference
  
  # for each Polygons object in MCP@polygons, create function that will randomly shift the original stopover polygon somewhere within a x km radius of the original MCP midpoint
  
  getnewMCP <- function(spdf) {
    
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

###-----------------------------------------------------------###
# EXTRACT CLIMATE VALUES FOR RANDOM STOPOVER POLYGON - CHECK FOR NA'S
###-----------------------------------------------------------###

buffer.metres <- c(25000,50000,75000,100000,125000,150000)

countouter <- 1

repeat { # OUTER REPEAT LOOP - regenerates random MCPs
  
  ### CREATE RANDOMIZED MCPS FOR PSEUDO-ABSENCES
  newMCPs <- sapply(allstopSPDFs@polygons,getnewMCP) # don't use lapply, because that effectively makes it a nested list, use sapply so that the returned object is a single list of Polygons objects from function
  
  randomstop <- SpatialPolygons(newMCPs)
  proj4string(randomstop) <- corine.crs
  randomstopSPDFs <- SpatialPolygonsDataFrame(randomstop, data.frame(mgroup=mgroupid, row.names=mgroupid))
  
  # transform from EPSG 3035 to same CRS as climate rasters (EPSG 4326)
  randomstop4326 <- spTransform(randomstopSPDFs, CRS=epsg4326)
  
  countinner <- 1
  
  extracted <- lapply(climate.data, extract.climate, randomstop4326)
  extract.clim.values <- convertclimate(extracted)
  
  # if all mean extracted climate values for all stopovers an individual makes have values (are not NA), then the climate data extraction process is complete, breaks the outer repeat loop
  # if some stopovers are NA, then enter repeat loop, first identifying the problem stopovers (the ones with NAs)
  
  if (all(!is.na(extract.clim.values$precip.mm))) {
    
    climate.complete <- extract.clim.values
    break # breaks the outer loop
    
  } else {
    
    # writes non-NA mgroups to compiled dataset
    climate.complete <- subset(extract.clim.values, !is.na(extract.clim.values$precip))
    
    # identify mgroups with NAs
    problemmgroups <- subset(extract.clim.values, is.na(extract.clim.values$precip.mm), mgroup, drop=TRUE)
    problemmgroups <- droplevels(problemmgroups)
    probmgroupids <- levels(problemmgroups)
    
    repeat { # INNER REPEAT LOOP - cycles through buffer levels
      
      # for mgroups with NAs, subset the allstopover SPDF to the problem mgroups, increase the buffer level by 1, and re-extract the climate raster
      probstops <- subset(randomstop4326, mgroup %in% probmgroupids)
      
      countinner <- countinner + 1
      
      # checks if the count has gone beyond the max buffer level, if so, breaks the inner loop which cycles through buffer levels, but does not break the outer loop and will repeat the stopover randomization
      if (countinner > length(buffer.metres)) {
        breakouter <- FALSE
        break # breaks the inner loop, but not the outer loop
      }
      
      extracted <- lapply(climate.data, extract.climate, probstops)
      extract.clim.values <- convertclimate(extracted)
      
      # if all these problem mgroups now have values (ie. all are non-NA), then add the climate data for these stopovers to climate.complete and break out of the inner AND outer repeat loops
      if (all(!is.na(extract.clim.values$precip.mm))) {
        climate.complete <- rbind(climate.complete, extract.clim.values)
        breakouter <- TRUE
        break # breaks the inner loop AND the outer loop
        
      } else {
        
        # if there is still > 1 mgroup with a problem, then need to continue loop; first, check if any mgroups that were a problem are now ok (ie. are non-NA) at the new buffer level. For ones that are ok, add its data to climate.complete, and repeat the loop with the mgroups that are still a problem
        okmgroups <- subset(extract.clim.values, !is.na(extract.clim.values$precip.mm))
        climate.complete <- rbind(climate.complete, okmgroups)
        
        # identify mgroups still with NAs
        problemmgroups <- subset(extract.clim.values, is.na(extract.clim.values$precip.mm), mgroup, drop=TRUE)
        problemmgroups <- droplevels(problemmgroups)
        probmgroupids <- levels(problemmgroups)
        
      } # end inner IF statement
      
    } # end INNER REPEAT loop
    
    if (breakouter) break # break outer repeat loop    
    
  } # end outer IF statement (if some mgroups are NA)
  
  countouter <- countouter + 1
  
} # end OUTER REPEAT loop

print(paste(newdat$name[1], " random stopovers generated ", countouter, " times", sep=""))
    
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