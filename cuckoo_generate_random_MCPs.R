##########################################################
#
#  SOURCE CODE - GENERATE RANDOM STOPOVER POLYGONS
#
#  Samantha Franks
#  11 April 2014
#
##########################################################


library(shape)
library(splancs)
  
  ############################ SOURCE STARTS HERE ############################
  
  ###----------------------------------------------------------------###
  # FUNCTION TO GENERATE RANDOM MCP WITHIN 50km's of ORIGINAL
  ###----------------------------------------------------------------###
  
  # draw a circle with centroid coordinates that are the centroid of the MCP using library(shape), then use csr() from library(splancs) to generate a random point somewhere within that circle polygon that will represent the new centroid coordinates of the randomized polygon
  
  # calculate the difference (in metres) between the new randomized centroid coordinates and the original stopover polygon's centroid coordinates, and elide the original polygon by the difference
  
  # for each Polygons object in MCP@polygons, create function that will randomly shift the original stopover polygon somewhere within a 50km radius of the original MCP midpoint
  
  getnewMCP <- function(spdf) {
    circleMCP <- getellipse(rx=randomradius, ry=randomradius, mid=spdf@labpt)
    
    randpts <- csr(circleMCP,1)
    centdiff <- spdf@labpt-randpts # difference in metres between original MCP centroid and new MCP centroid
    
    origMCP <- SpatialPolygons(list(spdf))
    proj4string(origMCP) <- corine.crs
    
    # use elide() and argument shift=c() in package(maptools) to shift coordinates of polygon by c(x,y) amount (centdiff) - c(x,y) is the difference between the centroid of original MCP and random pt in 50km radius circle around centroid of original MCP
    newMCP <- elide(origMCP, shift=c(centdiff)) # this outputs a SpatialPolygons object, I need a list of Polygons objects so I can stick all newMCP polygons together into a SpatialPolygons object  
    
    return(newMCP@polygons) # returns a list of Polygons objects
    
  }
  
  
  ###----------------------------------------------------------------###
  # EXTRACT CORINE VALUES FOR RANDOM STOPOVER POLYGON - CHECK FOR NA's, WATER
  ###----------------------------------------------------------------###
  
  # NOTE: for MCPs close to the coast or with points over water, inevitably I will get a randomized polygon (or several) with no corine values (NAs), or with the polygon mainly over sea grid cells
  # if any of the randomMCPs are comprised of > 50% NAs or values 39-44 (intertidal, coastal, sea and ocean, etc), then repeat polygon randomization function and raster value extraction until all randomMCPs sample predominantly from land

  
  repeat {
    
    ### CREATE RANDOMIZED MCPS FOR PSEUDO-ABSENCES
    newMCPs <- sapply(MCP@polygons,getnewMCP) # don't use lapply, because that effectively makes it a nested list, use sapply so that the returned object is a single list of Polygons objects from function
    
    newMCPs2 <- SpatialPolygons(newMCPs)
    randomMCP <- SpatialPolygonsDataFrame(newMCPs2, MCP@data)
    proj4string(randomMCP) <- corine.crs
    
    # plot randomized stopover polygons with original stopover polygons
    #   windows(12,12)
    #   plot(MCP[i,], col='black',xlim=c(MCPbox[[i]][1]-50000, MCPbox[[i]][2]+50000), ylim=c(MCPbox[[i]][3]-50000, MCPbox[[i]][4]+50000))
    #   plot(randomMCP[i,], col='blue',add=TRUE,xlim=c(MCPbox[[i]][1]-50000, MCPbox[[i]][2]+50000), ylim=c(MCPbox[[i]][3]-50000, MCPbox[[i]][4]+50000))
    
    # NOT ADDING A BUFFER ANYMORE bufferMCP <- gBuffer(MCP, width=50, byid=TRUE) # buffer each MCP by 50m and create new SPDF object
    
    corine.values[[a]] <- extract(r, randomMCP)
    
    # check that > 50% of corine values in each random stopover polygon are not from sea
    check <- logical()
    
    for (i in 1:length(corine.values[[a]])) {
      if ( (length(which(corine.values[[a]][[i]]==44 | corine.values[[a]][[i]]==39 | corine.values[[a]][[i]]==40 | corine.values[[a]][[i]]==41 | corine.values[[a]][[i]]==42 | corine.values[[a]][[i]]==43 | is.na(corine.values[[a]][[i]])))/length(corine.values[[a]][[i]]) < 0.5 ) ) {
        check[i] <- TRUE
      } else {
        check[i] <- FALSE
      }
    }
    
    if (all(check==TRUE)) break
    
  }
  
  ############################ SOURCE ENDS HERE ############################