###################################################
#
#   Check distribution of absences for degree of randomness
#
###################################################

library(shape)
library(splancs)

# function chooses a random point within x km of original stopover polygon centroid
randpts <- list()
circleMCP <- getellipse(rx=50, ry=50, mid=c(0,0))

npoints <- c(1000)

for (i in 1:npoints) {randpts[[i]] <- csr(circleMCP,1)}

testrandom <- data.frame(do.call(rbind, randpts))

coordinates(testrandom) <- c("xc","yc")

plot(testrandom)

c1 <- getellipse(rx=10, ry=10, mid=c(0,0))
c2 <- getellipse(rx=20, ry=20, mid=c(0,0))
c3 <- getellipse(rx=30, ry=30, mid=c(0,0))
c4 <- getellipse(rx=40, ry=40, mid=c(0,0))
c5 <- getellipse(rx=50, ry=50, mid=c(0,0))

r1 <- rbind(c1, c1[1, ])  # join
P1 <- Polygon(r1)
Ps1 <- Polygons(list(P1), ID = "a")

r2 <- rbind(c2, c2[1, ])  # join
P2 <- Polygon(r2)
Ps2 <- Polygons(list(P2), ID = "b")

r3 <- rbind(c3, c3[1, ])  # join
P3 <- Polygon(r3)
Ps3 <- Polygons(list(P3), ID = "c")

r4 <- rbind(c4, c4[1, ])  # join
P4 <- Polygon(r4)
Ps4 <- Polygons(list(P4), ID = "d")

r5 <- rbind(c5, c5[1, ])  # join
P5 <- Polygon(r5)
Ps5 <- Polygons(list(P5), ID = "e")

# Spatial Polygons Data Frame
SPs <- SpatialPolygons(list(Ps1, Ps2, Ps3, Ps4, Ps5))
plot(SPs)
plot(testrandom, add=TRUE)

x <- over(SPs, testrandom, returnList=TRUE)

y <- c(length(x[[1]]), length(x[[2]])-length(x[[1]]), length(x[[3]])-length(x[[2]]), length(x[[4]])-length(x[[3]]), length(x[[5]])-length(x[[4]]))
plot(y, main=paste(npoints, "points", sep=" "))
