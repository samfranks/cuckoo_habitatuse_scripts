##########################################################
#
#  ADD DATA to cuckoo resampled points and land cover value extractions
#
#  Samantha Franks
#	3 Dec 2013
# 17 Dec 2013 - modified for new resampled points with corrected distances & mgroups
#
#
#                         SOURCE CODE
#
##########################################################

###-------------------------------------------------------------###
#   Add age, UK breeding sites movement #, SE/SW migratory strategy, and failed/successful Sahara crossing to each bird
#   
#   Do each bird individually
#   
#   All birds have only breeding or autumn stopover sites, EXCEPT for Lyster, who has a spring migration stopover that needs taking out (mgroup29)
#
###-------------------------------------------------------------###

corinewithgroupvar <- list()

### --- BB --- ###

dd <- corine.values[[1]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young","adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1 | dd$mgroup==24)] <- "Y"
breedsite[which(dd$mgroup!=1 & dd$mgroup!=24)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[1]] <- newdd

### --- Chance --- ###

dd <- corine.values[[2]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young","adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1 | dd$mgroup==28)] <- "Y"
breedsite[which(dd$mgroup!=1 & dd$mgroup!=28)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[2]] <- newdd

### Chris

dd <- corine.values[[3]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(levels=c("young","adult"))
age[which(dd$year=="2011")] <- "young"
age[which(dd$year!="2011")] <- "adult"

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1 | dd$mgroup==23 | dd$mgroup==49)] <- "Y"
breedsite[which(dd$mgroup!=1 & dd$mgroup!=23 & dd$mgroup!=49)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[3]] <- newdd

### --- Clement --- ###

dd <- corine.values[[4]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[4]] <- newdd

### --- Dart --- ###

dd <- corine.values[[5]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("N", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[5]] <- newdd

### --- David --- ###

dd <- corine.values[[6]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- rep("adult", nrow(dd))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1 | dd$mgroup==31)] <- "Y"
breedsite[which(dd$mgroup!=1 & dd$mgroup!=31)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[6]] <- newdd

### --- Derek --- ###

dd <- corine.values[[7]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[7]] <- newdd

### --- Idemili --- ###

dd <- corine.values[[8]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("NA"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("NA", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("NA", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[8]] <- newdd

### --- Indy --- ###

dd <- corine.values[[9]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[9]] <- newdd

### --- Iolo --- ###

dd <- corine.values[[10]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("N", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[10]] <- newdd

### --- John --- ###

dd <- corine.values[[11]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1 | dd$mgroup==2)] <- "Y"
breedsite[which(dd$mgroup!=1 & dd$mgroup!=2)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("N", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[11]] <- newdd

### --- Karma --- ###

dd <- corine.values[[12]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1 | dd$mgroup==2)] <- "Y"
breedsite[which(dd$mgroup!=1 & dd$mgroup!=2)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("NA", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("N", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[12]] <- newdd

### --- Kasper --- ###

dd <- corine.values[[13]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[13]] <- newdd

### --- Ken --- ###

dd <- corine.values[[14]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[14]] <- newdd

### --- Livingstone --- ###

dd <- corine.values[[15]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[15]] <- newdd

### --- Lloyd --- ###

dd <- corine.values[[16]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[16]] <- newdd

### --- Lyster --- ###

dd <- corine.values[[17]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young","adult"))

# add variable denoting whether mgroup is from the breeding site
# mgroup 32 is a spring migration movement
breedsite <- factor(levels=c("N","Y","spring"))
breedsite[which(dd$mgroup==1 | dd$mgroup==33)] <- "Y"
breedsite[which(dd$mgroup==29)] <- "spring"
breedsite[which(dd$mgroup!=1 & dd$mgroup!=33)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- factor(dd$year, labels=c("Y","N"))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[17]] <- newdd

### --- Martin --- ###

dd <- corine.values[[18]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[18]] <- newdd

### --- Mungo --- ###

dd <- corine.values[[19]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[19]] <- newdd

### --- Nelson --- ###

dd <- corine.values[[20]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("N", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[20]] <- newdd

### --- Nick --- ###

dd <- corine.values[[21]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[21]] <- newdd

### --- Patch --- ###

dd <- corine.values[[22]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[22]] <- newdd

### --- Reacher --- ###

dd <- corine.values[[23]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("N", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[23]] <- newdd

### --- Roy --- ###

dd <- corine.values[[24]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
# mgroup 3 is a clear fall migration movement
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[24]] <- newdd

### --- Ryder --- ###

dd <- corine.values[[25]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("N", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[25]] <- newdd

### --- Skinner --- ###

dd <- corine.values[[26]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[26]] <- newdd

### --- Sussex --- ###

dd <- corine.values[[27]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[27]] <- newdd

### --- Tor --- ###

dd <- corine.values[[28]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[28]] <- newdd

### --- Wallace --- ###

dd <- corine.values[[29]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("N", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[29]] <- newdd

### --- Waller --- ###

dd <- corine.values[[30]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("young"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SE", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[30]] <- newdd

### --- Whortle --- ###

dd <- corine.values[[31]]
levels(dd$name)

# add age of bird in each year of data
levels(as.factor(dd$year))
age <- factor(dd$year, labels=c("adult"))

# add variable denoting whether mgroup is from the breeding site
breedsite <- factor(levels=c("N","Y"))
breedsite[which(dd$mgroup==1)] <- "Y"
breedsite[which(dd$mgroup!=1)] <- "N"

# add variable denoting whether SE or SW migration strategy
strategy <- rep("SW", nrow(dd))

# add variable denoting successful or failed Sahara crossing in each year
Sahara.success <- rep("Y", nrow(dd))

newdd <- data.frame(dd,age,breedsite,strategy,Sahara.success)
corinewithgroupvar[[31]] <- newdd


