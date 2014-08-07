###############################################################
#
# SOURCE FILE TO COMPILE STOPOVER INFORMATION LOOK-UP DATASET
#
###############################################################

###-----------------------------------------------------------###
#         SET WORKING DIRECTORIES
###-----------------------------------------------------------###
2

###-----------------------------------------------------------###
#         EXTRACT STOPOVER INFORMATION
###-----------------------------------------------------------###

stopoverinfo <- list()

# for loop for each cuckoo file starts here
for (a in 1:31) {
  
  ###--- LOAD DATA, ADD PROJECTION INFO, TRANSFORM ---###
  
  if (original) setwd(paste(datawd,origwd,sep=""))
  if (!original) setwd(paste(datawd,resampwd,sep=""))
  
  dataset <- read.csv(list.files()[a], header=T)
  
  ### check that dataset is not an individual to exclude (Idemili & Karma); break out of current loop run and continue to next loop level
  if (dataset$name[1] == "Karma" | dataset$name[1] =="Idemili") {
    next
  }
  
  stopovers <- subset(dataset, stopoversite=="Y")
  stopovers <- droplevels(stopovers)
  
  stopovers$mgroup <- as.factor(stopovers$mgroup)
  
  if (nrow(stopovers) == 0){
    next
  } else {
    
    eachstopover <- list()
    
    for (i in 1:length(levels(stopovers$mgroup))){
      
      eachstopover[[i]] <- stopovers[stopovers$mgroup==levels(stopovers$mgroup)[i], c("name","mgroup","year", "laststop","strategy","Sahara.success")][1,]
    }
    
    stopoverinfo[[a]] <- do.call(rbind,eachstopover)
    
  }
  
}

allbirds.stopoverinfo <- do.call(rbind, stopoverinfo)

setwd(datawd)
write.csv(allbirds.stopoverinfo, "cuckoo stopover data - mgroup_year_laststop_etc.csv", row.names=FALSE)
