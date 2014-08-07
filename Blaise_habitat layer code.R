library(raster)
library(lme4)
library(rgdal)
library(xlsReadWrite)
library(shapefiles)

rm(list=ls(all=TRUE))
setwd("//btodomain/Files/UNIXArchive/bbs")

lcmDM<-read.csv("N:/BBS habitats/land_elev_uk.csv")
    attach(lcmDM);names(lcmDM)

    pt = data.frame(x=easting,y=northing)
    coordinates(pt)=~x+y
    proj4string(pt)=CRS("+init=epsg:27700")
    translationUK<-spTransform(pt,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
    translationUK1<-as.data.frame(translationUK)

#################################################
#broadwood<-c(11)
wetland<-c(1,4,19)
heath<-c(8,9)
othersemi<-c(7,10,16)
broadwood<-c(11)

allgood<-c(1,4,19,8,9,7,10,16)
habtypes<-c("wetland","heath","othersemi","broadwood")

hab<-c(1,4,19,8,9,7,10,16)

#type<-matrix(0,length(lcm4),length(habtype1))
#for(i in 1:length(habtype1)){type[,i]<-as.vector(get(paste("lcm",habtype1[i],"try",sep="")))}
#GBhab<-rowSums(type)
GBhab2<-matrix(0,1076108,length(hab))
for(i in 1:length(hab)){
lcmhab1<-raster(paste("N:/BBS habitats/LCM2000/LCM2000GB_SUB_CLASS_PC_",hab[i],".tif",sep=""))
    GBhab1<-as.numeric(as.vector(lcmhab1))
    GBhab2[,i]<-GBhab1;    GBhab2[is.na(GBhab2)]<-0
}

GBhab<-rowSums(GBhab2)
    x1 <- xFromCell(lcmhab,1:ncell(lcmhab))
    y1 <- yFromCell(lcmhab,1:ncell(lcmhab))

lcmhab<-setValues(lcmhab1,GBhab)

countycutoff<-10
lcmDMhab<-extract(lcmhab,translationUK)
countyhab<-factor(lcmDM$county[lcmDMhab>25])

x25<-translationUK1[lcmDMhab>25,1]
y25<-translationUK1[lcmDMhab>25,2]

centrex<-tapply(x25,factor(countyhab),mean)[table(countyhab)>countycutoff]
centrey<-tapply(y25,factor(countyhab),mean)[table(countyhab)>countycutoff]

length(levels(county))
pal<-colorRampPalette(c('#f0f3ff','#0033BB'))(length(levels(county)))
plot(easting,northing,col=pal[match(county,names(centrex))],cex=0.4,pch=16,main="all semi-natural")



##############################################################################

yr<-1994
repeat{
    file.name<-paste("ebird",substr(yr,3,4),sep="")
    folder.name<-paste("data",substr(yr,3,4),sep="")
    file.loc<-paste(folder.name,"/",file.name,sep="")
    dat1<-readLines(file.loc)


siteBBS1<-sapply(1:length(dat1),function(i){substr(dat1[i],5,10)})
 
    ## MISMATCH ALERT
NsiteBBS1<-translationUK1[match(siteBBS1,lcmDM$site),2]
NsiteBBS<-NsiteBBS1[!is.na(NsiteBBS1)]
EsiteBBS<-translationUK1[match(siteBBS1,lcmDM$site),1][!is.na(NsiteBBS1)]
dat2<-dat1[!is.na(NsiteBBS1)]
siteBBS<-siteBBS1[!is.na(NsiteBBS1)]

pt = data.frame(x=EsiteBBS,y=NsiteBBS)
coordinates(pt)=~x+y
BBSsitehab<-extract(lcmhab,pt)
BBSsitehab[is.na(BBSsitehab)]<-0
habitatdata<-dat2[BBSsitehab>25]


    spec<-levels(factor(sapply(1:length(habitatdata),function(i){substr(habitatdata[i],44,45)})))
    habsites<-levels(factor(siteBBS[BBSsitehab>25]))
    BBShabsites1<-cbind(habsites,rep(yr,length(habsites)))
    
    if(yr==1994) {species<-spec; BBShabsites<-BBShabsites1}
    if(yr>1994) {species<-levels(factor(c(species,spec))); BBShabsites<-rbind(BBShabsites,BBShabsites1)}
    yr<-yr+1
    if(yr>2012) break
}


BBShabsitestab<-table(BBShabsites[,1],BBShabsites[,2])
BBShabsiteslab<-levels(factor(BBShabsites[,1]))
spno<-length(species)
################################################################################

yr<-1994
repeat{
  file.name<-paste("ebird",substr(yr,3,4),sep="")
  folder.name<-paste("data",substr(yr,3,4),sep="")
  file.loc<-paste(folder.name,"/",file.name,sep="")
  dat1<-readLines(file.loc)
    
  siteBBS1<-sapply(1:length(dat1),function(i){substr(dat1[i],5,10)})
  
  ## MISMATCH ALERT
  NsiteBBS1<-translationUK1[match(siteBBS1,lcmDM$site),2]
  NsiteBBS<-NsiteBBS1[!is.na(NsiteBBS1)]
  EsiteBBS<-translationUK1[match(siteBBS1,lcmDM$site),1][!is.na(NsiteBBS1)]
  dat2<-dat1[!is.na(NsiteBBS1)]
  siteBBS<-siteBBS1[!is.na(NsiteBBS1)]
  
  pt = data.frame(x=EsiteBBS,y=NsiteBBS)
  coordinates(pt)=~x+y
  BBSsitehab<-extract(lcmhab,pt)
  BBSsitehab[is.na(BBSsitehab)]<-0
  habitatdata<-dat2[BBSsitehab>25]
  
  squares<-c()
  for(i in 1:length(species)){
  spec1<-sapply(1:length(habitatdata),function(i){substr(habitatdata[i],44,45)})
  dat3<-(habitatdata[spec1==species[i]])
  if(length(dat3)>0){
    site1<-sapply(1:length(dat3),function(i){substr(dat3[i],5,10)})
  squares[i]<-length(table(site1))}else{squares[i]<-0}
}

if(yr==1994) {squares1<-squares}
if(yr>1994) {squares1<-cbind(squares1,squares)}
yr<-yr+1
if(yr>2012) break
}

result<-cbind(species,squares1)

###

# species is the list of all species on selected habitat with 25% cut-off
# squares1 is how many squares are they on each year within habitat
# next: probability of counting in each sqaure of target habitat
#       how that varies with N & E
#       likelihood in each county with over 10 squares

############################################################################


##### below get centre point for all squares within counties over 
#    habitat cut-off and at least 10 of these squares in county
### across whole country, not just BBS

N<-translationUK1[match(BBShabsiteslab,lcmDM$site),2]
E<-translationUK1[match(BBShabsiteslab,lcmDM$site),1]

################

# probability of counting in each sqaure of target habitat

chancecount<-array(rep(BBShabsitestab,spno),dim=c(length(BBShabsiteslab),19,spno))
yr<-1994
repeat{
  file.name<-paste("ebird",substr(yr,3,4),sep="")
  folder.name<-paste("data",substr(yr,3,4),sep="")
  file.loc<-paste(folder.name,"/",file.name,sep="")
  dat1<-readLines(file.loc)
  
  siteBBS1<-sapply(1:length(dat1),function(i){substr(dat1[i],5,10)})
  
  ## MISMATCH ALERT
  NsiteBBS1<-translationUK1[match(siteBBS1,lcmDM$site),2]
  NsiteBBS<-NsiteBBS1[!is.na(NsiteBBS1)]
  EsiteBBS<-translationUK1[match(siteBBS1,lcmDM$site),1][!is.na(NsiteBBS1)]
  dat2<-dat1[!is.na(NsiteBBS1)]
  siteBBS<-siteBBS1[!is.na(NsiteBBS1)]
  
  pt = data.frame(x=EsiteBBS,y=NsiteBBS)
  coordinates(pt)=~x+y
  BBSsitehab<-extract(lcmhab,pt)
  BBSsitehab[is.na(BBSsitehab)]<-0
  habitatdata<-dat2[BBSsitehab>25]
  
  for(i in 1:spno){
    spec1<-sapply(1:length(habitatdata),function(i){substr(habitatdata[i],44,45)})
    dat3<-(habitatdata[spec1==species[i]])
    
    if(length(dat3)>0){
      site1<-levels(factor(sapply(1:length(dat3),function(i){substr(dat3[i],5,10)})))}else{
        site1<-NA}
    
    chancecount[match(site1,BBShabsiteslab),yr-1993,i]<-chancecount[match(site1,BBShabsiteslab),yr-1993,i]+1    
   }

  yr<-yr+1
  if(yr>2012) break
}



#######################################


warninglist<-c()
predict<-c()
value<-c()
for(i in 1:spno){
  
  chanceofsquare<-matrix(0,dim(chancecount)[1],2)
  for(isite in 1:dim(chancecount)[1]){
    chanceofsquare[isite,]<-c(table(chancecount[isite,,i])[c(3,2)])
  }
  chanceofsquare[is.na(chanceofsquare)]<-0
  
  localWarnings <- list()
  value <- withCallingHandlers(glm(chanceofsquare~N+E,family=binomial), 
                               warning = function(w) {
                                 localWarnings[[length(localWarnings)+1]] <<- w
                                 invokeRestart("muffleWarning")
                               })
  warninglist<-cbind(warninglist,as.vector(list(i, localWarnings)))
  
  chance<-glm(chanceofsquare~N+E,family=binomial)
  predict[i]<-sum(1/(1+(1/exp(coef(summary(chance))[1]+
                                coef(summary(chance))[2]*centrey+
                                coef(summary(chance))[3]*centrex))))
  
}

hist(predict)
# totalpredict is the number of extra squares for each species if one
# square in each county (above cut-off) is monitored.
# The number of squares will therefore vary with habitat
# speciespredicted is divided by the number of counties above threshold

totalpredicted<-predict
for(i in 1:spno) {if(!warninglist[2,i]=="list()") totalpredicted[i]<-NA}
speciespredicted<-totalpredicted/length(table(countyhab)[table(countyhab)>countycutoff])

#}

result1<-cbind(species,rowMeans(squares1[,c(1:7,9:19)]),speciespredicted*50)

hist(speciespredicted)
############################################################################
############################################################################