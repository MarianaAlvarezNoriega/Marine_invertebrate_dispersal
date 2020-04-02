
closeAllConnections() 
rm(list=ls())


library(viridis)
library(ggplot2)
library(gridExtra)
library(boot)


# Set working directory
#setwd(~/...)

dat<-read.csv('mercator_data_shifted_latitude_range.csv', sep=',', header=T)
dat<-dat[complete.cases(dat),]
dat<-dat[!dat$speed>1000,]
dat2<-read.csv('speed_grid.csv',sep=',')
lat<-as.numeric(unique(dat2$Lat))
lat<-lat[!lat<(-55)]
lat<-lat[!lat>(55)]
lat<-sort(lat, decreasing=T)
lon<-as.numeric(unique(dat2$Lon))
lon<-sort(lon, decreasing=T)

#M<-matrix(NA, nrow=length(lat), ncol=length(lon))
#for (i in 2:length(lat)){
#  for (j in 2:length(lon)){
#    Drif<-dat2[dat2$Lat==lat[i] & dat2$Lon==lon[j],]
#    M[i,j]=Drif$vel
#  }
#}

#colnames(M)<-lon; row.names(M)<-lat

#saveRDS(M, 'Grid_drifter_data.rds')




M2<-matrix(NA, nrow=length(lat), ncol=length(lon))


for (i in 2:length(lat)){
  for (j in 2:length(lon)){
    Mer<-dat[dat$latitude>lat[i],]
    Mer<-Mer[Mer$latitude<lat[i-1],]
    Mer<-Mer[Mer$longitude>lon[j],]
    Mer<-Mer[Mer$longitude<lon[j-1],]
    M2[i,j]=mean(Mer$speed, na.rm=T)
  }
}

saveRDS(M2, 'Grid_Mercator_data.rds')
