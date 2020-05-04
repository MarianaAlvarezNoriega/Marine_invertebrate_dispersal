closeAllConnections() 
rm(list=ls())
dev.off()

library(viridis)
library(ggplot2)
library(gridExtra)
library(boot)
library(egg)
library(dplyr)
library(SpatialEpi)

# Set working directory
#setwd(~/...)

dat<-read.csv('mercator_data_shifted_latitude_range.csv', sep=',', header=T)
dat$speed[dat$speed>1000]<-NA
dat<-dat[complete.cases(dat),]

  
Dat<-read.csv('Marine_invertebrate_data.csv', sep=',')

mat1<-dat[,c(3,1)]
mat2<-Dat[,c(9:10)]

library(SpatialEpi)
km1<-latlong2grid(mat1)
colnames(km1)<-c('km_lat', 'km_lon')
dat<-cbind(dat, km1)

km2<-latlong2grid(mat2)
colnames(km2)<-c('km_lat', 'km_lon')
Dat<-cbind(Dat, km2)
Dat$mercator_speed<-c()

for (i in 1:nrow(Dat)){
  dat$difference<-sqrt((dat$km_lat-Dat$km_lat[i])^2+(dat$km_lon-Dat$km_lon[i])^2)
  n<-which(dat$difference==min(dat$difference))
  
  if(dat$difference[n]<20){
    Dat$mercator_speed[i]<-dat$speed[n]
  }
  
  else {Dat$mercator_speed[i]<-NA}
  print(i/nrow(Dat))
}


model1<-brm(mercator_speed~s(Latitude), data=Dat, file='Speed vs Latitude_Mercator2',
            control = list(max_treedepth = 13, adapt_delta=0.99), chains=3,thin = 5, iter=30000,
            warmup=5000)

bayes_R2(model1)


## Speed annual is from the drifter data (Laurindo et al 2017)
model2<-brm(speed_annual~s(Latitude), data=dat, file='Speed vs Latitude',
           control = list(max_treedepth = 13, adapt_delta=0.99), chains=3,thin = 5, iter=30000,
           warmup=5000)

bayes_R2(model2)
