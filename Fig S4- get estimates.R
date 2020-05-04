#setwd('~...')
dat<-read.csv('mercator_data_shifted_latitude_range.csv', sep=',', header=T)
dat<-dat[!dat$speed>1000,]

library(SpatialEpi)


# Get matrix of latitudes and longitudes
mat<-dat[,c(3,1)]

# Get grid of latitudes and longitudes in terms of km 
km<-latlong2grid(mat)

# Add latitude and longitude distance in km to the data
dat$lat_km<-km[,1]
dat$lon_km<-km[,2]


# Create vector to save correlations
corr<-c()

# Specify maximum distance between two points in km
# (for our Figure S6 distance=10, 100 or 1000)
distance<-100

# Start sampling loop

for (s in 1:50){
  
  # Create vectors to be used later
  x1<-c(); x2<-c(); lat1<-c(); lat2<-c();
  lon1<-c(); lon2<-c(); diff<-c()
  
  # Sample random rows from the data set (1000 rows)
  samp<-sample( c(1:nrow(dat)),1000, replace=F)
  
  # Subset of data that only includes the random rows
  dat2<-dat[samp,]
  dat2<-dat2[complete.cases(dat2),]
  
  
  # Get 10 random current speeds within a specific distance of the location of the selected row
  for (i in 1:10){
      nd<-dat
      
      # Fix distances near the edges of longitudes 180 and -180 
      if (dat2$longitude[i]<(-175)){
        nd2<-nd[nd$longitude>175,]
        nd<-nd[!nd$longitude>175,]
        nd2$longitude<--180+abs(180-nd2$longitude)
        nd<-rbind(nd,nd2)
      }
       if (dat2$longitude[i]>(175)){
        nd2<-nd[nd$longitude<(-175),]
        nd<-nd[!nd$longitude<(-175),]
        nd2$longitude<-(abs(180-nd2$longitude))+180
        nd<-rbind(nd,nd2)
       }
      
      # Get latitudes and longitudes in km
      mat<-nd[,c(3,1)]
      km<-latlong2grid(mat)
      
      nd$lat_km<-km[,1]
      nd$lon_km<-km[,2]
      
      # Estimate Euclidean distances
    nd$distance<-sqrt((nd$lat_km-dat2$lat_km[i])^2+(nd$lon_km-dat2$lon_km[i])^2)
    
    # Exclude points farther away than the specified distance
    nd<-nd[nd$distance<distance,]
    
    # Exclude points in identical locations
    nd<-nd[!nd$distance==0,]
    
    # Get speed, latitude, and longitude of the pair of points 
    if (nrow(nd)>0){
      n<-sample(c(1:nrow(nd)),1)
      x1<-c(x1,dat2$speed[i])
      x2<-c(x2, nd$speed[n])
      lat1<-c(lat1, dat2$latitude[i])
      lat2<-c(lat2,nd$latitude[n])
      lon1<-c(lon1,dat2$longitude[i])
      lon2<-c(lon2,nd$longitude[n])
    }
    
    
    
  }
  
  # Get data frame with all pairs of points
  d<-data.frame(speed1=x1, speed2=x2, lat1=lat1, lat2=lat2,
                lon1=lon1,lon2=lon2)
  
  # Eliminate duplicated combinations of points
  d$loc1<-paste(d$lat1, d$lon1)
  d$loc2<-paste(d$lat2, d$lon2)
  d<-d[!duplicated(d),]
  d<-d[!d$loc1==d$loc2,]
  
  # Estimate the correlation coefficient, and save it in a vector
  corr<-c(corr, cor.test(d$speed1, d$speed2, na.rm=T)$estimate)
  print(s)
}


#write.table(corr, 'correlation_100km.csv', sep=',', row.names = F)
