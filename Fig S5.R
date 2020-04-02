#setwd('~/...')

# Drifter matrix (data from Laurindo et al. 2017)
grid1<-readRDS('Grid_drifter_data.rds')

# Mercator model prediction matrix (in 0.25 degree resolution)
grid2<-readRDS('Grid_Mercator_data.rds')

grid2[grid2>1000]<-NA

# Data frame with speed from Laurindo et al. 2017
dat2<-read.csv('speed_grid.csv',sep=',')

# Get latitude and longitude values for the matrices
lat<-as.numeric(unique(dat2$Lat))
lat<-lat[!lat<(-55)]
lat<-lat[!lat>(55)]
lat<-sort(lat, decreasing=T)
lon<-as.numeric(unique(dat2$Lon))
lon<-sort(lon, decreasing=T)

row.names(grid1)<-lat
colnames(grid1)<-lon

row.names(grid2)<-lat
colnames(grid2)<-lon

library(reshape)

# Convert matrices into data frames
drifter<-melt(grid1)
mercator<-melt(grid2)

colnames(drifter)<-c('latitude', 'longitude', 'speed')
colnames(mercator)<-c('latitude', 'longitude', 'speed')

# Merge the two data sets (from the drifter data and the Mercator model)
d<-merge(drifter, mercator, by=c('latitude', 'longitude'))


# Compute the correlation between current speeds
correlation<-cor.test(d$speed.x, d$speed.y)


FigS5<-ggplot(d, aes(x=speed.x, y=speed.y))+
  geom_point(fill='grey', pch=21, col='black', alpha=0.1)+
  theme_classic()+
  xlab(expression(paste('Current speed from the drifter data (m', s^-1, ')')))+
  ylab(expression(paste('Current speed from the Mercator model (m', s^-1, ')')))+
  theme(text=element_text(size=15))+
  geom_abline(aes(slope=1, intercept=0), col='red')
FigS5
            
