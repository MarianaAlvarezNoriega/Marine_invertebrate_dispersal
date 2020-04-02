#setwd('~/...')

# Load the drifter data and the predictions from the Mercator model
# The predictions of the Mercator model are averaged at the same scale as 
# the resolution of the drifter data (1/4 degrees)

grid1<-readRDS('Grid_drifter_data.rds')
grid2<-readRDS('Grid_Mercator_data.rds')

# Exclude unreliable predictions
grid2[grid2>1000]<-NA

# Compute the difference between current speeds between data sets
grid3<-grid1-grid2

# Get latitude and longitude values
dat2<-read.csv('speed_grid.csv',sep=',')
lat<-as.numeric(unique(dat2$Lat))
lat<-lat[!lat<(-55)]
lat<-lat[!lat>(55)]
lat<-sort(lat, decreasing=T)
lon<-as.numeric(unique(dat2$Lon))
lon<-sort(lon, decreasing=T)

row.names(grid3)<-lat
colnames(grid3)<-lon

library(reshape)

# Convert matrix into data frame
diff<-melt(grid3)
colnames(diff)<-c('latitude', 'longitude', 'difference')


world<-map_data(map='world')
world<-map_data("world") %>% 
  filter(region != "Antarctica")


FigS3<-ggplot(diff, aes(x=longitude, y=latitude, fill=difference))+
  geom_tile()+
  scale_fill_gradient2(high='red', low='blue', midpoint=0,mid='white', na.value='grey', limits=c(-1,1),
                       name=expression(paste('m', s^-1)))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text=element_text(size=15))
FigS3





# Convert matrices to data frames
row.names(grid1)<-lat
colnames(grid1)<-lon
dat1<-melt(grid1)

row.names(grid2)<-lat
colnames(grid2)<-lon
dat2<-melt(grid2)

colnames(all)<-c('latitude', 'longitude', 'drifter', 'mercator')
all<-merge(dat1,dat2, by=c('latitude', 'longitude'))

cor.test(all$drifter, all$mercator)

# Create a correlation vector to save correlations across latitudes
corr<-c()

# Compute correlations for each latitude
for (i in 1:length(lat)){
  d<-all[all$latitude==lat[i],]
  d<-d[complete.cases(d),]
  if (nrow(d)>1){
    corr[i]<-cor.test(d$drifter, d$mercator, na.rm=T)$estimate
  }
  
  else {corr[i]<-NA}
}


FigS4<-ggplot()+
  geom_line(aes(x=lat, y=corr))+
  theme_classic()+
  theme(text=element_text(size=18))+
  xlab('Latitude')+ylab("Pearson's correlation coefficient")+
  scale_x_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  scale_y_continuous(limits=c(0,1))
FigS4


