
closeAllConnections() 
rm(list=ls())
dev.off()

library(viridis)
library(ggplot2)
library(gridExtra)
library(boot)
library(egg)

# Set working directory
#setwd(~/...)

dat<-read.csv('mercator_data_shifted_latitude_range.csv', sep=',', header=T)
dat<-dat[!dat$speed>1000,]
dat2<-read.csv('speed_grid.csv', sep=',', header=T)

# Get map
world<-map_data(map='world')
world<-map_data("world") %>% 
  filter(region != "Antarctica")
  
dat2<-dat2[complete.cases(dat2),]

ED5<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_point(data=dat2[dat2$Lat>-55 & dat2$Lat<55,], 
                                                                               aes(x=Lon, y=Lat, col=vel))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  scale_colour_viridis(option='viridis',  name='Speed')+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
       labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text = element_text(size=18), plot.title=element_text(hjust=0))+
  theme( axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'))
ED5

ggsave("ED5.jpg",
       plot=ED5,
       device="jpg", dpi=300)

# Plot drifter data
ED6<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_raster(data=dat2[dat2$Lat>-55 & dat2$Lat<55,], 
                                                                               aes(x=Lon, y=Lat, fill=vel))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  scale_fill_viridis(option='viridis',  name='Speed')+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
       labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text = element_text(size=18), plot.title=element_text(hjust=0))+
  theme( axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'))
ED6


ggsave("ED6.jpg",
       plot=ED6,
       device="jpg", dpi=300)
