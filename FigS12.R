
closeAllConnections() 
rm(list=ls())
dev.off()

library(viridis)
library(ggplot2)
library(gridExtra)
library(boot)

# Set working directory
#setwd('~/...')

dat<-read.csv('speed_grid.csv', sep=',', header=T)

# Get map
world<-map_data(map='world')
world<-map_data("world") %>% 
  filter(region != "Antarctica")


Dat<-read.csv('Marine_invertebrate_data.csv', sep=',')

# Get posterior distributions
feeding<-read.csv('feeding_posterior.csv', sep=',')
egg<-read.csv('egg_posterior.csv', sep=',')
pld<-read.csv('PLD_posterior.csv', sep=',')


# Create new data frame to make predictions
new<-data.frame('Lat'=unique(dat$Lat))
new$absLatitude<-abs(new$Lat)
new$Hemisphere<-ifelse(new$Lat<0,'S','N')
new<-new[new$absLatitude<=55,]


new$fit<-c()
new$lwr<-c()
new$upr<-c()

new2<-new; new3<-new; new4<-new; new5<-new


#Select random estimates from the posterior distribution
rand_n=3000

rand_feeding<-base::sample(1:nrow(feeding), rand_n, replace=F)
rand_egg<-base::sample(1:nrow(egg), rand_n, replace=F)
rand_pld<-base::sample(1:nrow(pld), rand_n, replace=F)


# Create matrices to save predictions
mat<-matrix(NA, nrow=rand_n, ncol=nrow(new))
mat2<-matrix(NA, nrow=rand_n, ncol=nrow(new)) # Developmental mode
mat3<-matrix(NA, nrow=rand_n, ncol=nrow(new)) # Egg size
mat4<-matrix(NA, nrow=rand_n, ncol=nrow(new)) # Planktonic duration
mat5<-matrix(NA, nrow=rand_n, ncol=nrow(new)) # Current speed


for (i in 1:rand_n){
  
  # Predict probability for each developmental mode
  pr_plank<-ifelse(new$Hemisphere=='N',
                   feeding$b_Intercept[rand_feeding[i]]+feeding$b_absLatitude[rand_feeding[i]]*new$absLatitude,
                   feeding$b_Intercept[rand_feeding[i]]+
                   (feeding$b_absLatitude[rand_feeding[i]]+feeding$b_absLatitude.HemisphereS[rand_feeding[i]])*new$absLatitude)
  prob_plank<-exp(pr_plank)/(exp(pr_plank)+1)
  
  prob_leci<-1-prob_plank
  
  

  # Predict egg size
  egg_P<-ifelse(new$Hemisphere=='N',
                egg$b_Intercept[rand_egg[i]]+egg$b_DevelopmentalModeP[rand_egg[i]]+
                  (egg$b_absLatitude[rand_egg[i]]+egg$b_absLatitude.DevelopmentalModeP[rand_egg[i]])*new$absLatitude,
                egg$b_Intercept[rand_egg[i]]+egg$b_DevelopmentalModeP[rand_egg[i]]+
                  (egg$b_absLatitude[rand_egg[i]]+egg$b_absLatitude.DevelopmentalModeP.HemisphereS[rand_egg[i]])*new$absLatitude)
  egg_L<-ifelse(new$Hemisphere=='N',egg$b_Intercept[rand_egg[i]]+
                  (egg$b_absLatitude[rand_egg[i]])*new$absLatitude,
                egg$b_Intercept[rand_egg[i]]+
                  (egg$b_absLatitude[rand_egg[i]]+egg$b_absLatitude.DevelopmentalModeL.HemisphereS[rand_egg[i]])*new$absLatitude)
  
  

  
  # Predict PLD 
  pld_P<-ifelse(new$Hemisphere=='N',
                (pld$b_Intercept[rand_pld[i]]+pld$b_DevelopmentalModeP[rand_pld[i]]+
                   (pld$b_absLatitude.DevelopmentalModeP[rand_pld[i]]+pld$b_absLatitude[rand_pld[i]])*new$absLatitude+
                   (pld$b_log10EggSize[rand_pld[i]]+pld$b_DevelopmentalModeP.log10EggSize[rand_pld[i]])*egg_P),
                (pld$b_Intercept[rand_pld[i]]+pld$b_DevelopmentalModeP[rand_pld[i]]+
                   (pld$b_absLatitude.DevelopmentalModeP[rand_pld[i]]+pld$b_absLatitude[rand_pld[i]]+
                      pld$b_absLatitude.HemisphereS[rand_pld[i]])*new$absLatitude+
                   (pld$b_log10EggSize[rand_pld[i]]+pld$b_DevelopmentalModeP.log10EggSize[rand_pld[i]])*egg_P))
  
  pld_L<-ifelse(new$Hemisphere=='N',
                (pld$b_Intercept[rand_pld[i]]+
                   (pld$b_absLatitude[rand_pld[i]])*new$absLatitude+
                   (pld$b_log10EggSize[rand_pld[i]])*egg_L),
                (pld$b_Intercept[rand_pld[i]]+
                   (pld$b_absLatitude[rand_pld[i]]+
                      pld$b_absLatitude.HemisphereS[rand_pld[i]])*new$absLatitude+
                   (pld$b_log10EggSize[rand_pld[i]])*egg_L))
  
  pld_P3<-ifelse(new$Hemisphere=='N',
                (pld$b_Intercept[rand_pld[i]]+pld$b_DevelopmentalModeP[rand_pld[i]]+
                   (pld$b_absLatitude.DevelopmentalModeP[rand_pld[i]]+pld$b_absLatitude[rand_pld[i]])*new$absLatitude+
                   (pld$b_log10EggSize[rand_pld[i]]+pld$b_DevelopmentalModeP.log10EggSize[rand_pld[i]])*1.1*egg_P),
                (pld$b_Intercept[rand_pld[i]]+pld$b_DevelopmentalModeP[rand_pld[i]]+
                   (pld$b_absLatitude.DevelopmentalModeP[rand_pld[i]]+pld$b_absLatitude[rand_pld[i]]+
                      pld$b_absLatitude.HemisphereS[rand_pld[i]])*new$absLatitude+
                   (pld$b_log10EggSize[rand_pld[i]]+pld$b_DevelopmentalModeP.log10EggSize[rand_pld[i]])*1.1*egg_P))
  
  pld_L3<-ifelse(new$Hemisphere=='N',
                (pld$b_Intercept[rand_pld[i]]+
                   (pld$b_absLatitude[rand_pld[i]])*new$absLatitude+
                   (pld$b_log10EggSize[rand_pld[i]])*1.1*egg_L),
                (pld$b_Intercept[rand_pld[i]]+
                   (pld$b_absLatitude[rand_pld[i]]+
                      pld$b_absLatitude.HemisphereS[rand_pld[i]])*new$absLatitude+
                   (pld$b_log10EggSize[rand_pld[i]])*1.1*egg_L))
  
  

  pld_P4<-ifelse(new$Hemisphere=='N',
                (pld$b_Intercept[rand_pld[i]]+pld$b_DevelopmentalModeP[rand_pld[i]]+
                   1.1*(pld$b_absLatitude.DevelopmentalModeP[rand_pld[i]]+pld$b_absLatitude[rand_pld[i]])*(new$absLatitude)+
                   (pld$b_log10EggSize[rand_pld[i]]+pld$b_DevelopmentalModeP.log10EggSize[rand_pld[i]])*egg_P),
                (pld$b_Intercept[rand_pld[i]]+pld$b_DevelopmentalModeP[rand_pld[i]]+
                   1.1*(pld$b_absLatitude.DevelopmentalModeP[rand_pld[i]]+pld$b_absLatitude[rand_pld[i]]+
                      pld$b_absLatitude.HemisphereS[rand_pld[i]])*(new$absLatitude)+
                   (pld$b_log10EggSize[rand_pld[i]]+pld$b_DevelopmentalModeP.log10EggSize[rand_pld[i]])*egg_P))
  
  pld_L4<-ifelse(new$Hemisphere=='N',
                (pld$b_Intercept[rand_pld[i]]+
                   1.1*(pld$b_absLatitude[rand_pld[i]])*(new$absLatitude)+
                   (pld$b_log10EggSize[rand_pld[i]])*egg_L),
                (pld$b_Intercept[rand_pld[i]]+
                   1.1*(pld$b_absLatitude[rand_pld[i]]+
                      pld$b_absLatitude.HemisphereS[rand_pld[i]])*(new$absLatitude)+
                   (pld$b_log10EggSize[rand_pld[i]])*egg_L))

  
  # PLD weighted by developmental mode probability
  mat[i,]<-(pld_P)*prob_plank+(pld_L)*prob_leci
  mat2[i,]<-(pld_P)*1.1*(prob_plank)+(pld_L)*(1-1.1*prob_plank)
  mat3[i,]<-(pld_P3)*prob_plank+(pld_L3)*prob_leci
  mat4[i,]<-(pld_P4)*prob_plank+(pld_L4)*prob_leci
  
}

# Compute median and 95% CIs
for (j in 1:nrow(new)){
  new$fit[j]<-median(mat[,j], na.rm=T)
  new$lwr[j]<-quantile(mat[,j], 0.025,na.rm=T)
  new$upr[j]<-quantile(mat[,j], 0.975, na.rm=T)
  new2$fit[j]<-median(mat2[,j], na.rm=T)
  new2$lwr[j]<-quantile(mat2[,j], 0.025,na.rm=T)
  new2$upr[j]<-quantile(mat2[,j], 0.975, na.rm=T)
  new3$fit[j]<-median(mat3[,j], na.rm=T)
  new3$lwr[j]<-quantile(mat3[,j], 0.025,na.rm=T)
  new3$upr[j]<-quantile(mat3[,j], 0.975, na.rm=T)
  new4$fit[j]<-median(mat4[,j], na.rm=T)
  new4$lwr[j]<-quantile(mat4[,j], 0.025,na.rm=T)
  new4$upr[j]<-quantile(mat4[,j], 0.975, na.rm=T)
}




# Merge predictions for each latitude with current speed data
j<-merge(dat, new, by='Lat', all.x=T)
j2<-merge(dat, new2, by='Lat', all.x=T)
j3<-merge(dat, new3, by='Lat', all.x=T)
j4<-merge(dat, new4, by='Lat', all.x=T)

# Convert speed from m/s to km/d
j$km_d<-j$vel*((3600*24)/1000)
j2$km_d<-j2$vel*((3600*24)/1000)
j3$km_d<-j3$vel*((3600*24)/1000)
j4$km_d<-j4$vel*((3600*24)/1000)

# Calculate dispersal distance 
j$dispersal<-log10(j$km_d*10^(j$fit))
j$dispersal_vel<-log10((1.1*j$km_d)*10^(j$fit))
j2$dispersal<-log10(j2$km_d*10^(j2$fit))
j3$dispersal<-log10(j3$km_d*10^(j3$fit))
j4$dispersal<-log10(j4$km_d*10^(j4$fit))

# Compare predictions 
j$dm_difference<-((10^j2$dispersal/10^j$dispersal)-1)
j$egg_difference<-((10^j3$dispersal/10^j$dispersal)-1)
j$pld_difference<-((10^j4$dispersal/10^j$dispersal)-1)
j$vel_difference<-((10^j$dispersal_vel/10^j$dispersal)-1)



# Plot estimated differences between predictions
plot1<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_raster(data=j[j$Lat>-55 & j$Lat<55,], 
                                                                               aes(x=Lon, y=Lat, fill=dm_difference))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
       labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text=element_text(size=15))+  scale_fill_gradient2(low='red', mid='white', high='blue', midpoint=0.1, name='')+#,
                                                         #  breaks=c(0.96,0.98,1,1.02,1.04))+
  ggtitle('a- Developmental mode')+
  theme(plot.title = element_text(hjust = 0.5))
plot1


plot2<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_raster(data=j[j$Lat>-55 & j$Lat<55,], 
                                                                                   aes(x=Lon, y=Lat, fill=egg_difference))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text=element_text(size=15))+
  scale_fill_gradient2(low='red', mid='white', high='blue', midpoint=0.1,name='')+ggtitle('b- Egg size')+
  theme(plot.title = element_text(hjust = 0.5))
plot2


plot3<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_raster(data=j[j$Lat>-55 & j$Lat<55,], 
                                                                                   aes(x=Lon, y=Lat, fill=pld_difference))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text=element_text(size=15))+
  scale_fill_gradient2(low='red', mid='white', high='blue', midpoint=0.1,name='')+ggtitle('c- Planktonic duration')+
  theme(plot.title = element_text(hjust = 0.5))
plot3


plot4<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_raster(data=j[j$Lat>-55 & j$Lat<55,], 
                                                                                   aes(x=Lon, y=Lat, fill=(vel_difference)))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text=element_text(size=15))+
  scale_fill_gradient2(low='red', mid='white', high='blue', midpoint=0.1,name='')+                                                  
  ggtitle('d- Current speed')+
  theme(plot.title = element_text(hjust = 0.5))
plot4

gA<-ggplotGrob(plot1)
gB<-ggplotGrob(plot2)
gC<-ggplotGrob(plot3)
gD<-ggplotGrob(plot4)

A<-gtable_cbind(gA,gB)
C<-gtable_cbind(gC,gD)
grid::grid.newpage(theme())
FigS12<-plot(gtable_rbind(A,C))
FigS12


