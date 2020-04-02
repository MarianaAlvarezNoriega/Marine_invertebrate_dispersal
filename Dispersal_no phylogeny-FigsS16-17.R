
closeAllConnections() 
rm(list=ls())
dev.off()

library(viridis)
library(ggplot2)
library(gridExtra)
library(boot)

# Set working directory
#setwd(~/...)

dat<-read.csv('speed_grid.csv', sep=',', header=T)

# Get map
world<-map_data(map='world')
world<-map_data("world") %>% 
  filter(region != "Antarctica")


Dat<-read.csv('Marine_invertebrate_data.csv', sep=',')

# Get posterior distributions
feeding<-read.csv('feeding_posterior_no_phylo.csv', sep=',')
egg<-read.csv('egg_posterior_no_phylo.csv', sep=',')
pld<-read.csv('PLD_posterior_no_phylo.csv', sep=',')


# Create new data frame to make predictions
new<-data.frame('Lat'=unique(dat$Lat))
new$absLatitude<-abs(new$Lat)
new$Hemisphere<-ifelse(new$Lat<0,'S','N')


new$fit<-c()
new$lwr<-c()
new$upr<-c()


#Select random estimates from the posterior distribution
rand_n=3000

rand_feeding<-base::sample(1:nrow(feeding), rand_n, replace=F)
rand_egg<-base::sample(1:nrow(egg), rand_n, replace=F)
rand_pld<-base::sample(1:nrow(pld), rand_n, replace=F)


# Create matrix to save predictions
mat<-matrix(NA, nrow=rand_n, ncol=nrow(new))

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
  

  
  # PLD weighted by developmental mode probability
  mat[i,]<-(pld_P)*prob_plank+(pld_L)*prob_leci
  
}

# Compute median and 95% CIs
for (j in 1:nrow(new)){
  new$fit[j]<-median(mat[,j], na.rm=T)
  new$lwr[j]<-quantile(mat[,j], 0.025,na.rm=T)
  new$upr[j]<-quantile(mat[,j], 0.975, na.rm=T)
}


# Merge predictions for each latitude with current speed data
j<-merge(dat, new, by='Lat', all.x=T)


# Convert speed from m/s to km/d
j$km_d<-j$vel*((3600*24)/1000)

# Calculate dispersal distance 
j$dispersal<-log10(j$km_d*10^(j$fit))


# Plot estimated dispersal distance
FigS17<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_raster(data=j[j$Lat>-55 & j$Lat<55,], 
                                                                               aes(x=Lon, y=Lat, fill=dispersal))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  scale_fill_viridis(option='viridis', breaks=c(-1,0,1,2,3), labels=c(0.1,1,10,expression(10^2), expression(10^3)), name='Dispersal')+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
       labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text=element_text(size=15))
FigS17


data2<-Dat[,c(7:13)]
data2<-data2[!(is.na(data2$DevelopmentalMode)==T &
                 is.na(data2$EggSize)==T &
                 is.na(data2$PlanktonicTime)==T) ,]


plotX<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_raster(data=j[j$Lat>-52 & j$Lat<55,], 
                                                                                aes(x=Lon, y=Lat, fill=fit))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  geom_count(data=data2, aes(x=Longitude,y=Latitude),  col='#CCCCCC', alpha=0.5, show.legend = T)+
  scale_fill_viridis(option='viridis', name='PD (days)', breaks=c(log10(5), log10(10), log10(15), log10(20),
                                                            log10(25),log10(30), log10(35)), 
                     limits=c(log10(5),log10(25)),
                     labels=c(5,10,'',20,'',30,35))+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text=element_text(size=15))

plotX



# Plot estimates of PLD for each random posterior distribution set of values (in grey lines),
# and the median and 95% CIs (in blue)

p<-ggplot()
i<-1

while(i <= rand_n){
  
  # Select predicted PLD for posterior distribution set of values i
  dataa<-data.frame('Lat'=new$Lat, 'fit'=mat[i,])
  p<-p+geom_line(data=dataa, aes(x=Lat, y=fit), col='#CCCCCC')
  i<-i+1
}

# Add median and CIs
FigS16<-p+
  geom_line(data=new,aes(x=Lat, y=fit), size=2, col='#0066FF')+
  geom_ribbon(data=new,aes(x=Lat, ymin=lwr,ymax=upr), alpha=0.3, fill='#0066FF')+
  theme_classic()+
  scale_y_continuous(breaks=c(log10(1),log10(2), log10(3), log10(4), log10(5), log10(6), log10(7), log10(8),log10(9),
                              log10(10),log10(20),log10(30),log10(40),log10(50),log10(60),log10(70),log10(80),
                              log10(90), log10(100), log10(200),log10(300),log10(400),log10(500),log10(600),log10(700),log10(800),
                              log10(900), log10(1000)),
                     labels=c(1,'','','','','','','','',10,'','','','','','','','',expression(10^2), '','','','','','','','',
                              expression(10^3)))+
  scale_x_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  ylab('Planktonic duration (days)')+
  theme(text = element_text(size=20), plot.title=element_text(hjust=1))+
  theme( axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'))+
  xlab('Latitude')
FigS16

