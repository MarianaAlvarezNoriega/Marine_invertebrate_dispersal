
closeAllConnections() 
rm(list=ls())
dev.off()

library(viridis)
library(ggplot2)
library(gridExtra)
library(boot)
library(egg)
library(dplyr)

# Set working directory
#setwd(~/...)

dat<-read.csv('mercator_data_shifted_latitude_range.csv', sep=',', header=T)

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
new<-data.frame('latitude'=unique(dat$latitude))
new$absLatitude<-abs(new$latitude)
new$Hemisphere<-ifelse(new$latitude<0,'S','N')


new$fit<-c()
new$lwr<-c()
new$upr<-c()


#Select random estimates from the posterior distribution
rand_n=1000

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
  saveRDS(mat, 'matrix_mercator.rds')
  print(i)
}

# Compute median and 95% CIs
for (j in 1:nrow(new)){
  new$fit[j]<-median(mat[,j], na.rm=T)
  new$lwr[j]<-quantile(mat[,j], 0.025,na.rm=T)
  new$upr[j]<-quantile(mat[,j], 0.975, na.rm=T)
}


# Merge predictions for each latitude with current speed data
j<-merge(dat, new, by='latitude', all.x=T)

#write.table(j, 'mercator_predictions.csv', sep=',', row.n=F)



# Convert speed from m/s to km/d
j$km_d<-j$vel*((3600*24)/1000)

# Calculate dispersal distance 
j$dispersal<-log10(j$km_d*10^(j$fit))

j<-j[!j$speed>1000,]

#j<-read.csv('mercator_predictions.csv', sep=',')

# Plot estimated dispersal distance
Fig4a<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_point(data=j, 
                                         aes(x=longitude, y=latitude, col=dispersal))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  scale_colour_viridis(option='viridis', breaks=c(-1,0,1,2,3), labels=c(0.1,1,10,expression(10^2), expression(10^3)), name='Dispersal')+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
       labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text = element_text(size=18), plot.title=element_text(hjust=0))+
  theme( axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'))+
  ggtitle('a')+
  coord_map(ylim=c(-60,60), xlim=c(-180,180))
Fig4a




data2<-Dat[,c(7:13)]
data2<-data2[!(is.na(data2$DevelopmentalMode)==T &
                 is.na(data2$EggSize)==T &
                 is.na(data2$PlanktonicTime)==T) ,]


Fig3<-ggplot()+theme(panel.background = element_rect(fill='#FFFFFF'))+geom_point(data=j, 
                                              aes(x=longitude, y=latitude, colour=fit))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill='#000000')+theme_classic()+
  geom_count(data=data2, aes(x=Longitude,y=Latitude),  col='#CCCCCC', alpha=0.5, show.legend = T)+
  scale_colour_viridis(option='viridis', name='PD (days)', breaks=c(log10(5), log10(10), log10(15), log10(20),
                                                            log10(25),log10(30), log10(35)), 
                    
                     labels=c(5,10,15,20,25,30,35))+
  scale_size_continuous(breaks=c(1,5,10,15,20,25,30))+
  xlab('')+ylab('')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme(text=element_text(size=15))+
  coord_map(ylim=c(-60,60), xlim=c(-180,180))+
  theme(legend.direction = "vertical", legend.box = "horizontal", legend.title=element_text(size=12),
        legend.text=element_text(size=12))
Fig3

library(Cairo)
ggsave("Fig3.pdf",
       plot=Fig3,
        device=cairo_pdf)


# Plot estimates of PLD for each random posterior distribution set of values (in grey lines),
# and the median and 95% CIs (in blue)

p<-ggplot()
i<-1

while(i <= rand_n){
  
  # Select predicted PLD for posterior distribution set of values i
  dataa<-data.frame('latitude'=new$latitude, 'fit'=mat[i,])
  p<-p+geom_line(data=dataa, aes(x=latitude, y=fit), col='#CCCCCC')
  i<-i+1
}

# Add median and CIs
ED4<-p+
  geom_line(data=new,aes(x=latitude, y=fit), size=2, col='#0066FF')+
  geom_ribbon(data=new,aes(x=latitude, ymin=lwr,ymax=upr), alpha=0.3, fill='#0066FF')+
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
ED4


library(Cairo)
ggsave("ED4.jpg",
       plot=ED4,
       device="jpg", dpi=300)


# Get GAM predictions 
mat2<-j[,1:2]

library(SpatialEpi)
km<-latlong2grid(mat2)
colnames(km)<-c('km_lat', 'km_lon')

j<-cbind(j,km)

lat_vals<-seq(signif(min(j$km_lat),2),signif(max(j$km_lat),2), by=500)
lon_vals<-seq(signif(min(j$km_lon),2),signif(max(j$km_lon),2), by=500)

df<-expand.grid('latitude'=lat_vals, 'longitude'=lon_vals)
df2<-df[1,]
df2$speed<-NA
df2$distance<-NA
df2$dispersal<-NA


nd<-j

for (i in 1:nrow(df)){
  nd$distance<-sqrt((nd$km_lat-df$latitude[i])^2+(nd$km_lon-df$longitude[i])^2)
  n<-which(nd$distance==min(nd$distance))
  
  if (min(nd$distance)<100){
    df2<-rbind(df2, data.frame('latitude'=nd[n,]$latitude,
                               'longitude'=nd[n,]$longitude,
                               'speed'=nd[n,]$speed,
                               'distance'=nd[n,]$distance,
                               'dispersal'=nd[n,]$dispersal))
  }
  
  
  print(i/nrow(df))
}
df2<-df2[-1,]

#write.table(df2, 'dispersal_distance_subset_at100km.csv', sep=',', row.n=F)

library(brms)


prior<-get_prior(bf(dispersal~s(latitude)), data=df2)
mod<-brm(bf(dispersal~s(latitude)), data=df2, chains=3, control=list(adapt_delta=0.99, max_treedepth=15),
         prior=prior,iter=4000, file='dispersal_gam_500km')


df<-data.frame('latitude'=seq(-55,55,by=1))
pr<-as.data.frame(fitted(mod, df))
df$fit<-pr$Estimate
df$lwr<-pr$Q2.5
df$upr<-pr$Q97.5


Fig4b<-ggplot()+
  geom_line(data=df, aes(x=latitude, y=fit),col='red' , size=1)+
  geom_ribbon(data=df, aes(x=latitude, ymin=lwr, ymax=upr),fill='red' ,alpha=0.2)+
  
  
  scale_y_continuous(breaks=c(log10(25),
                              log10(50),log10(75),
                              log10(100),log10(125), log10(150),log10(175),log10(200),log10(225),
                              log10(250),log10(275),log10(300),log10(325),log10(350)),
                     labels=c('',50,'',100,'',150,'',200,'',250,'',300,'',350), limits=c(log10(20), log10(360)))+
  scale_x_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme_classic()+
  ylab('Dispersal (km)')+
  theme(text = element_text(size=18), plot.title=element_text(hjust=0))+
  theme( axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'))+
  ggtitle('b')+
  xlab('')
Fig4b


# Get distribution of predicted dispersal distances per degree of latitude
m<-matrix(NA, ncol=length(seq(from=-54, to=55, by=1)), nrow=length(seq(-1.188971,3.478839, l=100)))
colnames(m)<-seq(from=-54, to=55, by=1)
row.names(m)<-seq(-1.188971,3.478839, l=100)



for (i in 2:ncol(m)){
  
  nd<-j[j$latitude<as.numeric(colnames(m)[i]) & j$latitude>(as.numeric(colnames(m)[i-1])),]
  n<-nrow(nd)
  
  for (k in 2:nrow(m)) {
    nd2<-nd[nd$dispersal<as.numeric(row.names(m)[k]) & nd$dispersal>(as.numeric(row.names(m)[k-1])),]
    
    if (nrow(nd2)>0){
      m[k,i]<-nrow(nd2)/nrow(nd)
    }
    else {
      m[k,i]<-0
    }
   
  }
  print(i)
  
}

m2<-m
for (i in 1:ncol(m2)){
  if ( max(m[,i], na.rm=T)>0 ){
    m2[,i]<-m[,i]/max(m[,i], na.rm=T)
  }
  else { m2[,i]=0}
 
}

library(reshape)

M<-melt(m)
M<-M[complete.cases(M),]

M2<-melt(m2)
M2<-M2[complete.cases(M2),]
colnames(M)<-c('dispersal', 'latitude', 'proportion')
colnames(M2)<-c('dispersal', 'latitude', 'proportion')




Fig4c<-ggplot()+
  geom_tile(data=M2, aes(x=latitude, y=dispersal, fill=proportion))+
  geom_tile(data=M2[M2$proportion==0,],aes(x=latitude, y=dispersal), fill='white' )+
  
  scale_y_continuous(breaks=c(log10(.1),log10(.2), log10(.3), log10(.4), log10(.5), log10(.6), log10(.7), log10(.8),log10(.9),
                              log10(1),log10(2), log10(3), log10(4), log10(5), log10(6), log10(7), log10(8),log10(9),
                              log10(10),log10(20),log10(30),log10(40),log10(50),log10(60),log10(70),log10(80),
                              log10(90), log10(100), log10(200),log10(300),log10(400),log10(500),log10(600),log10(700),log10(800),
                              log10(900), log10(1000)),
                     labels=c(0.1,'','','','','','','','',1,'','','','','','','','',10,'','','','','','','','',
                              expression(10^2), '','','','','','','','',
                              expression(10^3)))+
  scale_x_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  theme_classic()+
  ylab('Dispersal (km)')+
  theme(text = element_text(size=18), plot.title=element_text(hjust=0))+
  theme( axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'))+
  ggtitle('c')+
  xlab('')+
  scale_fill_viridis(option='plasma', name='')
Fig4c





Fig4<-ggarrange(Fig4a+ theme(text = element_text(size=15)),Fig4b+ theme(text = element_text(size=15)), 
                   Fig4c+ theme(text = element_text(size=15)), ncol=3)


#saveRDS(Fig4, 'Fig4.rds')
ggsave("Fig4.pdf",
       plot=Fig4,
       device=cairo_pdf, dpi=300)

